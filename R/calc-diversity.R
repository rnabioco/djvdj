#' Calculate repertoire diversity
#'
#' @param input Single cell object or data.frame containing V(D)J data. If a
#' data.frame is provided, the cell barcodes should be stored as row names.
#' @param data_col meta.data column containing values to use for
#' calculating diversity
#' @param cluster_col meta.data column containing cluster IDs to use for
#' grouping cells when calculating diversity. If cluster_col is omitted,
#' diversity index will be calculated using all cells.
#' @param method Method to use for calculating diversity. A named list can also
#' be passed to use multiple methods. The names should specify names for the
#' output columns.
#' @param downsample Downsample clusters to the same size when calculating
#' diversity metrics
#' @param n_boots Number of bootstrap replicates
#' @param chain Chain to use for calculating diversity. Set to NULL to include
#' all chains.
#' @param chain_col meta.data column containing chains for each cell
#' @param prefix Prefix to add to new columns
#' @param return_df Return results as a data.frame. If FALSE, results will be
#' added to the input object.
#' @param sep Separator used for storing per-chain V(D)J data for each cell
#' @return Single cell object or data.frame with diversity metrics
#' @importFrom abdiv simpson
#' @importFrom broom tidy
#' @importFrom boot boot
#'
#' @examples
#' # Calculate diversity using all cells
#' res <- calc_diversity(
#'   vdj_so,
#'   data_col = "clonotype_id",
#'   method   = abdiv::simpson
#' )
#'
#' head(res@meta.data, 1)
#'
#' # Group cells based on meta.data column before calculating diversity
#' res <- calc_diversity(
#'   vdj_sce,
#'   data_col    = "clonotype_id",
#'   cluster_col = "orig.ident"
#' )
#'
#' head(res@colData, 1)
#'
#' # Add a prefix to the new columns
#' # this is useful if multiple diversity calculations are stored in the
#' # meta.data
#' res <- calc_diversity(
#'   vdj_so,
#'   data_col = "clonotype_id",
#'   prefix   = "bcr_"
#' )
#'
#' head(res@meta.data, 1)
#'
#' # Calculate multiple metrics
#' res <- calc_diversity(
#'   vdj_sce,
#'   data_col = "clonotype_id",
#'   method = list(
#'     simpson = abdiv::simpson,
#'     shannon = abdiv::shannon
#'   )
#' )
#'
#' head(res@colData, 1)
#'
#' # Return a data.frame instead of adding the results to the input object
#' res <- calc_diversity(
#'   vdj_so,
#'   data_col  = "clonotype_id",
#'   return_df = TRUE
#' )
#'
#' head(res, 1)
#'
#' @export
calc_diversity <- function(input, data_col, cluster_col = NULL,
                           method = abdiv::simpson, downsample = FALSE,
                           n_boots = 1, chain = NULL, chain_col = "chains",
                           prefix = paste0(data_col, "_"), return_df = FALSE,
                           sep = ";") {

  # Check input values
  if (length(method) > 1 && is.null(names(method))) {
    stop("Must include names if using a list of methods.")
  }

  if (length(method) == 1 && is.null(names(method))) {
    nm <- as.character(substitute(method))
    nm <- last(nm)

    method <- purrr::set_names(list(method), nm)
  }

  if (n_boots < 1) stop("n_boots must be an integer >=1.")

  # Format input data
  names(method) <- paste0(prefix, names(method))

  meta <- vdj <- .get_meta(input)

  if (!is.null(cluster_col) && is.null(meta[[cluster_col]])) {
    stop(cluster_col, " not found in object, provide a different cluster_col")
  }

  # Add column containing original data_col values to use when merging back
  # with meta.data. This is important since data_col will be modified when
  # filtering chains.
  orig_data_col <- paste0(".", data_col)

  vdj <- dplyr::mutate(vdj, !!sym(orig_data_col) := !!sym(data_col))

  # filter chains if provided
  # if all chains removed, NA will be returned
  if (!is.null(chain)) {
    vdj <- mutate_meta(
      input,
      dplyr::mutate, !!sym(orig_data_col) := !!sym(data_col)
    )

    vdj <- fetch_vdj(
      vdj,
      vdj_cols      = c(data_col, chain_col),
      clonotype_col = data_col,
      filter_cells  = TRUE,
      unnest        = FALSE,
      sep           = sep
    )

    vdj <- .filter_chains(
      vdj,
      vdj_cols   = data_col,
      chain      = chain,
      chain_col  = chain_col,
      col_names  = "{.col}",
      allow_dups = FALSE
    )
  }

  # Get parameters for bootstrapped sampling
  # specify minimum and maximum number of cells for samples
  sam      <- dplyr::filter(vdj, !is.na(!!sym(data_col)))
  sam_rng <- c(10, Inf)
  sam_sz  <- nrow(sam)

  if (!is.null(cluster_col)) {
    sam_sz <- table(sam[[cluster_col]])
    sam_sz <- sam_sz[sam_sz >= sam_rng[1]]

    if (purrr::is_empty(sam_sz)) {
      stop(
        "To calculate diversity metrics using cluster_col, at least one ",
        "cluster must have at least ", sam_rng[1], " cells. Rerun with ",
        "cluster_col set to NULL to calculate diversity metrics for the ",
        "entire population."
      )
    }

    clsts  <- names(sam_sz)
    sam_sz <- min(sam_sz)

    sam <- dplyr::filter(sam, !!sym(cluster_col) %in% clsts)
  }

  if (sam_sz > sam_rng[2]) sam_sz <- sam_rng[2]

  # Downsample clusters to equal size
  if (!is.null(cluster_col)) sam <- dplyr::group_by(sam, !!sym(cluster_col))
  if (downsample)            sam <- dplyr::slice_sample(sam, n = sam_sz)

  # Calculate diversity
  .calc_div <- function(x, met, n_bts = n_boots) {
    fn <- function(x, i) met(table(x[i]))

    bt <- boot::boot(x, fn, R = n_bts)

    broom::tidy(bt)
  }

  div <- purrr::imap_dfr(method, ~ {
    dplyr::summarize(
      sam,
      met       = .y,
      diversity = list(.calc_div(!!sym(data_col), met = .x)),
      stderr    = purrr::map_dbl(.data$diversity, pull, "std.error"),
      diversity = purrr::map_dbl(.data$diversity, pull, "statistic")
    )
  })

  div_cols <- "diversity"

  if (n_boots > 1) div_cols <- c(div_cols, "stderr")
  else             div <- dplyr::select(div, -.data$stderr)

  div <- tidyr::pivot_longer(div, all_of(div_cols))
  div <- tidyr::unite(div, "name", .data$met, .data$name)
  div <- tidyr::pivot_wider(div)

  # Format results
  vdj_cols <- c(data_col, cluster_col)

  vdj <- dplyr::mutate(vdj, !!sym(data_col) := !!sym(orig_data_col))
  vdj <- dplyr::filter(vdj, !is.na(!!sym(data_col)))
  vdj <- dplyr::distinct(vdj, !!!syms(vdj_cols))

  if (!is.null(cluster_col)) vdj <- dplyr::left_join(vdj, div, by = cluster_col)
  else                       vdj <- dplyr::bind_cols(vdj, div)

  res <- dplyr::left_join(meta, vdj, by = vdj_cols)

  if (return_df) input <- meta

  res <- .add_meta(input, meta = res)

  res
}


#' Plot repertoire diversity
#'
#' @param input Single cell object or data.frame containing V(D)J data. If a
#' data.frame is provided, the cell barcodes should be stored as row names.
#' @param data_col meta.data column containing values to use for
#' calculating diversity
#' @param cluster_col meta.data column containing cluster IDs to use for
#' grouping cells when calculating clonotype abundance
#' @param group_col meta.data column to use for grouping IDs present in
#' cluster_col
#' @param method Function to use for calculating diversity. A named list of
#' functions can be passed to plot multiple diversity metrics,
#' e.g. list(simpson = abdiv::simpson, shannon = abdiv::shannon)
#' @param downsample Downsample clusters to the same size when calculating
#' diversity metrics
#' @param n_boots Number of bootstrap replicates to use for calculating standard error
#' @param chain Chain to use for calculating diversity. Set to NULL to include
#' all chains.
#' @param chain_col meta.data column containing chains for each cell
#' @param plot_colors Character vector containing colors for plotting
#' @param plot_lvls Character vector containing levels for ordering
#' @param facet_rows The number of facet rows for final plot, use this argument
#' if a list of functions is passed to method
#' @param facet_scales This passes a scales specification to
#' ggplot2::facet_wrap, can be 'fixed', 'free', 'free_x', or 'free_y'. 'fixed'
#' will cause plot facets to share the same scales
#' @param sep Separator used for storing per-chain V(D)J data for each cell
#' @param ... Additional arguments to pass to ggplot2, e.g. color, fill, size,
#' linetype, etc.
#' @return ggplot object
#' @importFrom abdiv simpson
#'
#' @examples
#' # Plot diversity using all cells
#' plot_diversity(vdj_so)
#'
#' # Specify method to use for calculating repertoire diversity
#' plot_diversity(
#'   vdj_sce,
#'   method = abdiv::shannon
#' )
#'
#' # Plot diversity separately for each cell cluster
#' plot_diversity(
#'   vdj_so,
#'   cluster_col = "orig.ident"
#' )
#'
#' # Plot multiple diversity metrics
#' plot_diversity(
#'   vdj_sce,
#'   cluster_col = "orig.ident",
#'   method = list(simpson = abdiv::simpson, shannon = abdiv::shannon)
#' )
#'
#' # Specify colors to use for each cell cluster
#' plot_diversity(
#'   vdj_so,
#'   cluster_col = "orig.ident",
#'   plot_colors = c(avid_2 = "green", avid_1 = "purple")
#' )
#'
#' # Specify order to use for plotting cell clusters
#' plot_diversity(
#'   vdj_sce,
#'   cluster_col = "orig.ident",
#'   plot_lvls = c("avid_2", "avid_1")
#' )
#'
#' # Specify how to organize facets when plotting multiple metrics
#'
#' mets <- list(
#'   simpson = abdiv::simpson,
#'   shannon = abdiv::shannon
#' )
#'
#' plot_diversity(
#'   vdj_so,
#'   method = mets,
#'   facet_rows = 2
#' )
#'
#' @export
plot_diversity <- function(input, data_col = "clonotype_id",
                           cluster_col = NULL, group_col = NULL,
                           method = abdiv::simpson, downsample = FALSE,
                           n_boots = 1, chain = NULL, chain_col = "chains",
                           plot_colors = NULL, plot_lvls = names(plot_colors),
                           facet_rows = 1, facet_scales = "free", sep = ";",
                           ...) {

  if (length(method) > 1 && is.null(names(method))) {
    stop("Must include names if providing a list of methods.")
  }

  .chk_group_cols(cluster_col, group_col)

  if (length(method) == 1 && is.null(names(method))) {
    nm <- as.character(substitute(method))
    nm <- dplyr::last(nm)

    method        <- list(method)
    names(method) <- nm
  }

  if (!is.null(group_col)) n_boots <- 1

  # Diversity columns
  div_cols <- paste0(names(method), "_", "diversity")
  err_cols <- paste0(names(method), "_", "stderr")

  all_div_cols <- div_cols

  if (n_boots > 1) all_div_cols <- c(all_div_cols, err_cols)

  plt_cols <- c(cluster_col, group_col, all_div_cols)

  # Calculate diversity
  # remove any existing diversity columns from plt_dat
  plt_dat <- .get_meta(input)
  plt_dat <- .prepare_meta(input = NULL, plt_dat)
  plt_dat <- dplyr::select(plt_dat, -any_of(all_div_cols))

  plt_dat <- calc_diversity(
    input       = plt_dat,
    cluster_col = cluster_col,
    method      = method,
    downsample  = downsample,
    n_boots     = n_boots,
    data_col    = data_col,
    chain       = chain,
    chain_col   = chain_col,
    prefix      = "",
    sep         = sep
  )

  # Format data for plotting
  plt_dat <- dplyr::filter(plt_dat, !is.na(!!sym(data_col)))
  plt_dat <- dplyr::distinct(plt_dat, !!!syms(plt_cols))
  plt_dat <- tidyr::pivot_longer(plt_dat, all_of(all_div_cols))

  re <- "^(.+)_(diversity|stderr)$"

  plt_dat <- tidyr::extract(plt_dat, .data$name, into = c("met", "type"), re)

  plt_dat <- tidyr::pivot_wider(
    plt_dat,
    names_from  = .data$type,
    values_from = .data$value
  )

  # Set plot levels
  lvls_col <- group_col
  lvls_col <- lvls_col %||% cluster_col

  plt_dat <- .set_lvls(plt_dat, lvls_col, plot_lvls)

  include_x_labs <- !is.null(cluster_col)

  if (!include_x_labs) cluster_col <- "met"

  # Plot arguments
  gg_args <- list(df_in = plt_dat, y = "diversity", clrs = plot_colors, ...)

  # Create grouped boxplot
  if (!is.null(group_col) && !is.null(cluster_col)) {
    gg_args$alpha         <- gg_args$alpha %||% 0.5
    gg_args$outlier.color <- gg_args$outlier.color %||% NA

    more_args <- list(x = group_col, .color = group_col, .fill = group_col)

    gg_args <- append(gg_args, more_args)

    res <- purrr::lift_dl(.create_boxes)(gg_args) +
      ggplot2::geom_jitter(
        position = ggplot2::position_jitterdodge(jitter.width = 0.05)
      ) +
      theme(legend.position = "right")

    # Create facets
    if (length(method) > 1) {
      res <- res +
        ggplot2::facet_wrap(~ met, nrow = facet_rows, scales = facet_scales) +
        ggplot2::theme(axis.title.y = ggplot2::element_blank())

    } else {
      res <- res +
        ggplot2::labs(y = names(method))
    }

    return(res)
  }

  # Create bar graphs
  # only add error bars if n_boots > 1
  if (is.null(gg_args$position)) {
    gg_args$position <- ggplot2::position_identity()
  }

  more_args <- list(x = cluster_col, .fill = cluster_col)

  gg_args <- append(gg_args, more_args)

  res <- purrr::lift_dl(.create_bars)(gg_args)

  if (n_boots > 1) {
    res <- res +
      ggplot2::geom_linerange(
        aes(
          !!sym(cluster_col),
          ymin = .data$diversity - .data$stderr,
          ymax = .data$diversity + .data$stderr
        )
      )
  }

  # Create facets
  if (length(method) > 1) {
    res <- res +
      ggplot2::facet_wrap(~ met, nrow = facet_rows, scales = "free") +
      ggplot2::theme(axis.title.y = ggplot2::element_blank())

  } else {
    res <- res +
      ggplot2::labs(y = names(method))
  }

  # Set theme
  res <- res +
    ggplot2::theme(legend.position = "none")

  if (!include_x_labs) {
    res <- res +
      ggplot2::theme(
        axis.text.x  = ggplot2::element_blank(),
        axis.ticks.x = ggplot2::element_blank()
      )
  }

  res
}

#' Create rarefaction curve
#'
#' @param input Single cell object or data.frame containing V(D)J data. If a
#' data.frame is provided, the cell barcodes should be stored as row names.
#' @param data_col meta.data column containing values to use for
#' calculating diversity
#' @param cluster_col meta.data column containing cluster IDs to use for
#' grouping cells when calculating clonotype abundance
#' @param method Method to use for calculating diversity
#' @param n_boots Number of bootstrap replicates to use for calculating
#' standard error
#' @param chain Chain to use for calculating diversity. Set to NULL to include
#' all chains.
#' @param chain_col meta.data column containing chains for each cell
#' @param plot_colors Character vector containing colors for plotting
#' @param plot_lvls Character vector containing levels for ordering
#' @param facet_rows The number of facet rows for final plot
#' @param facet_scales This passes a scales specification to
#' ggplot2::facet_wrap, can be 'fixed', 'free', 'free_x', or 'free_y'. 'fixed'
#' will cause plot facets to share the same scales.
#' @param sep Separator used for storing per-chain V(D)J data for each cell
#' @param ... Additional arguments to pass to ggplot2, e.g. color, fill, size,
#' linetype, etc.
#' @return ggplot object
#' @importFrom iNEXT iNEXT
#' @importFrom stringr str_to_lower
#' @export
plot_rarefaction <- function(input, data_col, cluster_col = NULL,
                             method = "richness", n_boots = 50,
                             chain = NULL, chain_col = "chains",
                             plot_colors = NULL,
                             plot_lvls = names(plot_colors), facet_rows = 1,
                             facet_scales = "free", sep = ";", ...) {

  # Set method
  mets <- c("richness" = 0, "shannon" = 1, "invsimpson" = 2)

  if (!all(method %in% names(mets))) {
    mets <- paste0(names(mets), collapse = ", ")

    stop("method must be one of: ", mets, ".")
  }

  method <- mets[method]

  met_labs <- purrr::set_names(names(method), as.character(unname(method)))

  # Fetch data
  vdj <- .get_meta(input)

  # filter chains if provided
  # if all chains removed, NA will be returned
  if (!is.null(chain)) {
    vdj <- fetch_vdj(
      vdj,
      vdj_cols      = c(data_col, chain_col),
      clonotype_col = data_col,
      filter_cells  = TRUE,
      unnest        = FALSE,
      sep           = sep
    )

    vdj <- .filter_chains(
      vdj,
      vdj_cols   = data_col,
      chain      = chain,
      chain_col  = chain_col,
      col_names  = "{.col}",
      allow_dups = FALSE
    )
  }

  vdj <- dplyr::filter(vdj, !is.na(!!sym(data_col)))

  # Split based on cluster_col
  if (!is.null(cluster_col)) vdj <- split(vdj, vdj[[cluster_col]])
  else                       vdj <- list(vdj)

  vdj <- purrr::map(vdj, ~ as.numeric(table(.x[[data_col]])))

  # Format rarefaction data
  plt_dat <- iNEXT::iNEXT(vdj, q = unname(method), nboot = n_boots)
  plt_dat <- plt_dat$iNextEst$size_based

  plt_dat <- dplyr::mutate(
    plt_dat,
    method = dplyr::recode(Method, "Observed" = "Rarefaction"),
    method = stringr::str_to_lower(method),
    Order.q = met_labs[as.character(Order.q)]
  )

  if (!is.null(cluster_col)) {
    plt_dat <- dplyr::rename(plt_dat, !!sym(cluster_col) := Assemblage)
  }

  plt_dat <- .set_lvls(plt_dat, cluster_col, plot_lvls)
  plt_dat <- .set_lvls(plt_dat, "Order.q", names(method))
  plt_dat <- .set_lvls(plt_dat, "method", c("rarefaction", "extrapolation"))

  # Plot standard error
  res <- ggplot2::ggplot(plt_dat, ggplot2::aes(m, qD, linetype = method))

  if (n_boots > 1) {
    gg_aes <- ggplot2::aes(x = m, ymin = qD.LCL, ymax = qD.UCL)

    if (!is.null(cluster_col))   gg_aes$fill <- sym(cluster_col)
    else if (length(method) > 1) gg_aes$fill <- sym("Order.q")

    res <- res +
      ggplot2::geom_ribbon(gg_aes, alpha = 0.25, color = NA)
  }

  # Add curves
  gg_aes  <- ggplot2::aes()
  gg_args <- list(...)

  if (!is.null(cluster_col))   gg_aes$colour  <- sym(cluster_col)
  else if (length(method) > 1) gg_aes$colour  <- sym("Order.q")
  else                         gg_args$colour <- gg_args$color %||% plot_colors

  gg_args$mapping <- gg_aes

  res <- res +
    purrr::lift_dl(ggplot2::geom_line)(gg_args) +
    ggplot2::scale_linetype_manual(values = c(1, 2)) +
    djvdj_theme() +
    labs(x = "sample size")

  if (!is.null(plot_colors)) {
    res <- res +
      ggplot2::scale_fill_manual(values = plot_colors) +
      ggplot2::scale_color_manual(values = plot_colors)
  }

  if (length(method) > 1) {
    res <- res +
      ggplot2::facet_wrap(~ Order.q, nrow = facet_rows, scales = facet_scales) +
      ggplot2::labs(y = "diversity")

    if (is.null(cluster_col)) {
      res <- res +
        ggplot2::guides(colour = "none", fill = "none")
    }

  } else {
    res <- res +
      labs(y = names(method))
  }

  res
}





