#' Calculate repertoire diversity
#'
#' @param input Single cell object or data.frame containing V(D)J data. If a
#' data.frame is provided, the cell barcodes should be stored as row names.
#' @param data_col meta.data column containing values to use for
#' calculating diversity, e.g. 'clonotype_id'
#' @param cluster_col meta.data column containing cluster IDs to use for
#' grouping cells when calculating diversity. If cluster_col is omitted,
#' diversity index will be calculated using all cells.
#' @param method Method to use for calculating diversity. A named list can also
#' be passed to use multiple methods. The names should specify names for the
#' output columns.
#' @param downsample Downsample clusters to the same size when calculating
#' diversity metrics
#' @param n_boots Number of bootstrap replicates for calculating standard error,
#' if n_boots is 0 this will be skipped.
#' @param chain Chain to use for calculating diversity. To calculate diversity
#' for a single chain, the column passed to the data_col argument must contain
#' per-chain data such as CDR3 sequences. Set to NULL to include all chains.
#' @param chain_col meta.data column containing chains for each cell
#' @param prefix Prefix to add to new columns
#' @param return_df Return results as a data.frame. If FALSE, results will be
#' added to the input object.
#' @param sep Separator used for storing per-chain V(D)J data for each cell
#' @return Single cell object or data.frame with diversity metrics
#' @importFrom abdiv simpson
#' @importFrom broom tidy
#' @importFrom boot boot
#' @seealso [plot_diversity()]
#'
#' @examples
#' # Calculate diversity using all cells
#' res <- calc_diversity(
#'   vdj_so,
#'   data_col = "clonotype_id",
#'   method   = abdiv::simpson
#' )
#'
#' head(slot(res, "meta.data"), 1)
#'
#' # Group cells based on meta.data column before calculating diversity
#' res <- calc_diversity(
#'   vdj_sce,
#'   data_col    = "clonotype_id",
#'   cluster_col = "orig.ident"
#' )
#'
#' head(slot(res, "colData"), 1)
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
#' head(slot(res, "meta.data"), 1)
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
#' head(slot(res, "colData"), 1)
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
                           n_boots = 0, chain = NULL, chain_col = "chains",
                           prefix = paste0(data_col, "_"), return_df = FALSE,
                           sep = ";") {

  # Check that columns are present in object
  .check_obj_cols(
    input, data_col, cluster_col, chain = chain, chain_col = chain_col
  )

  # Check input classes
  ARG_CLASSES$method <- list(
    arg = "method", Class = "function", len_one = FALSE
  )

  .check_args(ARG_CLASSES, environment())

  # Check input values
  if (length(method) > 1 && is.null(names(method))) {
    stop("Must include names if using a list of methods.")
  }

  if (length(method) == 1 && is.null(names(method))) {
    nm <- as.character(substitute(method))
    nm <- last(nm)

    method <- purrr::set_names(list(method), nm)
  }

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
      data_cols     = c(data_col, chain_col),
      clonotype_col = data_col,
      filter_cells  = TRUE,
      unnest        = FALSE,
      sep           = sep
    )

    vdj <- .filter_chains(
      vdj,
      data_cols  = data_col,
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
  div <- purrr::imap_dfr(method, ~ {
    dplyr::summarize(
      sam,
      met       = .y,
      diversity = list(.calc_div(!!sym(data_col), met = .x, n_bts = n_boots)),
      stderr    = purrr::map_dbl(.data$diversity, pull, "std.error"),
      diversity = purrr::map_dbl(.data$diversity, pull, "statistic")
    )
  })

  div_cols <- "diversity"

  if (n_boots > 1) div_cols <- c(div_cols, "stderr")
  else             div <- dplyr::select(div, -"stderr")

  div <- tidyr::pivot_longer(div, all_of(div_cols))
  div <- tidyr::unite(div, "name", "met", "name")
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

# Calculate diversity
# n_boots should be 0 to just calculate stat w/ no reps
# n_boots must be >1 to calculate standard error
# t0 is the calculated stat from the original data
.calc_div <- function(x, met, n_bts) {
  fn <- function(x, i) {
    res <- met(table(x[i]))

    if (is.na(res)) {
      stop(
        "The diversity metric returned NA. Some metrics will return NA if a ",
        "cluster is composed entirely of a single clonotype, ",
        "check your clusters and rerun."
      )
    }

    res
  }

  bt  <- boot::boot(x, fn, R = n_bts)
  res <- broom::tidy(bt)

  if (n_bts == 0) res <- tibble::tibble(statistic = bt$t0, std.error = NA)

  res
}


#' Plot repertoire diversity
#'
#' @param input Single cell object or data.frame containing V(D)J data. If a
#' data.frame is provided, the cell barcodes should be stored as row names.
#' @param data_col meta.data column containing values to use for
#' calculating diversity, e.g. 'clonotype_id'
#' @param cluster_col meta.data column containing cluster IDs to use for
#' grouping cells when calculating clonotype abundance
#' @param group_col meta.data column to use for grouping clusters present in
#' cluster_col
#' @param method Function to use for calculating diversity, e.g.
#' abdiv::simpson. A named list of functions can be passed to plot multiple
#' diversity metrics, e.g.
#' list(simpson = abdiv::simpson, shannon = abdiv::shannon)
#' @param downsample Downsample clusters to the same size when calculating
#' diversity metrics
#' @param n_boots Number of bootstrap replicates for calculating standard error,
#' if n_boots is 0 this will be skipped.
#' @param chain Chain to use for calculating diversity. To calculate diversity
#' for a single chain, the column passed to the data_col argument must contain
#' per-chain data such as CDR3 sequences. Set to NULL to include all chains.
#' @param chain_col meta.data column containing chains for each cell
#' @param plot_colors Character vector containing colors for plotting
#' @param plot_lvls Character vector containing levels for ordering
#' @param panel_nrow The number of rows to use for arranging plot panels
#' @param panel_scales Should scales for plot panels be fixed or free. This
#' passes a scales specification to ggplot2::facet_wrap, can be 'fixed', 'free',
#' 'free_x', or 'free_y'. 'fixed' will cause panels to share the same scales.
#' @param n_label Include a label showing the number of cells plotted
#' @param label_params Named list providing additional parameters to modify
#' clonotype and n label aesthetics, e.g. list(size = 4, color = "red")
#' @param sep Separator used for storing per-chain V(D)J data for each cell
#' @param ... Additional arguments to pass to ggplot2, e.g. color, fill, size,
#' linetype, etc.
#' @return ggplot object
#' @importFrom abdiv simpson
#' @seealso [calc_diversity()], [plot_rarefaction()]
#'
#' @examples
#' # Plot diversity using all cells
#' plot_diversity(
#'   vdj_so,
#'   data_col = "clonotype_id"
#')
#'
#' # Specify method to use for calculating repertoire diversity
#' plot_diversity(
#'   vdj_sce,
#'   data_col = "clonotype_id",
#'   method   = abdiv::shannon
#' )
#'
#' # Plot diversity separately for each cell cluster
#' plot_diversity(
#'   vdj_so,
#'   data_col    = "clonotype_id",
#'   cluster_col = "orig.ident"
#' )
#'
#' # Plot multiple diversity metrics
#' mets <- list(
#'   simpson = abdiv::simpson,
#'   shannon = abdiv::shannon
#' )
#'
#' plot_diversity(
#'   vdj_sce,
#'   data_col    = "clonotype_id",
#'   cluster_col = "orig.ident",
#'   method      = mets
#' )
#'
#' # Specify colors to use for each cell cluster
#' plot_diversity(
#'   vdj_so,
#'   data_col    = "clonotype_id",
#'   cluster_col = "orig.ident",
#'   plot_colors = c(avid_2 = "green", avid_1 = "purple")
#' )
#'
#' # Specify order to use for plotting cell clusters
#' plot_diversity(
#'   vdj_sce,
#'   data_col    = "clonotype_id",
#'   cluster_col = "orig.ident",
#'   plot_lvls   = c("avid_2", "avid_1")
#' )
#'
#' # Specify how to organize panels when plotting multiple metrics
#' plot_diversity(
#'   vdj_so,
#'   data_col   = "clonotype_id",
#'   method     = mets,
#'   panel_nrow = 2
#' )
#'
#' @export
plot_diversity <- function(input, data_col, cluster_col = NULL,
                           group_col = NULL, method = abdiv::simpson,
                           downsample = FALSE, n_boots = 0, chain = NULL,
                           chain_col = "chains", plot_colors = NULL,
                           plot_lvls = names(plot_colors), panel_nrow = NULL,
                           panel_scales = "free", n_label = TRUE,
                           label_params = list(), sep = ";", ...) {

  # Check that columns are present in object
  .check_obj_cols(
    input,
    data_col, cluster_col, group_col, chain = chain, chain_col = chain_col
  )

  # Check input classes
  ARG_CLASSES$method <- list(
    arg = "method", Class = "function", len_one = FALSE
  )

  .check_args(ARG_CLASSES, environment())

  # Check input values
  if (length(method) > 1 && is.null(names(method))) {
    stop("Must include names if providing a list of methods.")
  }

  .chk_group_cols(cluster_col, group_col, input)

  if (length(method) == 1 && is.null(names(method))) {
    nm <- as.character(substitute(method))
    nm <- dplyr::last(nm)

    method <- list(method)
    names(method) <- nm
  }

  if (!is.null(group_col)) n_boots <- 0

  # Diversity columns
  div_cols <- paste0(names(method), "_", "diversity")
  err_cols <- paste0(names(method), "_", "stderr")

  all_div_cols <- div_cols

  if (n_boots > 1) all_div_cols <- c(all_div_cols, err_cols)

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

  # Calculate number of cells for label
  .n <- nrow(plt_dat)

  keep_cols <- .get_matching_clmns(plt_dat, c(data_col, cluster_col))
  keep_cols <- c(cluster_col, data_col, keep_cols)
  plt_dat   <- dplyr::distinct(plt_dat, !!!syms(keep_cols))

  plt_dat <- tidyr::pivot_longer(plt_dat, all_of(all_div_cols))

  re <- "^(.+)_(diversity|stderr)$"

  plt_dat <- tidyr::extract(plt_dat, .data$name, into = c("met", "type"), re)

  plt_dat <- tidyr::pivot_wider(
    plt_dat,
    names_from  = "type", values_from = "value"
  )

  # Set plot levels
  lvls_col <- group_col %||% cluster_col

  plt_dat <- .set_lvls(plt_dat, lvls_col, plot_lvls)
  plt_dat <- .set_lvls(plt_dat, "met", names(method))

  include_x_labs <- !is.null(cluster_col)

  if (!include_x_labs) cluster_col <- "met"

  # Plot arguments
  gg_args <- list(df_in = plt_dat, y = "diversity", clrs = plot_colors, ...)

  if (n_label) gg_args$y_exp <- c(0.05, 0.1)

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
        ggplot2::facet_wrap(~ met, nrow = panel_nrow, scales = panel_scales) +
        ggplot2::theme(axis.title.y = ggplot2::element_blank())

    } else {
      res <- res +
        ggplot2::labs(y = names(method))
    }

  # Create bar graphs
  # only add error bars if n_boots > 1
  } else {
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
        ggplot2::facet_wrap(~ met, nrow = panel_nrow, scales = "free") +
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
  }

  if (n_label) res <- .add_n_label(res, label_params, .n)

  res
}

#' Plot rarefaction curves
#'
#' This method involves calculating species diversity for different sized
#' samples generated by randomly downsampling each cluster. By default the
#' bootstrapped 95% confidence interval will also be plotted.
#'
#' @param input Single cell object or data.frame containing V(D)J data. If a
#' data.frame is provided, the cell barcodes should be stored as row names.
#' @param data_col meta.data column containing values to use for
#' calculating diversity
#' @param cluster_col meta.data column containing cluster IDs to use for
#' grouping cells when calculating clonotype abundance
#' @param method Method to use for calculating diversity, available methods
#' include:
#'
#' - 'richness', species richness, this is equivalent to the calculation
#'   performed by [abdiv::richness()]
#' - 'shannon', the exponential of Shannon entropy
#' - 'invsimpson', the inverse Simpson index, this is equivalent to the
#'   calculation performed by [abdiv::invsimpson()]
#'
#' @param n_boots Number of bootstrap replicates for calculating standard error,
#' if n_boots is 0 this will be skipped.
#' @param chain Chain to use for calculating diversity. To calculate diversity
#' for a single chain, the column passed to the data_col argument must contain
#' per-chain data such as CDR3 sequences. Set to NULL to include all chains.
#' @param chain_col meta.data column containing chains for each cell
#' @param plot_colors Character vector containing colors for plotting
#' @param plot_lvls Character vector containing levels for ordering
#' @param panel_nrow The number of rows to use for arranging plot panels
#' @param panel_scales Should scales for plot panels be fixed or free. This
#' passes a scales specification to ggplot2::facet_wrap, can be 'fixed', 'free',
#' 'free_x', or 'free_y'. 'fixed' will cause panels to share the same scales.
#' @param ci_alpha Transparency to use when plotting 95% confidence interval
#' @param sep Separator used for storing per-chain V(D)J data for each cell
#' @param ... Additional arguments to pass to ggplot2, e.g. color, fill,
#' linetype, etc.
#' @return ggplot object
#' @importFrom iNEXT iNEXT
#' @importFrom stringr str_to_lower
#' @seealso [calc_diversity()], [plot_diversity()]
#'
#' @examples
#' # Plot rarefaction curve using all cells
#' plot_rarefaction(
#'   vdj_so,
#'   data_col = "clonotype_id"
#')
#'
#' # Specify method to use for calculating repertoire diversity
#' plot_rarefaction(
#'   vdj_sce,
#'   data_col = "clonotype_id",
#'   method   = "shannon"
#' )
#'
#' # Plot separate curves for each cell cluster
#' plot_rarefaction(
#'   vdj_so,
#'   data_col    = "clonotype_id",
#'   cluster_col = "orig.ident"
#' )
#'
#' # Specify colors to use for each cell cluster
#' plot_rarefaction(
#'   vdj_so,
#'   data_col    = "clonotype_id",
#'   cluster_col = "orig.ident",
#'   plot_colors = c(avid_1 = "green", avid_2 = "purple")
#' )
#'
#' # Specify order to use for plotting cell clusters
#' plot_rarefaction(
#'   vdj_sce,
#'   data_col    = "clonotype_id",
#'   cluster_col = "orig.ident",
#'   plot_lvls   = c("avid_2", "avid_1")
#' )
#'
#' @export
plot_rarefaction <- function(input, data_col, cluster_col = NULL,
                             method = "richness", n_boots = 50,
                             chain = NULL, chain_col = "chains",
                             plot_colors = NULL,
                             plot_lvls = names(plot_colors), panel_nrow = NULL,
                             panel_scales = "free", ci_alpha = 0.15, sep = ";",
                             ...) {

  # Check that columns are present in object
  .check_obj_cols(
    input,
    data_col, cluster_col, chain = chain, chain_col = chain_col
  )

  # Check input classes
  ARG_CLASSES$method <- list(arg = "method", len_one = FALSE)

  .check_args(ARG_CLASSES, environment())

  # Check input values
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
      data_cols     = c(data_col, chain_col),
      clonotype_col = data_col,
      filter_cells  = TRUE,
      unnest        = FALSE,
      sep           = sep
    )

    vdj <- .filter_chains(
      vdj,
      data_cols  = data_col,
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
    method  = dplyr::recode(.data$Method, "Observed" = "Rarefaction"),
    method  = stringr::str_to_lower(method),
    Order.q = met_labs[as.character(.data$Order.q)]
  )

  if (!is.null(cluster_col)) {
    plt_dat <- dplyr::rename(plt_dat, !!sym(cluster_col) := "Assemblage")
  }

  plt_dat <- .set_lvls(plt_dat, cluster_col, plot_lvls)
  plt_dat <- .set_lvls(plt_dat, "Order.q", names(method))
  plt_dat <- .set_lvls(plt_dat, "method", c("rarefaction", "extrapolation"))

  # Plot standard error
  res <- ggplot2::ggplot(
    plt_dat,
    ggplot2::aes(.data$m, .data$qD, linetype = method)
  ) +
    ggplot2::guides(linetype = ggplot2::guide_legend(title = NULL))

  if (n_boots > 1) {
    gg_aes <- ggplot2::aes(
      x = .data$m, ymin = .data$qD.LCL, ymax = .data$qD.UCL
    )

    if (!is.null(cluster_col))   gg_aes$fill <- sym(cluster_col)
    else if (length(method) > 1) gg_aes$fill <- sym("Order.q")

    res <- res +
      ggplot2::geom_ribbon(gg_aes, alpha = ci_alpha, color = NA)
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
      ggplot2::facet_wrap(~ Order.q, nrow = panel_nrow, scales = panel_scales) +
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

