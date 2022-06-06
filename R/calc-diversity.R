#' Calculate repertoire diversity
#'
#' @param input Single cell object or data.frame containing V(D)J data. If a
#' data.frame is provided, the cell barcodes should be stored as row names.
#' @param data_col meta.data column containing clonotype IDs to use for
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
#' @param chain Chain to use for calculating gene usage. Set to NULL to include
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
#'   method = abdiv::simpson
#' )
#'
#' head(res@meta.data, 1)
#'
#' # Group cells based on meta.data column before calculating diversity
#' res <- calc_diversity(
#'   vdj_sce,
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
#'   prefix = "bcr_"
#' )
#'
#' head(res@meta.data, 1)
#'
#' # Calculate multiple metrics
#' res <- calc_diversity(
#'   vdj_sce,
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
#'   return_df = TRUE
#' )
#'
#' head(res, 1)
#'
#' @export
calc_diversity <- function(input, data_col, cluster_col = NULL,
                           method = abdiv::simpson, downsample = FALSE,
                           n_boots = 100, chain = NULL, chain_col = "chains",
                           prefix = "", return_df = FALSE, sep = ";") {

  if (length(method) > 1 && is.null(names(method))) {
    stop("Must include names if using a list of methods.")
  }

  if (length(method) == 1 && is.null(names(method))) {
    nm <- as.character(substitute(method))
    nm <- last(nm)

    method <- purrr::set_names(list(method), nm)
  }

  # Format input data
  div_cols <- names(method) <- paste0(prefix, names(method))
  div_cols <- paste0(div_cols, "_", c("diversity", "stderr"))

  vdj_cols <- c(data_col, cluster_col)

  # filter chains if provided
  # if all chains removed, NA will be returned
  if (is.null(chain)) {
    meta <- .get_meta(input)
    vdj  <- dplyr::filter(meta, !is.na(!!sym(data_col)))

  } else {
    vdj <- fetch_vdj(
      input,
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
      allow_dups = FALSE,
      sep        = sep
    )

    vdj <- dplyr::filter(vdj, !is.na(!!sym(data_col)))
  }

  if (!is.null(cluster_col) && is.null(meta[[cluster_col]])) {
    stop(cluster_col, " not found in object, provide a different cluster_col")
  }

  # Get parameters for bootstrapped sampling
  sam     <- vdj
  sam_rng <- c(20, Inf)
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

  if (sam_sz > sam_rng[2]) {
    sam_sz <- sam_rng[2]
  }

  # Downsample clusters to equal size
  if (!is.null(cluster_col)) {
    sam <- dplyr::group_by(sam, !!sym(cluster_col))
  }

  if (downsample) {
    sam <- dplyr::slice_sample(sam, n = sam_sz)
  }

  # Calculate diversity
  .calc_div <- function(x, met, n_bts = n_boots) {
    fn <- function(x, i) {
      met(table(x[i]))
    }

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

  div <- tidyr::pivot_longer(div, c(.data$diversity, .data$stderr))
  div <- tidyr::unite(div, "name", .data$met, .data$name)
  div <- tidyr::pivot_wider(div)

  # Format results
  vdj <- dplyr::distinct(vdj, !!!syms(vdj_cols))

  if (!is.null(cluster_col)) {
    vdj <- dplyr::left_join(vdj, div, by = cluster_col)

  } else {
    vdj <- dplyr::bind_cols(vdj, div)
  }

  res <- dplyr::left_join(meta, vdj, by = vdj_cols)

  if (return_df) input <- meta

  res <- .add_meta(input, meta = res)

  res
}


#' Plot repertoire diversity
#'
#' @param input Single cell object or data.frame containing V(D)J data. If a
#' data.frame is provided, the cell barcodes should be stored as row names.
#' @param cluster_col meta.data column containing cluster IDs to use for
#' grouping cells when calculating clonotype abundance
#' @param method Function to use for calculating diversity. A named list of
#' functions can be passed to plot multiple diversity metrics,
#' e.g. list(simpson = abdiv::simpson, shannon = abdiv::shannon)
#' @param downsample Downsample clusters to the same size when calculating
#' diversity metrics
#' @param n_boots Number of bootstrap replicates
#' @param data_col meta.data column containing clonotype IDs to use for
#' calculating repertoire diversity
#' @param group_col meta.data column to use for grouping IDs present in
#' cluster_col
#' @param plot_colors Character vector containing colors for plotting
#' @param plot_lvls Character vector containing levels for ordering
#' @param facet_rows The number of facet rows, use this argument if a list of
#' functions is passed to method
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
plot_diversity <- function(input, cluster_col = NULL, group_col = NULL,
                           method = abdiv::simpson, downsample = FALSE,
                           n_boots = 100, data_col = "clonotype_id",
                           plot_colors = NULL, plot_lvls = NULL,
                           facet_rows = 1, sep = ";", ...) {

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

  # Diversity columns
  div_cols <- paste0(names(method), "_", "diversity")
  err_cols <- paste0(names(method), "_", "stderr")

  plt_cols <- c(cluster_col, group_col, div_cols, err_cols)

  # Calculate diversity
  plt_dat <- .get_meta(input)
  plt_dat <- .prepare_meta(input = NULL, plt_dat)
  plt_dat <- dplyr::select(plt_dat, -any_of(c(div_cols, err_cols)))

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
  plt_dat <- tidyr::pivot_longer(plt_dat, all_of(c(div_cols, err_cols)))

  re <- "^(.+)_(diversity|stderr)$"

  plt_dat <- tidyr::extract(plt_dat, .data$name, into = c("met", "type"), re)

  plt_dat <- tidyr::pivot_wider(
    plt_dat,
    names_from  = .data$type,
    values_from = .data$value
  )

  # Set plot levels
  plt_dat <- .set_lvls(plt_dat, cluster_col, plot_lvls)

  include_x_labs <- !is.null(cluster_col)

  if (is.null(cluster_col)) {
    cluster_col <- "met"
  }

  # Create grouped boxplot
  if (!is.null(group_col) && !is.null(cluster_col)) {
    plt_dat <- .set_lvls(plt_dat, group_col, plot_lvls)

    res <- .create_boxes(
      plt_dat,
      x      = group_col,
      y      = "diversity",
      .color = group_col,
      .fill  = group_col,
      alpha  = 0.5,
      clrs   = plot_colors,
      outlier.color = NA,
      ...
    ) +
      ggplot2::geom_jitter(
        position = ggplot2::position_jitterdodge(jitter.width = 0.05)
      ) +
      ggplot2::theme(legend.position = "right")

    # Create facets
    if (length(method) > 1) {
      res <- res +
        ggplot2::facet_wrap(~ met, nrow = facet_rows, scales = "free") +
        ggplot2::theme(axis.title.y = ggplot2::element_blank())

    } else {
      res <- res +
        ggplot2::labs(y = names(method))
    }

    return(res)
  }

  # Create bar graphs
  res <- ggplot2::ggplot(
    plt_dat,
    ggplot2::aes(!!sym(cluster_col), .data$diversity, fill = !!sym(cluster_col))
  ) +
    ggplot2::geom_col(...) +

    ggplot2::geom_linerange(
      aes(
        !!sym(cluster_col),
        ymin = .data$diversity - .data$stderr,
        ymax = .data$diversity + .data$stderr
      )
    )

  # Create facets
  if (length(method) > 1) {
    res <- res +
      ggplot2::facet_wrap(~ met, nrow = facet_rows, scales = "free") +
      ggplot2::theme(axis.title.y = ggplot2::element_blank())

  } else {
    res <- res +
      ggplot2::labs(y = names(method))
  }

  # Set plot colors
  if (!is.null(plot_colors)) {
    res <- res +
      ggplot2::scale_fill_manual(values = plot_colors)
  }

  # Set theme
  res <- res +
    djvdj_theme() +
    ggplot2::theme(
      legend.position = "none",
      axis.title.x    = ggplot2::element_blank()
    )

  if (!include_x_labs) {
    res <- res +
      ggplot2::theme(
        axis.text.x  = ggplot2::element_blank(),
        axis.ticks.x = ggplot2::element_blank()
      )
  }

  res
}

