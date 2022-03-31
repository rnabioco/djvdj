#' Calculate repertoire diversity
#'
#' @param input Single cell object or data.frame containing V(D)J data. If a
#' data.frame is provided, the cell barcodes should be stored as row names.
#' @param cluster_col meta.data column containing cluster IDs to use for
#' grouping cells when calculating diversity. If cluster_col is omitted,
#' diversity index will be calculated using all cells.
#' @param method Method to use for calculating diversity. A named list can also
#' be passed to use multiple methods. The names should specify names for the
#' output columns.
#' @param clonotype_col meta.data column containing clonotype IDs to use for
#' calculating diversity
#' @param prefix Prefix to add to new columns
#' @param return_df Return results as a data.frame. If set to FALSE, results
#' will be added to the input object.
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
calc_diversity <- function(input, cluster_col = NULL, method = abdiv::simpson,
                           clonotype_col = "clonotype_id", n_boots = 100, prefix = "",
                           return_df = FALSE) {

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

  vdj_cols <- c(clonotype_col, cluster_col)

  meta <- .get_meta(input)
  vdj  <- dplyr::filter(meta, !is.na(!!sym(clonotype_col)))

  if (!is.null(cluster_col) && is.null(meta[[cluster_col]])) {
    stop(cluster_col, " not found in object, provide a different cluster_col")
  }

  # Get parameters for bootstrapped sampling
  sam     <- vdj
  sam_rng <- c(50, 5000)
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

  sam <- dplyr::slice_sample(sam, n = sam_sz)

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
      diversity = list(.calc_div(!!sym(clonotype_col), met = .x)),
      stderr    = purrr::map_dbl(diversity, pull, "std.error"),
      diversity = purrr::map_dbl(diversity, pull, "statistic")
    )
  })

  div <- tidyr::pivot_longer(div, c(diversity, stderr))
  div <- tidyr::unite(div, name, met, name)
  div <- tidyr::pivot_wider(div)

  # Format results
  vdj <- dplyr::distinct(vdj, !!!syms(vdj_cols))

  if (!is.null(cluster_col)) {
    vdj <- dplyr::left_join(vdj, div, by = cluster_col)

  } else {
    vdj <- dplyr::bind_cols(vdj, div)
  }

  res <- dplyr::left_join(meta, vdj, by = vdj_cols)

  if (return_df) {
    input <- meta
  }

  res <- .add_meta(input, meta = res)

  res
}

# .calc_div <- function(clone_col, met) {
#   cnt <- table(clone_col)
#
#   met(cnt)
# }
#
# res <- purrr::map_dfr(1:n_boots, ~ {
#   idx <- .x
#   sam <- dplyr::slice_sample(vdj, n = sam_sz, replace = TRUE)
#
#   method %>%
#     purrr::imap_dfr(~ {
#       dplyr::summarize(
#         sam,
#         div = .calc_div(!!sym(clonotype_col), .x),
#         met = .y,
#         idx = idx
#       )
#     })
# })

# calc_diversity <- function(input, cluster_col = NULL, method = abdiv::simpson,
#                            clonotype_col = "clonotype_id", n_boots = 100, prefix = "",
#                            return_df = FALSE) {
#
#   if (length(method) > 1 && is.null(names(method))) {
#     stop("Must include names if using a list of methods.")
#   }
#
#   if (length(method) == 1 && is.null(names(method))) {
#     nm <- as.character(substitute(method))
#     nm <- last(nm)
#
#     method <- purrr::set_names(list(method), nm)
#   }
#
#   # Format input data
#   meta <- .get_meta(input)
#   vdj  <- dplyr::filter(meta, !is.na(!!sym(clonotype_col)))
#
#   # Count clonotypes
#   vdj_cols <- clonotype_col
#
#   if (!is.null(cluster_col)) {
#     vdj_cols <- c(cluster_col, vdj_cols)
#     vdj   <- dplyr::group_by(vdj, !!sym(cluster_col))
#   }
#
#   vdj <- dplyr::group_by(vdj, !!sym(clonotype_col), .add = TRUE)
#
#   vdj <- dplyr::summarize(
#     vdj,
#     .n      = dplyr::n_distinct(!!sym(CELL_COL)),
#     .groups = "drop"
#   )
#
#   # Calculate diversity
#   if (!is.null(cluster_col)) {
#     vdj <- dplyr::group_by(vdj, !!sym(cluster_col))
#   }
#
#   if (length(method) == 1 && is.null(names(method))) {
#     nm <- as.character(substitute(method))
#     nm <- dplyr::last(nm)
#
#     method        <- list(method)
#     names(method) <- nm
#   }
#
#   div_cols      <- paste0(prefix, names(method))
#   names(method) <- div_cols
#
#   vdj <- purrr::imap_dfr(method, ~ {
#     dplyr::mutate(
#       vdj,
#       diversity = .x(.data$.n),
#       met       = .y
#     )
#   })
#
#   vdj <- tidyr::pivot_wider(
#     vdj,
#     names_from  = "met",
#     values_from = "diversity"
#   )
#
#   vdj <- dplyr::ungroup(vdj)
#   vdj <- dplyr::select(vdj, all_of(c(vdj_cols, div_cols)))
#
#   # Format results
#   res <- dplyr::left_join(meta, vdj, by = vdj_cols)
#
#   if (return_df) {
#     input <- meta
#   }
#
#   res <- .add_meta(input, meta = res)
#
#   res
# }


#' Plot repertoire diversity
#'
#' @param input Single cell object or data.frame containing V(D)J data. If a
#' data.frame is provided, the cell barcodes should be stored as row names.
#' @param cluster_col meta.data column containing cluster IDs to use for
#' grouping cells when calculating clonotype abundance
#' @param method Function to use for calculating diversity. A named list of
#' functions can be passed to plot multiple diversity metrics,
#' e.g. list(simpson = abdiv::simpson, shannon = abdiv::shannon)
#' @param clonotype_col meta.data column containing clonotype IDs to use for
#' calculating repertoire diversity
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
plot_diversity <- function(input, cluster_col = NULL, method = abdiv::simpson, clonotype_col = "clonotype_id",
                           plot_colors = NULL, plot_lvls = NULL, facet_rows = 1, ...) {

  browser()

  if (length(method) > 1 && is.null(names(method))) {
    stop("Must include names if providing a list of methods.")
  }

  if (length(method) == 1 && is.null(names(method))) {
    nm <- as.character(substitute(method))
    nm <- dplyr::last(nm)

    method        <- list(method)
    names(method) <- nm
  }

  # Calculate diversity
  plt_dat <- calc_diversity(
    input         = input,
    cluster_col   = cluster_col,
    method        = method,
    clonotype_col = clonotype_col,
    prefix        = "",
    return_df     = TRUE
  )

  div_cols <- names(method)
  div_cols <- paste0(div_cols, "_", c("diversity", "stderr"))

  plt_dat <- dplyr::filter(plt_dat, !is.na(!!sym(clonotype_col)))
  plt_dat <- dplyr::distinct(plt_dat, !!!syms(c(cluster_col, div_cols)))

  plt_dat <- tidyr::pivot_longer(plt_dat, all_of(div_cols))

  plt_dat <- tidyr::separate(plt_dat, name, into = c("met", "type"))  # THIS CAN NOT DEPEND ON '_' SEPARATOR

  # Format data for plotting
  include_strips <- TRUE

  if (is.null(cluster_col)) {
    plt_dat <- tibble::tibble(
      div  = unlist(plt_dat, use.names = FALSE),
      name = names(plt_dat),
    )

    include_strips <- FALSE
    cluster_col <- "name"

  } else {
    plt_dat <- tidyr::pivot_longer(
      plt_dat,
      cols      = all_of(names(method)),
      values_to = "div"
    )

    plt_dat <- .set_lvls(plt_dat, cluster_col, plot_lvls)
  }

  # Set plot levels
  plt_dat <- .set_lvls(plt_dat, "name", plt_dat$name)

  # Create bar graphs
  res <- ggplot2::ggplot(
    plt_dat,
    ggplot2::aes(!!sym(cluster_col), .data$div, fill = !!sym(cluster_col))
  ) +
    ggplot2::geom_col(...)

  # Create facets
  if (length(method) > 1) {
    res <- res +
      ggplot2::facet_wrap(~ name, nrow = facet_rows, scales = "free") +
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

  if (!include_strips) {
    res <- res +
      ggplot2::theme(strip.text = ggplot2::element_blank())
  }

  res
}



