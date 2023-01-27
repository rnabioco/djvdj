#' Plot numerical data
#'
#' @param input Single cell object or data.frame, if a
#' data.frame is provided, cell barcodes should be stored as row names.
#' @param data_col meta.data column containing data to plot
#' @param cluster_col meta.data column containing cluster IDs to use for
#' grouping cells for plotting
#' @param group_col meta.data column to use for grouping clusters into separate
#' panels
#' @param method Method to use for plotting, possible values include:
#'
#' - 'histogram'
#' - 'density'
#' - 'boxplot'
#' - 'violin'
#'
#' @param units Units to use for y-axis when method is 'histogram'. Use
#' 'frequency' to show the number of values or 'percent' to show the percentage
#' of total values.
#'
#' ## Plot aesthetics
#'
#' @param plot_colors Character vector specifying colors to use for cell
#' clusters specified by cluster_col. When cluster_col is NULL, plot colors can
#' be directly modified with the ggplot2 parameters color and fill,
#' e.g. fill = "red", color = "black"
#' @param plot_lvls Character vector containing order to use for plotting cell
#' clusters specified by cluster_col
#' @param trans Transformation to use when plotting data, e.g. 'log10'. By
#' default values are not transformed, refer to [ggplot2::continuous_scale()]
#' for more options.
#' @param panel_nrow The number of rows to use for arranging plot panels
#' @param panel_scales Should scales for plot panels be fixed or free. This
#' passes a scales specification to ggplot2::facet_wrap, can be 'fixed', 'free',
#' 'free_x', or 'free_y'. 'fixed' will cause panels to share the same scales.
#' @param n_label Location on plot where n label should be added, this can be
#' any combination of the following:
#'
#' - 'corner', display the total number of cells plotted in the top right
#'   corner, the position of the label can be modified by passing `x` and `y`
#'   specifications with the `label_params` argument
#' - 'axis', display the number of cells plotted for each group shown on the
#'   x-axis
#' - 'legend', display the number of cells plotted for each group shown in the
#'   plot legend
#' - 'none', do not display the number of cells plotted
#'
#' @param label_params Named list providing additional parameters to modify
#' n label aesthetics, e.g. list(size = 4, color = "red")
#' @param ... Additional arguments to pass to ggplot2, e.g. color, fill, size,
#' linetype, etc.
#'
#' ## VDJ arguments
#'
#' @param per_chain If `TRUE` values will be plotted for each chain, i.e. each
#' data point represents a chain. If `FALSE` values will be summarized for each
#' cell using `summary_fn` before plotting, i.e. each data point represents a
#' cell.
#' @param chain Chain(s) to use for filtering data before plotting. If NULL
#' data will not be filtered based on chain.
#' @param chain_col meta.data column containing chains for each cell
#' @param summary_fn Function to use for summarizing values when `per_chain` is
#' `FALSE`, can be either a function, e.g. `mean`, or a purrr-style
#' lambda, e.g. `~ mean(.x, na.rm = TRUE)` where `.x` refers to the column. If
#' `NULL`, the mean will be calculated.
#' @param sep Separator used for storing per-chain V(D)J data for each cell
#'
#' @seealso [summarize_vdj()] for more examples on how per-chain data can be
#' summarized for each cell
#' @return ggplot object
#' @name plot_numerical
NULL

#' @rdname plot_numerical
#' @export
plot_histogram <- function(input, data_col, cluster_col = NULL,
                           group_col = NULL, method = "histogram",
                           units = "frequency", plot_colors = NULL,
                           plot_lvls = names(plot_colors), trans = "identity",
                           panel_nrow = NULL, panel_scales = "fixed",
                           n_label = NULL, label_params = list(), ...,
                           per_chain = TRUE, chain = NULL,
                           chain_col = global$chain_col, summary_fn = mean,
                           sep = global$sep) {

  # Check input classes
  .check_args()

  # Check input values
  .check_group_cols(cluster_col, group_col, input)

  # Format plot data
  plt_dat <- .format_plot_data(
    input,
    data_col    = data_col,
    cluster_col = cluster_col,
    plot_lvls   = plot_lvls,
    chain       = chain,
    chain_col   = chain_col,
    per_chain   = per_chain,
    summary_fn  = summary_fn,
    sep         = sep,

    # additional variables to fetch from object
    group_col
  )

  # Check for sep
  # Need this info to set y-axis title
  has_sep <- .detect_sep(
    .get_meta(input),
    data_col, sep, n_rows = 1000, return_names = FALSE
  )

  if (!has_sep) per_chain <- FALSE

  # Create histogram
  gg_args <- list(
    df_in        = plt_dat,
    x            = data_col,
    grp          = group_col,
    .color       = cluster_col,
    .fill        = cluster_col,
    method       = method,
    units        = units,
    trans        = trans,
    nrow         = panel_nrow,
    scales       = panel_scales,
    n_label      = n_label,
    label_params = label_params,
    ...
  )

  gg_args$alpha <- gg_args$alpha %||% 0.5

  gg_args$y_ttl <- .get_axis_label(
    units, sfx = ifelse(per_chain, "chains", "cells")
  )

  res <- lift(.create_hist)(gg_args)

  res
}

#' @rdname plot_numerical
#' @export
plot_violin <- function(input, data_col, cluster_col = NULL, group_col = NULL,
                        method = "violin", plot_colors = NULL,
                        plot_lvls = names(plot_colors), trans = "identity",
                        panel_nrow = NULL, panel_scales = "free_x",
                        n_label = NULL, label_params = list(), ...,
                        per_chain = TRUE, chain = NULL,
                        chain_col = global$chain_col, summary_fn = mean,
                        sep = global$sep) {

  # Check input classes
  .check_args()

  # Check input values
  .check_group_cols(cluster_col, group_col, input)

  # Format plot data
  plt_dat <- .format_plot_data(
    input,
    data_col    = data_col,
    cluster_col = cluster_col,
    plot_lvls   = plot_lvls,
    chain       = chain,
    chain_col   = chain_col,
    per_chain   = per_chain,
    summary_fn  = summary_fn,
    sep         = sep,

    # additional variables to fetch from object
    group_col
  )

  # Create violins
  gg_args <- list(
    df_in        = plt_dat,
    x            = cluster_col,
    y            = data_col,
    grp          = group_col,
    .color       = cluster_col,
    .fill        = cluster_col,
    method       = method,
    trans        = trans,
    nrow         = panel_nrow,
    scales       = panel_scales,
    show_points  = identical(method, "violin"),
    n_label      = n_label,
    label_params = label_params,
    ...
  )

  gg_args$point.size <- gg_args$point.size %||% 1
  gg_args$alpha      <- gg_args$alpha %||% 0.5

  res <- lift(.create_boxes)(gg_args)

  res
}

#' Create scatter plot
#'
#' @export
plot_scatter <- function(input, data_col = NULL, x = "UMAP_1", y = "UMAP_2",
                         group_col = NULL, summary_fn = NULL, chain = NULL,
                         chain_col = global$chain_col, top = NULL,
                         plot_colors = NULL, plot_lvls = names(plot_colors),
                         trans = "identity", min_q = NULL, max_q = NULL,
                         panel_nrow = NULL, panel_scales = "fixed",
                         na_color = "grey80", other_label = "other",
                         n_label = NULL, label_params = list(),
                         sep = global$sep, ...) {

  # Check input classes
  .check_args()

  # Format plot data
  plt_dat <- .format_plot_data(
    input,
    data_col    = data_col,
    cluster_col = group_col,
    plot_lvls   = plot_lvls,
    per_chain   = FALSE,
    chain       = chain,
    chain_col   = chain_col,
    summary_fn  = summary_fn,
    sep         = sep,

    # additional variables to fetch from object
    x, y
  )

  # Set plot arguments
  gg_args <- list(
    fn           = ggplot2::geom_point,
    x            = x,
    y            = y,
    .color       = data_col,
    grp          = group_col,
    na_clr       = na_color,
    nrow         = panel_nrow,
    scales       = panel_scales,
    trans_clr    = trans,
    n_label      = n_label,
    label_params = label_params,
    ...
  )

  rev_lgnd <- FALSE

  set_other_color <- !other_label %in% names(plot_colors)

  # If no data_col provided, return scatter plot
  if (is.null(data_col)) {
    gg_args$df_in <- plt_dat

    res <- lift(.create_scatter)(gg_args)

    return(res)
  }

  # Adjust values based on min_q and max_q
  is_num <- is.numeric(plt_dat[[data_col]])

  if (is_num && (!is.null(min_q) || !is.null(max_q))) {
    plt_dat <- .set_lims(plt_dat, data_col, min_q, max_q)
  }

  # Convert logical to character
  if (is.logical(plt_dat[[data_col]])) {
    plt_dat <- purrr::modify_at(plt_dat, data_col, as.character)

    plot_lvls <- plot_lvls %||% c("TRUE", "FALSE")
  }

  # Set data_col colors and levels, use ggplot2 defaults if plot_colors is NULL
  if (is_num) {
    plot_colors <- plot_colors %||% c("#132B43", "#56B1F7")

  } else {
    # Set other group
    if (!is.null(top)) {
      plt_dat <- .set_other_grps(
        plt_dat, data_col, top = top, other_label = other_label
      )

      if (is.null(plot_lvls)) {
        rnk_lvls  <- levels(plt_dat[[data_col]])
        plot_lvls <- stats::na.omit(rev(rnk_lvls))
        rev_lgnd  <- TRUE

      } else {
        plot_lvls <- c(other_label, plot_lvls)
      }
    }

    # Set default levels if plot_lvls is NULL
    if (is.null(plot_lvls)) {
      dat       <- plt_dat[[data_col]]
      plot_lvls <- levels(dat) %||% sort(unique(dat))
    }

    # Set color names, need to do this so na_color is not overridden if the
    # user provides extra values for plot_colors
    plot_colors <- plot_colors %||% scales::hue_pal()(length(plot_lvls))

    if (is.null(names(plot_colors))) {
      lvls <- stats::na.omit(plot_lvls)

      plot_colors <- plot_colors[seq_along(lvls)]
      plot_colors <- stats::na.omit(plot_colors)

      names(plot_colors) <- lvls[seq_along(plot_colors)]
    }

    if (set_other_color) plot_colors[[other_label]] <- "grey60"
  }

  plt_dat <- .set_lvls(plt_dat, data_col, plot_lvls)

  # Sort data so largest values are plotted on top
  # Do not use arrange() since NAs always get sorted to bottom
  # res <- arrange(plt_dat, !!sym(data_col))
  plt_dat <- plt_dat[order(plt_dat[[data_col]], na.last = FALSE), ]

  gg_args$df_in <- plt_dat
  gg_args$clrs  <- plot_colors

  res <- lift(.create_scatter)(gg_args)

  # Adjust legend point size
  if (!is_num) {
    lgnd_args <- list(override.aes = list(size = 3))

    if (rev_lgnd) lgnd_args$reverse <- TRUE

    res <- res +
      ggplot2::guides(color = lift(ggplot2::guide_legend)(lgnd_args))
  }

  res
}



#' Format data for plotting
#'
#' @param input Single cell object or data.frame, if a
#' data.frame is provided, cell barcodes should be stored as row names.
#' @param data_col meta.data column containing data to plot
#' @param ... Names of additional columns to check for in input
#' @param cluster_col meta.data column containing cluster IDs to use for
#' grouping cells for plotting
#' @param plot_lvls Character vector containing order to use for plotting cell
#' clusters specified by cluster_col
#' @param chain Chain(s) to use for filtering data before plotting. If NULL
#' data will not be filtered based on chain.
#' @param chain_col meta.data column containing chains for each cell
#' @param per_chain If `TRUE` values will be plotted for each chain, i.e. each
#' data point represents a chain. If `FALSE` values will be summarized for each
#' cell using `summary_fn` before plotting, i.e. each data point represents a
#' cell.
#' @param summary_fn Function to use for summarizing values when `per_chain` is
#' `FALSE`, can be either a function, e.g. `mean`, or a purrr-style
#' lambda, e.g. `~ mean(.x, na.rm = TRUE)` where `.x` refers to the column. If
#' `NULL`, the mean will be calculated.
#' @param filter_cells If `TRUE` cells not containing V(D)J data (represented
#' as `NA`s) will be removed.
#' @param slot Slot to pull data from when `input` is a Seurat object
#' @param sep Separator used for storing per-chain V(D)J data for each cell
#' @return data.frame formatted for plotting
#' @noRd
.format_plot_data <- function(input, ...) {

  UseMethod(".format_plot_data", input)
}

.format_plot_data.default <- function(input, data_col, ..., cluster_col = NULL,
                                      plot_lvls = NULL, chain = NULL,
                                      chain_col = global$chain_col,
                                      per_chain = TRUE, summary_fn = mean,
                                      filter_cells = per_chain,
                                      sep = global$sep) {

  # Check that columns are present in object
  .check_obj_cols(input, data_col, ..., chain = chain, chain_col = chain_col)

  # Format plot data
  res <- .get_meta(input)

  has_sep <- .detect_sep(
    res, data_col, sep, n_rows = 1000, return_names = FALSE
  )

  # Summarize per-chain data
  # Only remove NAs when plotting per-chain data
  if (has_sep) {
    if (per_chain) {
      if (is.null(chain)) fetch_cols <- data_col
      else                fetch_cols <- c(data_col, chain_col)

      res <- fetch_vdj(
        res,
        data_cols     = fetch_cols,
        clonotype_col = NULL,
        unnest        = TRUE
      )

      if (!is.null(chain)) {
        res <- dplyr::filter(res, !!sym(chain_col) %in% chain)
      }

    } else {
      res <- summarize_vdj(
        res,
        data_cols = data_col,
        fn        = summary_fn,
        chain     = chain,
        chain_col = chain_col,
        sep       = sep,
        return_df = TRUE
      )
    }
  }

  # Remove cells with missing V(D)J data
  if (filter_cells) res <- dplyr::filter(res, !is.na(!!sym(data_col)))

  # Order clusters based on plot_lvls
  res <- .set_lvls(res, cluster_col, plot_lvls)

  res
}

.format_plot_data.Seurat <- function(input, data_col, ..., cluster_col = NULL,
                                     plot_lvls = NULL, chain = NULL,
                                     chain_col = global$chain_col,
                                     per_chain = TRUE, summary_fn = mean,
                                     filter_cells = per_chain,
                                     slot = "data", sep = global$sep) {

  # Fetch variables and add to meta.data
  # want input data to include meta.data and any features from FetchData
  plt_vars <- c(data_col, ...)

  plt_dat <- Seurat::FetchData(
    input,
    vars = unique(plt_vars),
    slot = slot
  )

  # Format input data
  meta <- .get_meta(input)

  plt_dat <- .merge_meta(meta, plt_dat)

  res <- .format_plot_data(
    plt_dat,
    data_col     = data_col,
    cluster_col  = cluster_col,
    plot_lvls    = plot_lvls,
    chain        = chain,
    chain_col    = chain_col,
    per_chain    = per_chain,
    summary_fn   = summary_fn,
    filter_cells = filter_cells,
    sep          = sep
  )

  res
}
