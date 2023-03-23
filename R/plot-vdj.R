#' Plot numerical single-cell data
#'
#' Visualize numerical single-cell data by creating a histogram, density plot,
#' violin plots, or boxplots. When plotting V(D)J data, values can be plotted
#' separately for each chain or summarized and plotted for each cell.
#'
#' @param input Single cell object or data.frame, if a
#' data.frame is provided, cell barcodes should be stored as row names.
#' @param data_col meta.data column containing data to plot
#' @param cluster_col meta.data column containing cluster IDs to use for
#' grouping cells for plotting
#' @param group_col meta.data column to use for grouping clusters into separate
#' panels
#' @param method Method to use for plotting, possible values include,
#' 'histogram', 'density', 'boxplot', or 'violin'
#'
#' @param units Units to use for y-axis when method is 'histogram'. Use
#' 'frequency' to show the number of values or 'percent' to show the percentage
#' of total values.
#' @param top To only show the top cell groups, provide
#' one of the following, all other cells will be labeled using the value
#' provided to the `other_label` argument. If `NULL` this will be automatically
#' set.
#'
#' - Integer specifying the number of top groups to show
#' - Vector specifying the names of cell groups to show
#'
#' @param other_label Label to use for 'other' cells when `top` is specified, if
#' `NULL` all cell groups will be shown.
#'
#' ## Aesthetics
#'
#' @param plot_colors Character vector specifying colors to use for cell
#' clusters specified by cluster_col. When cluster_col is `NULL`, plot colors
#' can be directly modified with the ggplot2 parameters `color` and `fill`,
#' e.g. `fill = "red"`, `color = "black"`
#' @param plot_lvls Character vector containing order to use for plotting cell
#' clusters specified by cluster_col
#' @param trans Transformation to use when plotting data, e.g. 'log10'. By
#' default values are not transformed, refer to [ggplot2::continuous_scale()]
#' for more options.
#' @param panel_nrow The number of rows to use for arranging plot panels
#' @param panel_scales Should scales for plot panels be fixed or free. This
#' passes a scales specification to [ggplot2::facet_wrap()], can be 'fixed',
#' 'free', 'free_x', or 'free_y'. 'fixed' will cause panels to share the same
#' scales.
#' @param na_color Color to use for missing values. If plotting V(D)J data,
#' cells lacking data will be plotted as `NA`s.
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
#' @seealso [summarize_vdj()] for more examples on how per-chain data can be
#' summarized for each cell
#' @return ggplot object
#' @name plot_numerical
NULL

#' @rdname plot_numerical
#' @export
plot_histogram <- function(input, data_col, cluster_col = NULL,
                           group_col = NULL, method = "histogram",
                           top = NULL, other_label = "other",
                           units = "frequency", plot_colors = NULL,
                           plot_lvls = NULL, trans = "identity",
                           panel_nrow = NULL, panel_scales = "fixed",
                           na_color = "grey80",
                           n_label = NULL, label_params = list(), ...,
                           per_chain = FALSE, chain = NULL,
                           chain_col = global$chain_col, summary_fn = mean,
                           sep = global$sep) {

  # Check input classes
  .check_args()

  # Format plot data
  plt_dat <- .format_plot_data(
    input,
    data_col    = data_col,
    cluster_col = cluster_col,
    group_col   = group_col,
    plot_lvls   = plot_lvls,
    top         = top,
    other_label = other_label,
    chain       = chain,
    chain_col   = chain_col,
    per_chain   = per_chain,
    summary_fn  = summary_fn,
    sep         = sep,

    # additional variables to fetch from object
    group_col
  )

  plot_colors <- .set_colors(plt_dat, cluster_col, plot_colors,
                             other_label = other_label)

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
    clrs         = plot_colors,
    nrow         = panel_nrow,
    scales       = panel_scales,
    na_clr       = na_color,
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
                        method = "violin", top = NULL, other_label = "other",
                        plot_colors = NULL, plot_lvls = NULL,
                        trans = "identity", panel_nrow = NULL,
                        panel_scales = "free_x", na_color = "grey80",
                        n_label = NULL, label_params = list(), ...,
                        per_chain = FALSE, chain = NULL,
                        chain_col = global$chain_col, summary_fn = mean,
                        sep = global$sep) {

  # Check input classes
  .check_args()

  # Format plot data
  plt_dat <- .format_plot_data(
    input,
    data_col    = data_col,
    cluster_col = cluster_col,
    group_col   = group_col,
    plot_lvls   = plot_lvls,
    top         = top,
    other_label = other_label,
    rank_method = "boxplot",
    chain       = chain,
    chain_col   = chain_col,
    per_chain   = per_chain,
    summary_fn  = summary_fn,
    sep         = sep,

    # additional variables to fetch from object
    group_col
  )

  plot_colors <- .set_colors(plt_dat, cluster_col, plot_colors,
                             other_label = other_label)

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
    clrs         = plot_colors,
    nrow         = panel_nrow,
    scales       = panel_scales,
    na_clr       = na_color,
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


#' Create 2D scatter plot
#'
#' Create a scatter plot with cells colored based on the provided feature.
#' This can be used to create dimensional reduction plots (e.g. UMAP, tSNE, PCA)
#' or to compare different cell attributes (e.g. CD4 vs CD8 expression).
#' When plotting V(D)J data, per-chain values will be summarized for each cell.
#'
#' @param input Single cell object or data.frame, if a
#' data.frame is provided, cell barcodes should be stored as row names.
#' @param data_col Name of meta.data column or other variable (e.g. gene name)
#' to use for coloring points
#' @param x,y Name of meta.data column or other variable to plot on x and
#' y-axis
#' @param group_col meta.data column to use for splitting plot into panels
#' @param data_slot Slot to pull data from when `input` is a Seurat object
#' @param top To only show the top cell groups, provide
#' one of the following, all other cells will be labeled using the value
#' provided to the `other_label` argument. If `NULL` this will be automatically
#' set.
#'
#' - Integer specifying the number of top groups to show
#' - Vector specifying the names of cell groups to show
#'
#' @param other_label Label to use for 'other' cells when `top` is specified, if
#' `NULL` all cell groups will be shown.
#'
#' ## Aesthetics
#'
#' @param plot_colors Character vector specifying colors to use for cell
#' clusters specified by cluster_col.
#' @param plot_lvls Character vector containing order to use for plotting cell
#' clusters specified by cluster_col.
#' @param outline Add an outline around each cluster, outline aesthetics can be
#' modified by passing arguments directly to [ggtrace::geom_point_trace()]
#' @param trans Transformation to use when plotting data, e.g. 'log10'. By
#' default values are not transformed, refer to [ggplot2::continuous_scale()]
#' for more options.
#' @param panel_nrow The number of rows to use for arranging plot panels
#' @param panel_scales Should scales for plot panels be fixed or free. This
#' passes a scales specification to [ggplot2::facet_wrap()], can be 'fixed',
#' 'free', 'free_x', or 'free_y'. 'fixed' will cause panels to share the same
#' scales.
#' @param min_q Minimum quantile cutoff for color scale.
#' @param max_q Maximum quantile cutoff for color scale.
#' @param na_color Color to use for missing values. If plotting V(D)J data,
#' cells lacking data will be plotted as `NA`s.
#' @param n_label Location on plot where n label should be added, this can be
#' any combination of the following:
#'
#' - 'corner', display the total number of cells plotted in the top right
#'   corner, the position of the label can be modified by passing `x` and `y`
#'   specifications with the `label_params` argument
#' - 'legend', display the number of cells plotted for each group shown in the
#'   plot legend
#' - 'none', do not display the number of cells plotted
#'
#' @param label_params Named list providing additional parameters to modify
#' n label aesthetics, e.g. list(size = 4, color = "red")
#' @param ... Additional arguments to pass to [ggplot2::geom_point()], or
#' [ggtrace::geom_point_trace()] if `outline = TRUE`, e.g. color, size, etc.
#'
#' ## VDJ arguments
#'
#' @param chain Chain(s) to use for filtering data before plotting. If NULL
#' data will not be filtered based on chain.
#' @param chain_col meta.data column containing chains for each cell
#' @param summary_fn Function to use for summarizing per-chain values for each
#' cell, can be either a function, e.g. `mean`, or a purrr-style
#' lambda, e.g. `~ mean(.x, na.rm = TRUE)` where `.x` refers to the column. If
#' `NULL`, the mean will be calculated for numeric values, non-numeric columns
#' will be combined into a single string.
#' @param sep Separator used for storing per-chain V(D)J data for each cell
#' @seealso [summarize_vdj()] for more examples on how per-chain data can be
#' summarized for each cell
#' @return ggplot object
#' @importFrom ggtrace geom_point_trace
#' @export
plot_scatter <- function(input, data_col = NULL, x = "UMAP_1", y = "UMAP_2",
                         group_col = NULL, data_slot = "data", top = NULL,
                         other_label = "other", plot_colors = NULL,
                         plot_lvls = NULL, outline = FALSE, trans = "identity",
                         panel_nrow = NULL, panel_scales = "fixed",
                         min_q = NULL, max_q = NULL, na_color = "grey80",
                         n_label = NULL, label_params = list(), ...,
                         chain = NULL, chain_col = global$chain_col,
                         summary_fn = NULL, sep = global$sep) {

  # Check input classes
  .check_args(data_col = list(allow_null = TRUE))

  # Format plot data
  plt_dat <- .format_plot_data(
    input,
    data_col    = data_col,
    cluster_col = data_col,
    group_col   = group_col,
    plot_lvls   = plot_lvls,
    top         = top,
    other_label = other_label,
    rev_lvls    = TRUE,
    per_chain   = FALSE,
    chain       = chain,
    chain_col   = chain_col,
    summary_fn  = summary_fn,
    slot        = data_slot,
    sep         = sep,

    # additional variables to fetch from object
    x, y
  )

  # Set plot arguments
  gg_args <- list(
    df_in        = plt_dat,
    fn           = ggplot2::geom_point,
    x            = x,
    y            = y,
    .color       = data_col,
    .fill        = data_col,
    grp          = group_col,
    clrs         = plot_colors,
    trans_clr    = trans,
    nrow         = panel_nrow,
    scales       = panel_scales,
    na_clr       = na_color,
    n_label      = n_label,
    label_params = label_params,
    ...
  )

  # Adjust plot arguments for geom_point_trace
  if (outline) {
    gg_args <- .standardize_aes(gg_args)

    gg_args$fn     <- ggtrace::geom_point_trace
    gg_args$colour <- gg_args$colour %||% "black"
  }

  # If no data_col provided, return scatter plot
  if (is.null(data_col)) {
    res <- lift(.create_scatter)(gg_args)

    return(res)
  }

  # Is data_col numeric
  is_num <- is.numeric(plt_dat[[data_col]])

  # Adjust values based on min_q and max_q
  if (is_num && (!is.null(min_q) || !is.null(max_q))) {
    plt_dat <- .set_lims(plt_dat, data_col, min_q, max_q)
  }

  # If not all levels are specified by plot_lvls, use ordering set by
  # .format_plot_data, reverse legend so top ranked levels are listed first
  # ACTUALLY SHOULD ALWAYS REVERSE LEGEND
  lvls <- stats::na.omit(levels(plt_dat[[data_col]]))

  rev_lgnd <- !identical(data_col, group_col)

  if (!identical(data_col, group_col)) lvls <- rev(lvls)

  # Set default colors
  plot_colors <- .set_colors(plt_dat, data_col, plot_colors, lvls,
                             other_label = other_label)

  # Create scatter plot
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
                                      group_col = NULL, plot_lvls = NULL,
                                      top = NULL, other_label = "other",
                                      rank_method = n, rev_lvls = FALSE,
                                      chain = NULL,
                                      chain_col = global$chain_col,
                                      per_chain = TRUE, summary_fn = mean,
                                      filter_cells = per_chain, slot = "data",
                                      sep = global$sep) {

  # Check that columns are present in object
  .check_obj_cols(input, data_col, ..., chain = chain, chain_col = chain_col)

  # Format plot data
  res <- .get_meta(input)

  if (is.null(data_col)) return(res)

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
  } else if (!is.null(chain)) {
    cli::cli_warn(
      "`data_col` does not contain per-chain data, chains were not filtered."
    )
  }

  # Remove cells with missing V(D)J data
  if (filter_cells) res <- dplyr::filter(res, !is.na(!!sym(data_col)))

  # Set column to use for ranking cluster_col and group_col
  if (identical(rank_method, "boxplot")) val_col <- data_col
  else                                   val_col <- NULL

  # Set cluster_col levels
  res <- .set_other_grps(
    res, cluster_col, val_col = val_col, plot_lvls = plot_lvls, top = top,
    other_label = other_label, method = rank_method, rev = rev_lvls
  )

  # Set group_col levels
  # If group_col is not same as cluster_col, use levels from df_in for
  # default ordering
  # If not, re-set group_col levels, since these should never be reversed
  # Only reverse levels when group_col and cluster_col are the same
  if (!is.null(group_col) && !identical(cluster_col, group_col)) {
    grp_lvls <- levels(res[[group_col]])

  } else {
    grp_lvls <- NULL
  }

  res <- .set_other_grps(
    res, group_col, val_col = val_col, plot_lvls = grp_lvls,
    method = rank_method
  )

  # Sort data so largest values are plotted on top
  # Do not use arrange() since NAs always get sorted to bottom
  # res <- arrange(plt_dat, !!sym(data_col))
  res <- res[order(res[[data_col]], na.last = FALSE), ]

  res
}

.format_plot_data.Seurat <- function(input, data_col, ..., cluster_col = NULL,
                                     group_col = NULL, plot_lvls = NULL,
                                     top = NULL, other_label = "other",
                                     rank_method = n, rev_lvls = FALSE,
                                     chain = NULL, chain_col = global$chain_col,
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
    top          = top,
    other_label  = other_label,
    rank_method  = rank_method,
    rev_lvls     = rev_lvls,
    chain        = chain,
    chain_col    = chain_col,
    per_chain    = per_chain,
    summary_fn   = summary_fn,
    filter_cells = filter_cells,
    sep          = sep
  )

  res
}
