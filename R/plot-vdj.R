#' Create 2D feature plot
#'
#' Create a scatter plot with cells colored based on the provided feature.
#' This can be used to create dimensional reduction plots
#' (e.g. UMAP, tSNE, PCA) or to compare different cell attributes
#' (e.g. CD4 vs CD8 expression).
#'
#' @export
plot_features <- function(input, ...) {

  UseMethod("plot_features", input)
}

#' @param input Single cell object or data.frame containing V(D)J data. If a
#' data.frame is provided, the cell barcodes should be stored as row names.
#' @param x Variable to plot on x-axis
#' @param y Variable to plot on y-axis
#' @param group_col meta.data column to use for splitting plot into panels
#' @param plot_colors Vector of colors to use for plotting
#' @param plot_lvls Levels to use for ordering feature
#' @param trans Transformation to use when coloring cells by a continuous
#' variable, e.g. 'log10'. By default values are not transformed, refer to
#' [ggplot2::continuous_scale()] for more options.
#' @param min_q Minimum quantile cutoff for color scale.
#' @param max_q Maximum quantile cutoff for color scale.
#' @param panel_nrow The number of rows to use for arranging plot panels
#' @param panel_scales Should scales for plot panels be fixed or free. This
#' passes a scales specification to ggplot2::facet_wrap, can be 'fixed', 'free',
#' 'free_x', or 'free_y'. 'fixed' will cause panels to share the same scales.
#' @param na_color Color to use for missing values
#' @param n_label Include a label showing the number of cells plotted
#' @param label_params Named list providing additional parameters to modify
#' clonotype and n label aesthetics, e.g. list(size = 4, color = "red")
#' @param ... Additional arguments to pass to [ggplot2::geom_point()], e.g.
#' color, fill, size, linetype, etc.
#' @return ggplot object
#'
#' @examples
#' # Create UMAP with cells colored based on sample name
#' plot_features(
#'   tiny_so,
#'   feature = "orig.ident"
#' )
#'
#' # Compare UMI counts for each cell with number of genes detected
#' plot_features(
#'   tiny_sce,
#'   feature = "orig.ident",
#'   x = "nFeature_RNA",
#'   y = "nCount_RNA"
#' )
#'
#' @export
#' @name plot_features
plot_features.default <- function(input, feature = NULL, x = "UMAP_1",
                                  y = "UMAP_2", group_col = NULL,
                                  plot_colors = NULL,
                                  plot_lvls = names(plot_colors),
                                  trans = "identity", min_q = NULL,
                                  max_q = NULL, panel_nrow = NULL,
                                  panel_scales = "fixed", na_color = "grey80",
                                  n_label = TRUE, label_params = list(), ...) {

  # Check arguments
  if (identical(x, y)) stop("'x' and 'y' must be different")

  plt_dat <- .get_meta(input)

  # Calculate number of cells for label
  .n <- nrow(plt_dat)

  plt_vars  <- c(x, y, feature)
  miss_vars <- plt_vars[!plt_vars %in% colnames(plt_dat)]

  if (length(miss_vars) > 0) {
    stop(paste(miss_vars, collapse = ", "), " not found in object")
  }

  gg_args <- list(
    x          = x,
    y          = y,
    feat       = feature,
    grp        = group_col,
    pt_args    = list(...),
    p_nrow     = panel_nrow,
    p_scales   = panel_scales,
    n_lab      = n_label,
    lab_params = label_params
  )

  # If no feature provided, return scatter plot
  if (is.null(feature)) {
    gg_args$dat <- plt_dat

    gg_args$pt_args$color <- gg_args$pt_args$color %||% plot_colors

    res <- purrr::lift_dl(.create_scatter)(gg_args)

    return(res)
  }

  # Adjust values based on min_q and max_q
  is_num <- is.numeric(plt_dat[[feature]])

  if ((!is.null(min_q) || !is.null(max_q)) && is_num) {
    plt_dat <- .set_lims(
      plt_dat,
      ft = feature,
      mn = min_q,
      mx = max_q
    )
  }

  # Default continuous colors
  if (is_num) plot_colors <- plot_colors %||% c("#132B43", "#56B1F7")

  # Convert logical to character
  if (is.logical(plt_dat[[feature]])) {
    plt_dat <- purrr::modify_at(plt_dat, feature, as.character)
  }

  # Set feature levels
  # Add NA to plot_lvls so missing values are plotted on the bottom
  if (is.null(plot_lvls) && !is_num) {
    dat <- plt_dat[[feature]]

    plot_lvls <- levels(dat) %||% sort(unique(dat))
  }

  plt_dat <- .set_lvls(plt_dat, feature, plot_lvls)

  # Sort data so largest values are plotted on top
  # Do not use arrange() since NAs always get sorted to bottom
  # res <- arrange(plt_dat, !!sym(feature))
  plt_dat <- plt_dat[order(plt_dat[[feature]], na.last = FALSE), ]

  gg_args$dat <- plt_dat

  res <- purrr::lift_dl(.create_scatter)(gg_args)

  # Set feature colors, use ggplot2 defaults if plot_colors is NULL
  if (is_num) {
    plot_colors <- plot_colors %||% c("#132B43", "#56B1F7")

    res <- res +
      ggplot2::scale_color_gradientn(
        colors = plot_colors,
        na.value = na_color,
        trans    = trans
      )

  } else {
    plot_colors <- plot_colors %||% scales::hue_pal()(length(plot_lvls))

    # Set color names, need to do this so na_color is not overridden if the
    # user provides extra values for plot_colors
    if (is.null(names(plot_colors))) {
      lvls <- stats::na.omit(plot_lvls)

      plot_colors <- plot_colors[seq_along(lvls)]
      plot_colors <- stats::na.omit(plot_colors)

      names(plot_colors) <- lvls[seq_along(plot_colors)]
    }

    res <- res +
      ggplot2::scale_color_manual(
        values   = plot_colors,
        na.value = na_color
      ) +
      ggplot2::guides(
        color = ggplot2::guide_legend(override.aes = list(size = 3))
      )
  }

  res
}

.create_scatter <- function(dat, x, y, feat = feature,
                            grp = group_col, p_nrow = panel_nrow,
                            p_scales = panel_scales, pt_args,
                            n_lab = TRUE, lab_params = list()) {

  # If don't use lift_dl() to set plt_aes, will get a warning when
  # manually adjusting colors: "Scale for 'colour' is already present..."
  plt_aes <- list(x = sym(x), y = sym(y))

  if (!is.null(feat)) plt_aes$color <- sym(feat)

  plt_aes <- purrr::lift_dl(ggplot2::aes)(plt_aes)

  res <- ggplot2::ggplot(dat, plt_aes) +
    purrr::lift_dl(ggplot2::geom_point)(pt_args) +
    djvdj_theme()

  if (!is.null(grp)) {
    res <- res +
      ggplot2::facet_wrap(
        stats::as.formula(paste("~", grp)),
        nrow = p_nrow, scales = p_scales
      )
  }

  if (n_lab) {
    res <- res +
      ggplot2::scale_y_continuous(expand = ggplot2::expansion(c(0.05, 0.1)))

    res <- .add_n_label(res, lab_params, nrow(dat))
  }

  res
}

#' @rdname plot_features
#' @param feature Variable to use for coloring points
#' @param data_slot Slot in the Seurat object to pull data
#' @importFrom Seurat FetchData
#' @export
plot_features.Seurat <- function(input, feature = NULL, x = "UMAP_1",
                                 y = "UMAP_2", group_col = NULL,
                                 data_slot = "data", plot_colors = NULL,
                                 plot_lvls = names(plot_colors),
                                 trans = "identity", min_q = NULL,
                                 max_q = NULL, panel_nrow = NULL,
                                 panel_scales = "fixed", na_color = "grey80",
                                 ...) {

  # Fetch variables and add to meta.data
  # want input data to include meta.data and any features from FetchData
  plt_vars <- c(x, y, feature, group_col)

  plt_dat <- Seurat::FetchData(
    input,
    vars = unique(plt_vars),
    slot = data_slot
  )

  # Format input data
  # keep rownames since default method will create rowname column
  plt_dat <- .merge_meta(input, plt_dat)
  plt_dat <- .get_meta(plt_dat)

  res <- plot_features(
    input        = plt_dat,
    x            = x,
    y            = y,
    feature      = feature,
    group_col    = group_col,
    plot_colors  = plot_colors,
    plot_lvls    = plot_lvls,
    trans        = trans,
    min_q        = min_q,
    max_q        = max_q,
    panel_nrow   = panel_nrow,
    panel_scales = panel_scales,
    na_color     = na_color,
    ...
  )

  res
}

#' @rdname plot_features
#' @export
plot_vdj_feature <- function(input, ...) {

  UseMethod("plot_vdj_feature", input)
}

#' @rdname plot_features
#' @description plot_vdj_feature() allows per-chain V(D)J data to be summarized
#' and plotted for each cell.
#' @param data_col meta.data column containing V(D)J data to use for coloring
#' cells
#' @param summary_fn Function to use for summarizing values for each cell,
#' possible values can be either a function, e.g. mean, or a purrr-style
#' lambda, e.g. ~ mean(.x, na.rm = TRUE) where ".x" refers to the column. If
#' NULL, the mean will be calculated for numeric values, non-numeric columns
#' will be combined into a single string.
#' @param chain Chain(s) to use for filtering data before plotting. If NULL
#' data will not be filtered based on chain.
#' @param chain_col meta.data column containing chains for each cell
#' @param sep Separator used for storing per-chain V(D)J data for each cell
#' @param ... Additional arguments to pass to ggplot2, e.g. color, fill, size,
#' linetype, etc.
#' @seealso [summarize_vdj()] for more examples on how per-chain data can be
#' summarized for each cell
#'
#' @examples
#' # Plot average CDR3 length for each cell for light chains
#' plot_vdj_feature(
#'   vdj_so,
#'   data_col = "cdr3_length",
#'   summary_fn = mean,
#'   chain = c("IGK", "IGL")
#' )
#'
#' # Plot median number of insertions for each cell
#' plot_vdj_feature(
#'   vdj_sce,
#'   data_col = "all_ins",
#'   summary_fn = stats::median
#' )
#'
#' # Using a lambda function to summarize values
#' # use ".x" to refer to values in the column
#' # this calculates the mean number of deletions for each cell and
#' # log10-transforms the result
#' plot_vdj_feature(
#'   vdj_so,
#'   data_col = "all_del",
#'   summary_fn = ~ log10(mean(.x) + 1)
#' )
#'
#' @export
plot_vdj_feature.default <- function(input, data_col, x = "UMAP_1",
                                     y = "UMAP_2", group_col = NULL,
                                     summary_fn = NULL,
                                     chain = NULL, plot_colors = NULL,
                                     plot_lvls = names(plot_colors),
                                     trans = "identity", min_q = NULL,
                                     max_q = NULL, panel_nrow = NULL,
                                     panel_scales = "fixed",
                                     na_color = "grey80", chain_col = "chains",
                                     sep = ";", ...) {

  plt_dat <- summarize_vdj(
    input,
    data_cols = data_col,
    fn        = summary_fn,
    chain     = chain,
    chain_col = chain_col,
    sep       = sep,
    return_df = TRUE
  )

  res <- plot_features(
    plt_dat,
    x            = x,
    y            = y,
    feature      = data_col,
    group_col    = group_col,
    plot_colors  = plot_colors,
    plot_lvls    = plot_lvls,
    min_q        = min_q,
    max_q        = max_q,
    panel_nrow   = panel_nrow,
    panel_scales = panel_scales,
    na_color     = na_color,
    trans        = trans,
    ...
  )

  res
}

#' @rdname plot_features
#' @param data_slot Slot in the Seurat object to pull data
#' @importFrom Seurat FetchData
#' @export
plot_vdj_feature.Seurat <- function(input, data_col, x = "UMAP_1",
                                    y = "UMAP_2", group_col = NULL,
                                    data_slot = "data", summary_fn = NULL,
                                    chain = NULL, plot_colors = NULL,
                                    plot_lvls = names(plot_colors),
                                    trans = "identity", min_q = NULL,
                                    max_q = NULL, panel_nrow = NULL,
                                    panel_scales = "fixed", na_color = "grey80",
                                    chain_col = "chains", sep = ";", ...) {

  # Fetch variables and add to meta.data
  # want input data to include meta.data and any features from FetchData
  # SHOULD WARNINGS BE SUPPRESSED??
  plt_vars <- c(x, y, data_col, group_col)

  plt_dat <- Seurat::FetchData(
    input,
    vars = unique(plt_vars),
    slot = data_slot
  )

  # Format input data
  # keep rownames since default method will create rowname column
  plt_dat <- .merge_meta(input, plt_dat)
  plt_dat <- .get_meta(plt_dat)

  res <- plot_vdj_feature(
    plt_dat,
    x            = x,
    y            = y,
    data_col     = data_col,
    group_col    = group_col,
    summary_fn   = summary_fn,
    chain        = chain,
    plot_colors  = plot_colors,
    plot_lvls    = plot_lvls,
    min_q        = min_q,
    max_q        = max_q,
    panel_nrow   = panel_nrow,
    panel_scales = panel_scales,
    na_color     = na_color,
    trans        = trans,
    chain_col    = chain_col,
    sep          = sep,
    ...
  )

  res
}


#' Plot continuous V(D)J data
#'
#' @param input Single cell object or data.frame containing V(D)J data. If a
#' data.frame is provided, cell barcodes should be stored as row names.
#' @param data_col meta.data column(s) containing continuous V(D)J data to
#' plot
#' @param cluster_col meta.data column containing cluster IDs to use for
#' grouping cells for plotting
#' @param group_col meta.dats column to use for grouping clusters present in
#' cluster_col
#' @param chain Chain(s) to use for filtering data before plotting. If NULL
#' data will not be filtered based on chain.
#' @param method Method to use for plotting, possible values include:
#'
#' - 'histogram'
#' - 'density'
#' - 'boxplot'
#' - 'violin'
#'
#' @param units Units to use for y-axis when method is set to 'histogram'. Use
#' 'frequency' to show number of values or 'percent' to show the percentage of
#' total values.
#'
#' @param plot_colors Character vector specifying colors to use for cell
#' clusters specified by cluster_col. When cluster_col is NULL, plot colors can
#' be directly modified with the ggplot2 parameters color and fill,
#' e.g. fill = "red", color = "black"
#' @param plot_lvls Character vector containing order to use for plotting cell
#' clusters specified by cluster_col
#' @param panel_nrow The number of rows to use for arranging plot panels
#' @param panel_scales Should scales for plot panels be fixed or free. This
#' passes a scales specification to ggplot2::facet_wrap, can be 'fixed', 'free',
#' 'free_x', or 'free_y'. 'fixed' will cause panels to share the same scales.
#' @param chain_col meta.data column containing chains for each cell
#' @param n_label Include a label showing the number of cells plotted
#' @param label_params Named list providing additional parameters to modify
#' clonotype and n label aesthetics, e.g. list(size = 4, color = "red")
#' @param sep Separator used for storing per-chain V(D)J data for each cell
#' @param ... Additional arguments to pass to ggplot2, e.g. color, fill, size,
#' linetype, etc.
#' @return ggplot object
#' @name plot_vdj
NULL

#' @rdname plot_vdj
#' @param per_cell Should values be plotted per cell, i.e. each data point
#' would represent one cell. If TRUE, values will be summarized for each cell
#' using summary_fn. If FALSE, values will be plotted for each chain.
#' @param summary_fn Function to use for summarizing values when per_cell is
#' TRUE, possible values can be either a function, e.g. mean, or a purrr-style
#' lambda, e.g. ~ mean(.x, na.rm = TRUE) where '.x' refers to the column. If
#' NULL, the mean will be calculated.
#' @param trans Transformation to use for plotting data, e.g. 'log10'. By
#' default values are not transformed, refer to [ggplot2::continuous_scale()]
#' for more options.
#' @seealso [summarize_vdj()] for more examples on how per-chain data can be
#' summarized for each cell
#'
#' @examples
#' # Create histogram
#' plot_vdj(
#'   vdj_sce,
#'   data_col = "reads"
#' )
#'
#' # Create boxplots
#' plot_vdj(
#'   vdj_sce,
#'   data_col = "reads",
#'   method = "boxplot"
#' )
#'
#' # Pass additional arguments to ggplot2
#' plot_vdj(
#'   vdj_so,
#'   data_col = "reads",
#'   color = "red",
#'   bins = 25
#' )
#'
#' # Compare cell clusters
#' plot_vdj(
#'   vdj_sce,
#'   data_col = "cdr3_length",
#'   cluster_col = "orig.ident",
#'   method = "violin"
#' )
#'
#' # log10 transform the axis
#' plot_vdj(
#'   vdj_so,
#'   data_col = "reads",
#'   cluster_col = "orig.ident",
#'   trans = "log10"
#' )
#'
#' # Express y-axis units as percent of total values
#' plot_vdj(
#'   vdj_sce,
#'   data_col = "umis",
#'   cluster_col = "orig.ident",
#'   units = "percent"
#' )
#'
#' # Only plot values for heavy chains
#' plot_vdj(
#'   vdj_so,
#'   data_col = "umis",
#'   chain = "IGH"
#' )
#'
#' # Plot the median number of reads for each cell
#' plot_vdj(
#'   vdj_sce,
#'   data_col = "reads",
#'   per_cell = TRUE,
#'   summary_fn = stats::median
#' )
#'
#' # Set colors for cell clusters
#' plot_vdj(
#'   vdj_so,
#'   data_col = "cdr3_length",
#'   cluster_col = "orig.ident",
#'   plot_colors = c(avid_1 = "red", avid_2 = "purple")
#' )
#'
#' # Set order to use for plotting cell clusters
#' plot_vdj(
#'   vdj_sce,
#'   data_col = "cdr3_length",
#'   cluster_col = "orig.ident",
#'   plot_lvls = c("avid_2", "avid_1"),
#'   method = "boxplot"
#' )
#'
#' @export
plot_vdj <- function(input, data_col, per_cell = FALSE, summary_fn = mean,
                     cluster_col = NULL, group_col = NULL, chain = NULL,
                     method = "histogram", units = "frequency",
                     plot_colors = NULL, plot_lvls = names(plot_colors),
                     trans = "identity", panel_nrow = NULL,
                     panel_scales = "free_x", chain_col = "chains",
                     n_label = TRUE, label_params = list(), sep = ";", ...) {

  # Format input data
  if (per_cell) {
    plt_dat <- summarize_vdj(
      input,
      data_cols = data_col,
      fn        = summary_fn,
      chain     = chain,
      chain_col = chain_col,
      sep       = sep,
      return_df = TRUE
    )

  } else {
    fetch_cols <- data_col

    if (!is.null(chain)) fetch_cols <- c(chain_col, fetch_cols)

    plt_dat <- fetch_vdj(
      input,
      data_cols     = fetch_cols,
      clonotype_col = NULL,
      unnest        = TRUE
    )

    if (!is.null(chain)) {
      plt_dat <- dplyr::filter(plt_dat, !!sym(chain_col) %in% chain)
    }
  }

  plt_dat <- dplyr::filter(plt_dat, !is.na(!!sym(data_col)))

  .n <- nrow(plt_dat)

  # Create plots
  gg_args <- list(
    df_in        = plt_dat,
    data_col     = data_col,
    cluster_col  = cluster_col,
    group_col    = group_col,
    method       = method,
    units        = units,
    plot_colors  = plot_colors,
    plot_lvls    = plot_lvls,
    trans        = trans,
    panel_nrow   = panel_nrow,
    panel_scales = panel_scales,
    per_cell     = per_cell,
    ...
  )

  if (n_label) gg_args$y_exp <- c(0.05, 0.1)

  res <- purrr::lift_dl(.plot_vdj)(gg_args)

  # Add n label
  if (n_label) res <- .add_n_label(res, label_params, .n)

  res
}

#' @rdname plot_vdj
#' @param df_in input data.frame
#' @noRd
.plot_vdj <- function(df_in, data_col, cluster_col = NULL, group_col = NULL,
                      method = "boxplot", units = "frequency",
                      plot_colors = NULL, plot_lvls = names(plot_colors),
                      trans = "identity", y_exp = c(0.05, 0.05),
                      panel_nrow = NULL, panel_scales = "fixed",
                      per_cell = FALSE, ...) {

  # Order clusters based on plot_lvls
  df_in <- .set_lvls(df_in, cluster_col, plot_lvls)

  # Set ggplot arguments
  gg_args <- list(
    df_in  = df_in,
    x      = data_col,
    grp    = group_col,
    .color = cluster_col,
    .fill  = cluster_col,
    clrs   = plot_colors,
    method = method,
    trans  = trans,
    y_exp  = y_exp,
    nrow   = panel_nrow,
    scales = panel_scales,
    ...
  )

  # Set default alpha to 0.5
  gg_args$alpha <- gg_args$alpha %||% 0.5

  # Create histogram
  if (method %in% c("histogram", "density")) {
    gg_args$units <- units

    gg_args$y_ttl <- .get_axis_label(
      units, sfx = ifelse(per_cell, "cells", "chains")
    )

    res <- purrr::lift_dl(.create_hist)(gg_args)

    return(res)
  }

  # Create boxplot plot
  gg_args$x <- cluster_col
  gg_args$y <- data_col

  res <- purrr::lift_dl(.create_boxes)(gg_args)

  res
}
