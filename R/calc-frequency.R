#' Calculate frequency of cell groups present in object
#'
#' Calculate the frequency of each cell label present in the provided meta.data
#' column. This is useful for comparing the proportion of cells belonging to
#' different samples, cell types, clonotypes, isotypes, etc.
#'
#' @param input Single cell object or data.frame containing V(D)J data. If a
#' data.frame is provided, the cell barcodes should be stored as row names.
#' @param data_col meta.data column containing cell labels to use for
#' calculating frequency. To calculate clonotype frequencies, provide the column
#' containing clonotype IDs, to calculate isotype frequencies provide the column
#' containing cell isotypes. By default the clonotype_id is used for
#' calculations.
#' @param cluster_col meta.data column containing cluster IDs to use for
#' grouping cells when calculating clonotype abundance
#' @param prefix Prefix to add to new columns
#' @param return_df Return results as a data.frame. If set to FALSE, results
#' will be added to the input object.
#' @return Single cell object or data.frame with clonotype frequencies
#' @seealso [plot_frequency()], [plot_clone_frequency()]
#'
#' @examples
#' # Calculate clonotype abundance using all cells
#' res <- calc_frequency(
#'   vdj_so,
#'   data_col = "clonotype_id"
#' )
#'
#' head(slot(res, "meta.data"), 1)
#'
#' # Group cells based on meta.data column before calculating abundance
#' res <- calc_frequency(
#'   vdj_sce,
#'   data_col = "clonotype_id",
#'   cluster_col = "orig.ident"
#' )
#'
#' head(slot(res, "colData"), 1)
#'
#' # Add a prefix to the new columns
#' # this is useful if multiple abundance calculations are stored in the
#' # meta.data
#' res <- calc_frequency(
#'   vdj_so,
#'   data_col = "clonotype_id",
#'   prefix = "bcr_"
#' )
#'
#' head(slot(res, "meta.data"), 1)
#'
#' # Return a data.frame instead of adding the results to the input object
#' res <- calc_frequency(
#'   vdj_sce,
#'   data_col = "clonotype_id",
#'   return_df = TRUE
#' )
#'
#' head(res, 1)
#'
#' @export
calc_frequency <- function(input, data_col, cluster_col = NULL, prefix = paste0
                           (data_col, "_"), return_df = FALSE) {

  # Check that columns are present in object
  .check_obj_cols(input, data_col, cluster_col)

  # Check arguments
  .check_args()

  # Format input data
  meta <- .get_meta(input)
  vdj  <- dplyr::filter(meta, !is.na(!!sym(data_col)))

  vdj <- dplyr::select(vdj, all_of(c(global$cell_col, data_col, cluster_col)))

  # Calculate clonotype abundance
  vdj <- .calc_freq(
    df_in      = vdj,
    cell_col   = global$cell_col,
    dat_col    = data_col,
    clust_col  = cluster_col,
    out_prefix = prefix
  )

  vdj <- dplyr::select(vdj, -all_of(data_col))

  freq_clmn <- paste0(prefix, "freq")
  grp_clmn  <- paste0(prefix, "grp")

  vdj <- dplyr::mutate(
    vdj,
    !!sym(grp_clmn) := .calc_freq_grp(!!sym(freq_clmn))
  )

  # Format results
  res <- dplyr::left_join(meta, vdj, by = global$cell_col)

  if (return_df) {
    input <- meta
  }

  res <- .add_meta(input, meta = res)

  res
}

#' Calculate frequency of a cell label
#'
#' @param df_in Input data.frame
#' @param cell_col Column containing cell IDs
#' @param data_col Column containing data for calculating abundance
#' (e.g. clonotype IDs)
#' @param clust_col Column containing cluster IDs to use for grouping cells
#' @param out_prefix Prefix to add to output columns
#' @return data.frame containing clonotype abundances
#' @noRd
.calc_freq <- function(df_in, cell_col, dat_col, clust_col = NULL,
                       out_prefix = "") {

  # Count number of cells in each group
  if (!is.null(clust_col)) {
    df_in <- dplyr::group_by(df_in, !!sym(clust_col))
  }

  df_in <- dplyr::mutate(
    df_in,
    .n_cells = dplyr::n_distinct(!!sym(cell_col))
  )

  # Calculate frequency
  res <- dplyr::group_by(df_in, !!sym(dat_col), .add = TRUE)

  res <- dplyr::mutate(
    res,
    .freq = dplyr::n_distinct(!!sym(cell_col)),
    .pct  = (.data$.freq / .data$.n_cells) * 100
  )

  # Identify shared labels
  if (!is.null(clust_col)) {
    res <- dplyr::group_by(res, !!sym(dat_col))

    res <- dplyr::mutate(
      res,
      .shared = dplyr::n_distinct(!!sym(clust_col)) > 1
    )
  }

  res <- dplyr::ungroup(res)

  # Rename output columns
  new_cols <- c("freq", "pct")

  if (!is.null(clust_col)) {
    new_cols <- c(new_cols, "shared")
  }

  new_cols <- purrr::set_names(
    paste0(".", new_cols),
    paste0(out_prefix, new_cols)
  )

  res <- dplyr::select(
    res,
    all_of(c(cell_col, dat_col)),
    !!!syms(new_cols)
  )

  res
}

#' Divide clonotypes into groups based on frequency
#'
#' @param x Numeric vector of clonotype frequencies
#' @param n_grps Number of groups to return
#' @return Named vector containing group labels
#' @noRd
.calc_freq_grp <- function(x, n_grps = 4) {

  n_grps <- n_grps - 1

  # Set unique clonotypes as their own group
  uniq_x <- unique(x)
  uniq_x <- uniq_x[uniq_x > 1]
  uniq_x <- sort(uniq_x)

  # Divide into groups
  labs <- tibble::tibble(
    x   = uniq_x,
    grp = dplyr::ntile(uniq_x, n_grps)
  )

  # Format labels
  labs <- dplyr::group_by(labs, .data$grp)

  labs <- dplyr::mutate(
    labs,
    lab = paste0(unique(range(x)), collapse = "-")
  )

  labs <- purrr::set_names(labs$lab, labs$x)
  labs <- c("1" = "1", labs)

  res <- unname(labs[as.character(x)])

  lvls <- c(NA, unique(unname(labs)))

  res <- factor(res, levels = lvls, exclude = NULL)

  res
}


#' Fetch top clonotypes
#'
#' WIP
#'
#' @param input Single cell object or data.frame containing V(D)J data. If a
#' data.frame is provided, the cell barcodes should be stored as row names.
#' @param data_cols meta.data column(s) to use for identifying top clonotypes
#' @param cluster_col meta.data column containing cluster IDs to use for
#' grouping cells.
#' @param n_clones Number of top clonotypes to identify.
#' @param sep Separator used for storing per-chain V(D)J data for each cell
#' @return data.frame containing information for top clonotypes
#' @noRd
fetch_top_clones <- function(input, data_cols, cluster_col = NULL,
                             n_clones = 10, sep = global$sep) {

  # Check that columns are present in object
  .check_obj_cols(input, data_cols, cluster_col)

  # Check input classes
  .check_args()

}


#' Plot clonotype frequency
#'
#' @param input Single cell object or data.frame containing V(D)J data. If a
#' data.frame is provided, the cell barcodes should be stored as row names.
#' @param data_col meta.data column containing clonotype IDs to use for
#' calculating clonotype abundance
#' @param cluster_col meta.data column containing cluster IDs to use for
#' grouping cells when calculating clonotype frequencies. Clonotypes will be
#' plotted separately for each cluster.
#' @param group_col meta.data column to use for grouping cluster IDs present in
#' cluster_col. This is useful when there are multiple replicates or patients
#' for each treatment condition. This is only applicable when `method` is
#' 'line'.
#' @param method Method to use for plotting, possible values include:
#'
#' - 'bar', create a bargraph
#' - 'line', create a rank-abundance plot
#'
#' @param units Units to plot on the y-axis, either 'frequency' or 'percent'
#' @param plot_colors Character vector containing colors for plotting
#' @param plot_lvls Levels to use for ordering clusters
#' @param trans Transformation to use for plotting data, e.g. 'log10'. By
#' default values are not transformed, refer to [ggplot2::continuous_scale()]
#' @param n_clones Number of top clonotypes to plot (default is 10). If method
#' is set to 'line', this will specify the number of clonotypes to label
#' (default is 3).
#' @param panel_nrow The number of rows to use for arranging plot panels, use
#' this when separate bar graphs are created for each cell cluster
#' @param panel_scales Should scales for plot panels be fixed or free. This
#' passes a scales specification to ggplot2::facet_wrap, can be 'fixed', 'free',
#' 'free_x', or 'free_y'. 'fixed' will cause panels to share the same scales.
#' Use this when separate bar graphs are created for each cell cluster.
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
#' clonotype and n label aesthetics, e.g. list(size = 4, color = "red")
#' @param ... Additional arguments to pass to ggplot2, e.g. color, fill, size,
#' linetype, etc.
#' @return ggplot object
#' @importFrom ggrepel geom_text_repel
#' @importFrom methods formalArgs
#' @seealso [calc_frequency()], [plot_frequency()]
#'
#' @examples
#' # Plot clonotype abundance using all cells
#' plot_clone_frequency(
#'   vdj_so,
#'   data_col = "clonotype_id"
#' )
#'
#' # Plot clonotype abundance separately for each cell cluster
#' plot_clone_frequency(
#'   vdj_sce,
#'   cluster_col = "orig.ident"
#' )
#'
#' # Plot the frequency of each clonotype instead of percentage
#' plot_clone_frequency(
#'   vdj_sce,
#'   cluster_col = "orig.ident",
#'   units = "frequency"
#' )
#'
#' # Specify colors to use for each cell cluster
#' plot_clone_frequency(
#'   vdj_so,
#'   cluster_col = "orig.ident",
#'   plot_colors = c(avid_1 = "blue", avid_2 = "red")
#' )
#'
#' # Specify order to use for plotting cell clusters
#' plot_clone_frequency(
#'   vdj_sce,
#'   cluster_col = "orig.ident",
#'   plot_lvls = c("avid_2", "avid_1")
#' )
#'
#' # Specify the number of top clonotypes to plot
#' plot_clone_frequency(
#'   vdj_so,
#'   n_clones = 5
#' )
#'
#' #' # Create line graph
#' # use n_clones to set the number of clonotypes to label
#' plot_clone_frequency(
#'   vdj_so,
#'   cluster_col = "orig.ident",
#'   method = "line",
#'   n_clones = 3
#' )
#'
#' @export
plot_clone_frequency <- function(input, data_col = global$clonotype_col,
                                 cluster_col = NULL, group_col = NULL,
                                 method = "bar", units = "percent",
                                 plot_colors = NULL,
                                 plot_lvls = names(plot_colors),
                                 trans = "identity", n_clones = NULL,
                                 panel_nrow = NULL, panel_scales = "free_x",
                                 n_label = "corner", label_params = list(), ...) {

  # Check that columns are present in object
  .check_obj_cols(input, data_col, cluster_col, group_col)

  # Check input classes
  .check_args()

  # Check input values
  .check_possible_values(
    method = c("bar", "line"),
    units  = c("frequency", "percent")
  )

  n_clones <- n_clones %||% switch(method, bar = 10, line = 3)

  if (identical(method, "bar") && n_clones <= 0) {
    cli::cli_abort("`n_clones` must be >0")
  }

  if (identical(method, "line") && n_clones < 0) {
    cli::cli_abort("`n_clones` must be >=0")
  }

  # For bargraph allow user to facet plot using group_col when cluster_col is
  # NULL
  if (identical(method, "bar") && is.null(cluster_col)) {
    cluster_col <- group_col
    group_col   <- NULL
  }

  .check_group_cols(cluster_col, group_col, input)

  abun_col <- switch(units, frequency = ".freq", percent = ".pct")
  y_lab    <- .get_axis_label(units)

  # Calculate clonotype abundance
  plt_dat <- calc_frequency(
    input       = input,
    cluster_col = cluster_col,
    data_col    = data_col,
    prefix      = ".",
    return_df   = TRUE
  )

  plt_dat <- tibble::as_tibble(plt_dat, rownames = global$cell_col)
  plt_dat <- dplyr::filter(plt_dat, !is.na(!!sym(data_col)))

  # Save data for n label
  n_lab_dat  <- plt_dat

  # Identify data columns that the user should have access to
  keep_cols <- .get_matching_clmns(plt_dat, c(data_col, cluster_col))
  keep_cols <- c(cluster_col, data_col, keep_cols)
  plt_dat   <- dplyr::distinct(plt_dat, !!!syms(keep_cols))

  # Add and format label column for plotting
  plt_dat <- dplyr::mutate(plt_dat, .lab = trim_lab(!!sym(data_col)))

  # Rank by abundance
  if (!is.null(cluster_col)) {
    plt_dat <- .set_lvls(plt_dat, cluster_col, plot_lvls)
    plt_dat <- dplyr::group_by(plt_dat, !!sym(cluster_col))

    plt_dat <- dplyr::mutate(
      plt_dat,
      !!sym(data_col) := paste0(!!sym(cluster_col), "_", !!sym(data_col))
    )
  }

  plt_dat <- dplyr::mutate(
    plt_dat, rank = dplyr::row_number(dplyr::desc(!!sym(abun_col)))
  )

  # Identify top clonotypes
  top_clones <- dplyr::slice_min(
    plt_dat,
    order_by  = rank,
    n         = n_clones,
    with_ties = FALSE
  )

  plt_dat    <- dplyr::ungroup(plt_dat)
  top_clones <- dplyr::ungroup(top_clones)

  # Create bar graph
  gg_args <- list(
    grp          = cluster_col,
    clrs         = plot_colors,
    nrow         = panel_nrow,
    scales       = panel_scales,
    n_label      = n_label,
    label_params = label_params,
    label_data   = n_lab_dat,
    ...
  )

  if (identical(method, "bar")) {
    plt_labs   <- purrr::set_names(top_clones$.lab, top_clones[[data_col]])
    top_clones <- dplyr::arrange(top_clones, desc(!!sym(abun_col)))
    lvls       <- unique(top_clones[[data_col]])

    top_clones <- .set_lvls(top_clones, data_col, lvls)

    gg_args$df_in <- top_clones
    gg_args$x     <- data_col
    gg_args$y     <- abun_col
    gg_args$y_ttl <- y_lab
    gg_args$.fill <- gg_args$.color <- cluster_col
    gg_args$ang   <- 45
    gg_args$hjst  <- 1
    gg_args$trans <- trans

    res <- lift(.create_bars)(gg_args)

    # Format clonotype labels
    res <- res +
      ggplot2::scale_x_discrete(labels = plt_labs)

  } else {
    gg_args$df_in     <- plt_dat
    gg_args$fn        <- ggplot2::geom_line
    gg_args$x         <- "rank"
    gg_args$y         <- abun_col
    gg_args$.color    <- cluster_col
    gg_args$grp       <- group_col
    gg_args$trans_y   <- trans
    gg_args$linewidth <- gg_args$linewidth %||% 1

    if (n_clones > 0) {
      gg_args$label_params <- .get_uniq_text_args(label_params, "geom_text")
    }

    res <- lift(.create_plot)(gg_args) +
      labs(y = .get_axis_label(units))

    # Add clonotype labels
    lab_args <- label_params

    if (n_clones > 0) {
      if (!all(n_label == "none")) {
        label_params <- .get_uniq_text_args(label_params, "geom_text_repel")
      }

      label_params$mapping        <- ggplot2::aes(label = .data$.lab)
      label_params$data           <- top_clones
      label_params$nudge_x        <- label_params$nudge_x %||% 500
      label_params$direction      <- label_params$direction %||% "y"
      label_params$segment.colour <- label_params$segment.colour %||% "black"
      label_params$segment.size   <- label_params$segment.size %||% 0.2
      label_params$segment.alpha  <- label_params$segment.alpha %||% 0.2
      label_params$show.legend    <- label_params$show.legend %||% FALSE

      if (!is.null(label_params$size)) {
        label_params$size <- label_params$size / ggplot2::.pt
      }

      res <- res +
        lift(ggrepel::geom_text_repel)(label_params)
    }
  }

  res
}


#' Plot frequency of cell groups present in object
#'
#' Plot the frequency of each cell label present in the provided meta.data
#' column. This is useful for comparing the proportion of cells belonging to
#' different samples, cell types, clonotypes, isotypes, etc.
#'
#' @param input Single cell object or data.frame containing V(D)J data. If a
#' data.frame is provided, the cell barcodes should be stored as row names.
#' @param data_col meta.data column containing cell labels to use for
#' calculating frequency, e.g. cell types, isotypes etc.
#' @param cluster_col meta.data column containing cluster IDs (or patients,
#' treatment conditions, etc.) to use when calculating frequency. Calculations
#' will be performed separately for each cluster.
#' @param group_col meta.data column to use for grouping cluster IDs present in
#' cluster_col. This is useful when there are multiple replicates or patients
#' for each treatment condition.
#' @param units Units to plot on the y-axis, either 'frequency' or 'percent'
#' @param stack If TRUE, stacked bargraphs will be generated, otherwise grouped
#' bargraphs will be generated
#' @param plot_colors Character vector containing colors for plotting
#' @param plot_lvls Levels to use for ordering clusters or groups
#' @param trans Transformation to use for plotting data, e.g. 'log10'. By
#' default values are not transformed, refer to [ggplot2::continuous_scale()]
#' for more options. Values can only be transformed when stack is `FALSE`
#' @param n_top Number of top cell labels present in data_col to show on plot,
#' other cells will be labeled based on the `other_label` argument. If `NULL`,
#' this will be automatically selected.
#' @param other_label Label to use for 'other' cells, if `NULL` all cell labels
#' present in data_col will be displayed on the plot.
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
#' @return ggplot object
#' @seealso [calc_frequency()], [plot_clone_frequency()]
#'
#' @examples
#' # Plot frequency of different isotypes
#' plot_frequency(
#'   vdj_so,
#'   data_col = "isotype"
#' )
#'
#' # Plot frequency separately for cell clusters
#' plot_frequency(
#'   vdj_sce,
#'   data_col    = "isotype",
#'   cluster_col = "orig.ident"
#' )
#'
#' # Create grouped bar graphs
#' plot_frequency(
#'   vdj_sce,
#'   data_col    = "isotype",
#'   cluster_col = "orig.ident",
#'   stack       = FALSE
#' )
#'
#' # Plot number of cells on the y-axis
#' plot_frequency(
#'   vdj_so,
#'   data_col    = "seurat_clusters",
#'   cluster_col = "orig.ident",
#'   units       = "frequency"
#' )
#'
#' @export
plot_frequency <- function(input, data_col, cluster_col = NULL,
                           group_col = NULL, units = "percent", stack = NULL,
                           plot_colors = NULL, plot_lvls = NULL,
                           trans = "identity", n_top = NULL,
                           other_label = "other", n_label = "axis",
                           label_params = list(), ...) {

  # Check that columns are present in object
  .check_obj_cols(input, data_col, cluster_col)

  # Check input classes
  .check_args()

  # Check input values
  .check_group_cols(cluster_col, group_col, input)

  .check_possible_values(units = c("frequency", "percent"))

  if (!identical(trans, "identity")) stack <- stack %||% FALSE
  else                               stack <- stack %||% !is.null(cluster_col)

  if (stack && !identical(trans, "identity")) {
    cli::cli_abort("Values can only be transformed when `stack` is `FALSE`")
  }

  abun_col <- switch(units, frequency = ".freq", percent = ".pct")

  y_lab <- .get_axis_label(units)

  # Set plot variables
  if (!is.null(group_col))        x_col <- data_col
  else if (!is.null(cluster_col)) x_col <- cluster_col
  else                            x_col <- data_col

  clr_col <- group_col %||% data_col

  # Calculate clonotype abundance
  plt_dat <- calc_frequency(
    input       = input,
    cluster_col = cluster_col,
    data_col    = data_col,
    prefix      = ".",
    return_df   = TRUE
  )

  # Format plot data
  plt_dat <- tibble::as_tibble(plt_dat, rownames = global$cell_col)
  plt_dat <- dplyr::filter(plt_dat, !is.na(!!sym(data_col)))

  # Set data for n label
  n_lab_dat  <- plt_dat

  # Identify data columns that the user should have access to
  keep_cols <- .get_matching_clmns(plt_dat, c(data_col, cluster_col))
  keep_cols <- c(cluster_col, data_col, keep_cols)
  plt_dat   <- dplyr::distinct(plt_dat, !!!syms(keep_cols))

  # Rank values in data_col
  .rank_values <- function(df_in, dat_clmn, val_clmn) {
    res <- dplyr::group_by(df_in, !!sym(dat_clmn))
    res <- dplyr::summarize(res, mn = mean(!!sym(val_clmn)))
    res <- dplyr::arrange(res, desc(.data$mn))
    res <- pull(res, dat_clmn)
    res
  }

  rnk <- .rank_values(plt_dat, data_col, abun_col)

  # Set other group based on top groups in data_col
  n_dat <- length(rnk)
  n_top <- n_top %||% ifelse(n_dat > 50, 10, 20)

  if (n_top < n_dat && !is.null(other_label)) {
    keep_dat <- rnk[seq_len(n_top)]

    plt_dat <- dplyr::mutate(plt_dat, !!sym(data_col) := ifelse(
      !!sym(data_col) %in% keep_dat,
      !!sym(data_col),
      other_label
    ))

    keep_cols <- .get_matching_clmns(plt_dat, c(data_col, cluster_col))
    keep_cols <- c(cluster_col, data_col, keep_cols)

    plt_dat <- dplyr::group_by(plt_dat, !!!syms(keep_cols))

    plt_dat <- dplyr::summarize(
      plt_dat, !!sym(abun_col) := sum(!!sym(abun_col))
    )

    rnk <- .rank_values(plt_dat, data_col, abun_col)

    if (!is.null(names(plot_colors)) && !other_label %in% names(plot_colors)) {
      plot_colors[[other_label]] <- "grey80"
    }
  }

  plt_dat <- .set_lvls(plt_dat, data_col, rnk)

  # Plot arguments
  gg_args <- list(
    y = abun_col, clrs = plot_colors, trans = trans,
    n_label = n_label, label_params = label_params,
    label_data = n_lab_dat, ...
  )

  # Create grouped boxplot
  if (!is.null(group_col)) {
    plot_lvls <- plot_lvls %||% names(plot_colors)
    plt_dat   <- .set_lvls(plt_dat, group_col, plot_lvls)

    gg_args$df_in         <- plt_dat
    gg_args$x             <- x_col
    gg_args$.fill         <- gg_args$.color <- clr_col
    gg_args$alpha         <- gg_args$alpha %||% 0.5
    gg_args$outlier.color <- gg_args$outlier.color %||% NA

    res <- lift(.create_boxes)(gg_args) +
      ggplot2::geom_jitter(
        position = ggplot2::position_jitterdodge(jitter.width = 0.05)
      ) +
      labs(y = y_lab) +
      theme(legend.position = "right")

  # Create bar graph
  } else {
    plt_dat <- .set_lvls(plt_dat, cluster_col, plot_lvls)

    gg_args$df_in <- plt_dat
    gg_args$x     <- x_col
    gg_args$y_ttl <- y_lab
    gg_args$.fill <- gg_args$.color <- clr_col
    gg_args$ang   <- 45
    gg_args$hjst  <- 1

    # When cluster_col is provided set default position to dodge
    if (!is.null(cluster_col)) {
      gg_pos <- ggplot2::position_dodge(preserve = "single")

      if (stack) gg_pos <- ggplot2::position_stack()

      gg_args$position <- gg_args$position %||% gg_pos
    }

    res <- lift(.create_bars)(gg_args)
  }

  res
}
