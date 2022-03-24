#' Calculate abundance of a cell label
#'
#' Count the frequency of the provided cell labels such as clonotypes or isotypes
#'
#' @param input Single cell object or data.frame containing V(D)J data. If a
#' data.frame is provided, the cell barcodes should be stored as row names.
#' @param data_col meta.data column containing cell labels to use for
#' calculating abundance. To calculate clonotype abundance, provide the column
#' containing clonotype IDs, to calculate isotype abundance provide the column
#' containing cell isotypes. By default the clonotype_id is used for
#' calculations.
#' @param cluster_col meta.data column containing cluster IDs to use for
#' grouping cells when calculating clonotype abundance
#' @param prefix Prefix to add to new columns
#' @param return_df Return results as a data.frame. If set to FALSE, results
#' will be added to the input object.
#' @return Single cell object or data.frame with clonotype abundance metrics
#'
#' @examples
#' # Calculate clonotype abundance using all cells
#' res <- calc_abundance(
#'   vdj_so,
#'   data_col = "clonotype_id"
#' )
#'
#' head(res@meta.data, 1)
#'
#' # Group cells based on meta.data column before calculating abundance
#' res <- calc_abundance(
#'   vdj_sce,
#'   cluster_col = "orig.ident"
#' )
#'
#' head(res@colData, 1)
#'
#' # Add a prefix to the new columns
#' # this is useful if multiple abundance calculations are stored in the
#' # meta.data
#' res <- calc_abundance(
#'   vdj_so,
#'   prefix = "bcr_"
#' )
#'
#' head(res@meta.data, 1)
#'
#' # Return a data.frame instead of adding the results to the input object
#' res <- calc_abundance(
#'   vdj_sce,
#'   return_df = TRUE
#' )
#'
#' head(res, 1)
#'
#' @export
calc_abundance <- function(input, data_col = "clonotype_id", cluster_col = NULL, prefix = "",
                           return_df = FALSE) {

  # Format input data
  meta <- .get_meta(input)
  vdj  <- dplyr::filter(meta, !is.na(!!sym(data_col)))

  vdj <- dplyr::select(
    vdj,
    !!sym(CELL_COL), all_of(c(data_col, cluster_col))
  )

  # Calculate clonotype abundance
  vdj <- .calc_abund(
    df_in      = vdj,
    cell_col   = CELL_COL,
    dat_col    = data_col,
    clust_col  = cluster_col,
    out_prefix = prefix
  )

  vdj <- dplyr::select(vdj, -all_of(data_col))

  # Format results
  res <- dplyr::left_join(meta, vdj, by = CELL_COL)

  if (return_df) {
    input <- meta
  }

  res <- .add_meta(input, meta = res)

  res
}

#' Calculate abundance of a cell label
#'
#' @param df_in Input data.frame
#' @param cell_col Column containing cell IDs
#' @param data_col Column containing data for calculating abundance
#' (e.g. clonotype IDs)
#' @param clust_col Column containing cluster IDs to use for grouping cells
#' @param out_prefix Prefix to add to output columns
#' @return data.frame containing clonotype abundances
#' @noRd
.calc_abund <- function(df_in, cell_col, dat_col, clust_col = NULL, out_prefix = "") {

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


#' Plot clonotype abundance
#'
#' @param input Single cell object or data.frame containing V(D)J data. If a
#' data.frame is provided, the cell barcodes should be stored as row names.
#' @param cluster_col meta.data column containing cluster IDs to use for
#' grouping cells when calculating clonotype abundance. Clonotypes will be
#' plotted separately for each cluster.
#' @param clonotype_col meta.data column containing clonotype IDs to use for
#' calculating clonotype abundance
#' @param type Type of plot to create, can be 'bar' or 'line'
#' @param yaxis Units to plot on the y-axis, either 'frequency' or 'percent'
#' @param plot_colors Character vector containing colors for plotting
#' @param plot_lvls Character vector containing levels for ordering
#' @param plot_n Number of top clonotypes to plot. If type is set to 'line',
#' this will specify the number of values to label.
#' @param label_aes Named list providing additional aesthetics (color, size,
#' etc.) for clonotype labels when creating line graph
#' @param facet_rows The number of facet rows, use this when separate bar
#' graphs are created for each cell cluster
#' @param facet_scales This passes a scales specification to
#' ggplot2::facet_wrap, can be 'fixed', 'free', 'free_x', or 'free_y'. 'fixed'
#' will cause plot facets to share the same scales. Use this when separate bar
#' graphs are created for each cell cluster.
#' @param ... Additional arguments to pass to ggplot2, e.g. color, fill, size,
#' linetype, etc.
#' @return ggplot object
#' @importFrom ggrepel geom_text_repel
#'
#' @examples
#' # Plot clonotype abundance using all cells
#' plot_abundance(
#'   vdj_so,
#'   data_col = "clonotype_id"
#' )
#'
#' # Plot clonotype abundance separately for each cell cluster
#' plot_abundance(
#'   vdj_sce,
#'   cluster_col = "orig.ident"
#' )
#'
#' # Plot the frequency of each clonotype instead of percentage
#' plot_abundance(
#'   vdj_sce,
#'   cluster_col = "orig.ident",
#'   yaxis = "frequency"
#' )
#'
#' # Specify colors to use for each cell cluster
#' plot_abundance(
#'   vdj_so,
#'   cluster_col = "orig.ident",
#'   plot_colors = c(avid_1 = "blue", avid_2 = "red")
#' )
#'
#' # Specify order to use for plotting cell clusters
#' plot_abundance(
#'   vdj_sce,
#'   cluster_col = "orig.ident",
#'   plot_lvls = c("avid_2", "avid_1")
#' )
#'
#' # Specify the number of top clonotypes to plot
#' plot_abundance(
#'   vdj_so,
#'   plot_n = 5
#' )
#'
#' #' # Create line graph
#' # use plot_n to set the number of clonotypes to label
#' plot_abundance(
#'   vdj_so,
#'   cluster_col = "orig.ident",
#'   type = "line",
#'   plot_n = 3
#' )
#'
#' @export
plot_clone_abundance <- function(input, cluster_col = NULL, clonotype_col = "clonotype_id", type = "bar",
                                 yaxis = "percent", plot_colors = NULL, plot_lvls = NULL, plot_n = 10,
                                 label_aes = list(), facet_rows = 1, facet_scales = "free_x", ...) {

  if (!yaxis %in% c("frequency", "percent")) {
    stop("yaxis must be either 'frequency' or 'percent'.")
  }

  if (!type %in% c("bar", "line")) {
    stop("type must be either 'bar' or 'line'.")
  }

  if (identical(type, "bar") && plot_n <= 0) {
    stop("If type is set to 'bar', plot_n must be >0.")
  }

  # Calculate clonotype abundance
  plt_dat <- calc_abundance(
    input       = input,
    cluster_col = cluster_col,
    data_col    = clonotype_col,
    prefix      = ".",
    return_df   = TRUE
  )

  abun_col <- ".pct"

  if (identical(yaxis, "frequency")) {
    abun_col <- ".freq"
  }

  plt_dat <- tibble::as_tibble(plt_dat, rownames = CELL_COL)
  plt_dat <- dplyr::filter(plt_dat, !is.na(!!sym(clonotype_col)))

  keep_cols <- c(cluster_col, clonotype_col, abun_col)

  plt_dat <- dplyr::distinct(plt_dat, !!!syms(keep_cols))

  # Add and format label column for plotting
  plt_dat <- dplyr::mutate(
    plt_dat,
    .lab = trim_lab(!!sym(clonotype_col))
  )

  # Rank by abundance
  if (!is.null(cluster_col)) {
    plt_dat <- .set_lvls(plt_dat, cluster_col, plot_lvls)
    plt_dat <- dplyr::group_by(plt_dat, !!sym(cluster_col))

    plt_dat <- dplyr::mutate(
      plt_dat,
      !!sym(clonotype_col) := paste0(!!sym(cluster_col), "_", !!sym(clonotype_col))
    )
  }

  plt_dat <- dplyr::mutate(
    plt_dat,
    rank = dplyr::row_number(dplyr::desc(!!sym(abun_col)))
  )

  # Identify top clonotypes
  top_clones <- dplyr::slice_min(
    plt_dat,
    order_by  = rank,
    n         = plot_n,
    with_ties = FALSE
  )

  plt_dat    <- dplyr::ungroup(plt_dat)
  top_clones <- dplyr::ungroup(top_clones)

  # Create bar graph
  if (identical(type, "bar")) {
    plt_labs <- purrr::set_names(
      top_clones$.lab,
      top_clones[[clonotype_col]]
    )

    top_clones <- dplyr::arrange(top_clones, desc(!!sym(abun_col)))

    lvls <- rev(unique(top_clones[[clonotype_col]]))

    top_clones <- .set_lvls(
      df_in = top_clones,
      clmn  = clonotype_col,
      lvls  = lvls
    )

    res <- .create_bars(
      df_in = top_clones,
      x     = clonotype_col,
      y     = abun_col,
      y_ttl = yaxis,
      .fill = cluster_col,
      clrs  = plot_colors,
      ang   = 45,
      hjst  = 1,
      ...
    )

    # Format clonotype labels
    res <- res +
      ggplot2::scale_x_discrete(labels = plt_labs)

    if (!is.null(cluster_col)) {
      res <- res +
        ggplot2::facet_wrap(
          stats::as.formula(paste0("~ ", cluster_col)),
          nrow   = facet_rows,
          scales = facet_scales
        )
    }

    return(res)
  }

  # Plot abundance vs rank
  plt_aes <- ggplot2::aes(rank, !!sym(abun_col))

  clr_aes <- plt_aes

  if (!is.null(cluster_col)) {
    clr_aes$colour <- sym(cluster_col)
  }

  res <- ggplot2::ggplot(plt_dat, plt_aes) +
    ggplot2::geom_line(clr_aes, ...) +
    ggplot2::labs(y = yaxis) +
    djvdj_theme()

  if (!is.null(plot_colors)) {
    res <- res +
      ggplot2::scale_color_manual(values = plot_colors)
  }

  # Add labels
  if (plot_n > 0) {
    res <- res +
      ggrepel::geom_text_repel(
        ggplot2::aes(label = .data$.lab),
        data          = top_clones,
        nudge_x       = 500,
        direction     = "y",
        segment.size  = 0.2,
        segment.alpha = 0.2,
        size          = 3
      )

    res <- .add_aes(res, label_aes, 2)
  }

  res
}

#' Plot isotype abundance
#'
#' @param input Single cell object or data.frame containing V(D)J data. If a
#' data.frame is provided, the cell barcodes should be stored as row names.
#' @param data_col meta.data column containing cell labels to use for
#' calculating abundance. To calculate clonotype abundance, provide the column
#' containing clonotype IDs, to calculate isotype abundance provide the column
#' containing cell isotypes. By default clonotype_id is used for calculations.
#' @param cluster_col meta.data column containing cluster IDs to use for
#' grouping cells when calculating clonotype abundance. Clonotypes will be
#' plotted separately for each cluster.
#' @param yaxis Units to plot on the y-axis, either 'frequency' or 'percent'
#' @param plot_colors Character vector containing colors for plotting
#' @param plot_lvls Character vector containing levels for ordering
#' @param facet_rows The number of facet rows, use this when separate bar
#' graphs are created for each cell cluster
#' @param ... Additional arguments to pass to ggplot2, e.g. color, fill, size,
#' linetype, etc.
#' @return ggplot object
#' @export
plot_abundance <- function(input, data_col = "clonotype_id", cluster_col = NULL, group_col = NULL,
                           yaxis = "percent", plot_colors = NULL, plot_lvls = NULL, facet_rows = NULL,
                           ...) {

  if (!yaxis %in% c("frequency", "percent")) {
    stop("yaxis must be either 'frequency' or 'percent'.")
  }

  if (!is.null(group_col) && is.null(cluster_col)) {
    stop("cluster_col must be provided when group_col is specified.")
  }

  # Calculate clonotype abundance
  plt_dat <- calc_abundance(
    input       = input,
    cluster_col = cluster_col,
    data_col    = data_col,
    prefix      = ".",
    return_df   = TRUE
  )

  # Set axis labels
  y_lab <- str_c(data_col, " ", yaxis)

  abun_col <- ".pct"

  if (identical(yaxis, "frequency")) {
    abun_col <- ".freq"
  }

  plt_dat <- tibble::as_tibble(plt_dat, rownames = CELL_COL)
  plt_dat <- dplyr::filter(plt_dat, !is.na(!!sym(data_col)))

  keep_cols <- c(cluster_col, group_col, data_col, abun_col)

  plt_dat <- dplyr::distinct(plt_dat, !!!syms(keep_cols))

  # Rank isotypes
  plt_dat <- dplyr::group_by(plt_dat, !!sym(data_col))

  rnk <- dplyr::summarize(plt_dat, mn = mean(!!sym(abun_col)))
  rnk <- dplyr::arrange(rnk, desc(mn))
  rnk <- pull(rnk, data_col)

  plt_dat <- dplyr::ungroup(plt_dat)

  # Create grouped boxplot
  if (!is.null(group_col)) {
    plt_dat <- .set_lvls(plt_dat, group_col, plot_lvls)
    plt_dat <- .set_lvls(plt_dat, data_col, rnk)

    res <- .create_boxes(
      plt_dat,
      x      = data_col,
      y      = abun_col,
      .color = group_col,
      .fill  = group_col,
      alpha  = 0.5,
      clrs   = plot_colors,
      outlier.color = NA,
      ...
    ) +
      geom_jitter(position = position_jitterdodge(jitter.width = 0.05)) +
      labs(y = y_lab) +
      theme(legend.position = "right")

    return(res)
  }

  # Create bar graph
  # .create_bars reverses level order
  x_col <- data_col

  if (!is.null(cluster_col)) {
    plt_dat <- .set_lvls(plt_dat, cluster_col, plot_lvls)

    x_col <- cluster_col
  }

  plt_dat <- .set_lvls(plt_dat, data_col, rev(rnk))

  new_args <- list(
    df_in = plt_dat,
    x     = x_col,
    y     = abun_col,
    y_ttl = y_lab,
    .fill = data_col,
    clrs  = plot_colors,
    ang   = 45,
    hjst  = 1,
    ...
  )

  # When cluster_col is provided set default position to dodge
  if (!is.null(cluster_col) && is.null(new_args$position)) {
    new_args$position <- position_dodge(preserve = "single")
  }

  res <- purrr::lift_dl(.create_bars)(new_args)

  res
}

