#' Calculate frequency of cell groups present in the object
#'
#' Calculate the frequency of each cell label present in the provided meta.data
#' column. This is useful for comparing the proportion of cells belonging to
#' different samples, cell types, clonotypes, isotypes, etc.
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
#' res <- calc_frequency(
#'   vdj_so,
#'   data_col = "clonotype_id"
#' )
#'
#' head(res@meta.data, 1)
#'
#' # Group cells based on meta.data column before calculating abundance
#' res <- calc_frequency(
#'   vdj_sce,
#'   data_col = "clonotype_id",
#'   cluster_col = "orig.ident"
#' )
#'
#' head(res@colData, 1)
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
#' head(res@meta.data, 1)
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
calc_frequency <- function(input, data_col, cluster_col = NULL, prefix = paste0(data_col, "_"),
                           return_df = FALSE) {

  # Format input data
  meta <- .get_meta(input)
  vdj  <- dplyr::filter(meta, !is.na(!!sym(data_col)))

  vdj <- dplyr::select(
    vdj,
    all_of(c(CELL_COL, data_col, cluster_col))
  )

  # Calculate clonotype abundance
  vdj <- .calc_freq(
    df_in      = vdj,
    cell_col   = CELL_COL,
    dat_col    = data_col,
    clust_col  = cluster_col,
    out_prefix = prefix
  )

  vdj <- dplyr::select(vdj, -all_of(data_col))

  freq_clmn <- paste0(prefix, "freq")
  grp_clmn <- paste0(prefix, "grp")

  vdj <- dplyr::mutate(
    vdj,
    !!sym(grp_clmn) := .calc_freq_grp(!!sym(freq_clmn))
  )

  # Format results
  res <- dplyr::left_join(meta, vdj, by = CELL_COL)

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
.calc_freq <- function(df_in, cell_col, dat_col, clust_col = NULL, out_prefix = "") {

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
  labs <- tibble(
    x   = uniq_x,
    grp = ntile(uniq_x, n_grps)
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
#' @param n_clones Number of top clonotypes to plot (default is 10). If type is set to 'line',
#' this will specify the number of clonotypes to label (default is 3).
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
#' plot_clonal_abundance(
#'   vdj_so,
#'   data_col = "clonotype_id"
#' )
#'
#' # Plot clonotype abundance separately for each cell cluster
#' plot_clonal_abundance(
#'   vdj_sce,
#'   cluster_col = "orig.ident"
#' )
#'
#' # Plot the frequency of each clonotype instead of percentage
#' plot_clonal_abundance(
#'   vdj_sce,
#'   cluster_col = "orig.ident",
#'   yaxis = "frequency"
#' )
#'
#' # Specify colors to use for each cell cluster
#' plot_clonal_abundance(
#'   vdj_so,
#'   cluster_col = "orig.ident",
#'   plot_colors = c(avid_1 = "blue", avid_2 = "red")
#' )
#'
#' # Specify order to use for plotting cell clusters
#' plot_clonal_abundance(
#'   vdj_sce,
#'   cluster_col = "orig.ident",
#'   plot_lvls = c("avid_2", "avid_1")
#' )
#'
#' # Specify the number of top clonotypes to plot
#' plot_clonal_abundance(
#'   vdj_so,
#'   n_clones = 5
#' )
#'
#' #' # Create line graph
#' # use n_clones to set the number of clonotypes to label
#' plot_clonal_abundance(
#'   vdj_so,
#'   cluster_col = "orig.ident",
#'   type = "line",
#'   n_clones = 3
#' )
#'
#' @export
plot_clonal_abundance <- function(input, cluster_col = NULL,
                                  clonotype_col = "clonotype_id", type = "bar",
                                  yaxis = "percent", plot_colors = NULL,
                                  plot_lvls = names(plot_colors), n_clones = NULL,
                                  label_aes = list(), facet_rows = 1,
                                  facet_scales = "free_x", ...) {

  if (!yaxis %in% c("frequency", "percent")) {
    stop("yaxis must be either 'frequency' or 'percent'.")
  }

  if (identical(type, "bar")) {
    n_clones <- n_clones %||% 10

    if (n_clones <= 0) stop("n_clones must be >0.")

  } else if (identical(type, "line")) {
    n_clones <- n_clones %||% 3

    if (n_clones < 0) stop("n_clones must be >=0.")

  } else {
    stop("type must be either 'bar' or 'line'.")
  }

  if (identical(type, "bar") && n_clones <= 0) {
    stop("If type is set to 'bar', n_clones must be >0.")
  }

  # Calculate clonotype abundance
  plt_dat <- calc_frequency(
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
    n         = n_clones,
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
  if (n_clones > 0) {
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

#' Plot frequency of cell groups/labels present in the object
#'
#' Plot the frequency of each cell label present in the provided meta.data
#' column. This is useful for comparing the proportion of cells belonging to
#' different samples, cell types, isotypes, etc. To compare clonotype
#' abundance, use the plot_clonal_abundance() function.
#'
#' @param input Single cell object or data.frame containing V(D)J data. If a
#' data.frame is provided, the cell barcodes should be stored as row names.
#' @param data_col meta.data column containing cell labels to use for
#' calculating frequency, e.g. cell types, isotypes, etc. This function is not
#' designed to plot clonal abundance, use the plot_clonal_abundance() function
#' for this purpose.
#' @param cluster_col meta.data column containing cluster IDs (or patients,
#' treatment conditions, etc.) to use when calculating frequency. Calculations
#' will be performed separately for each cluster.
#' @param group_col meta.data column to use for grouping cluster IDs present in
#' cluster_col. This is useful when there are multiple replicates or patients
#' for each treatment condition.
#' @param yaxis Units to plot on the y-axis, either 'frequency' or 'percent'
#' @param plot_colors Character vector containing colors for plotting
#' @param plot_lvls Character vector containing levels for ordering
#' @param ... Additional arguments to pass to ggplot2, e.g. color, fill, size,
#' linetype, etc.
#' @return ggplot object
#' @export
plot_frequency <- function(input, data_col, cluster_col = NULL,
                           group_col = NULL, yaxis = "percent",
                           plot_colors = NULL, plot_lvls = names(plot_colors),
                           ...) {

  if (!yaxis %in% c("frequency", "percent")) {
    stop("yaxis must be either 'frequency' or 'percent'.")
  }

  .chk_group_cols(cluster_col, group_col)

  # Calculate clonotype abundance
  plt_dat <- calc_frequency(
    input       = input,
    cluster_col = cluster_col,
    data_col    = data_col,
    prefix      = ".",
    return_df   = TRUE
  )

  # Set axis labels
  y_lab <- paste0(data_col, " ", yaxis)

  abun_col <- ".pct"

  if (identical(yaxis, "frequency")) {
    abun_col <- ".freq"
  }

  plt_dat <- tibble::as_tibble(plt_dat, rownames = CELL_COL)
  plt_dat <- dplyr::filter(plt_dat, !is.na(!!sym(data_col)))

  keep_cols <- c(cluster_col, group_col, data_col, abun_col)

  plt_dat <- dplyr::distinct(plt_dat, !!!syms(keep_cols))

  # Rank values in data_col
  plt_dat <- dplyr::group_by(plt_dat, !!sym(data_col))

  rnk <- dplyr::summarize(plt_dat, mn = mean(!!sym(abun_col)))
  rnk <- dplyr::arrange(rnk, desc(.data$mn))
  rnk <- pull(rnk, data_col)

  plt_dat <- dplyr::ungroup(plt_dat)

  # Plot arguments
  gg_args <- list(y = abun_col, clrs = plot_colors, ...)

  # Create grouped boxplot
  if (!is.null(group_col)) {
    plt_dat <- .set_lvls(plt_dat, group_col, plot_lvls)
    plt_dat <- .set_lvls(plt_dat, data_col, rnk)

    gg_args$alpha         <- gg_args$alpha %||% 0.5
    gg_args$outlier.color <- gg_args$outlier.color %||% NA

    more_args <- list(
      df_in  = plt_dat,
      x      = data_col,
      .color = group_col,
      .fill  = group_col
    )

    gg_args <- append(gg_args, more_args)

    res <- purrr::lift_dl(.create_boxes)(gg_args) +
      ggplot2::geom_jitter(
        position = ggplot2::position_jitterdodge(jitter.width = 0.05)
      ) +
      labs(y = y_lab) +
      theme(legend.position = "right")

    return(res)
  }

  # Create bar graph
  # .create_bars reverses level order
  x_col <- data_col

  if (!is.null(cluster_col)) {
    plt_dat <- .set_lvls(plt_dat, cluster_col, plot_lvls)
    x_col   <- cluster_col
  }

  plt_dat <- .set_lvls(plt_dat, data_col, rev(rnk))

  more_args <- list(
    df_in = plt_dat,
    x     = x_col,
    y_ttl = y_lab,
    .fill = data_col,
    ang   = 45,
    hjst  = 1
  )

  gg_args <- append(gg_args, more_args)

  # When cluster_col is provided set default position to dodge
  if (!is.null(cluster_col)) {
    gg_pos <- ggplot2::position_dodge(preserve = "single")

    gg_args$position <- gg_args$position %||% gg_pos
  }

  res <- purrr::lift_dl(.create_bars)(gg_args)

  res
}

