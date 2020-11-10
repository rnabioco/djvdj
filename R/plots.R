#' Create two dimensional scatter plot
#'
#' @param sobj_in Seurat object or data.frame containing data for plotting
#' @param x Variable to plot on x-axis
#' @param y Variable to plot on y-axis
#' @param feature Variable to use for coloring points
#' @param data_slot Slot in the Seurat object to pull data
#' @param pt_size Point size
#' @param plot_colors Vector of colors to use for plot
#' @param feat_lvls Levels to use for ordering feature
#' @param facet_col Variable to use for splitting plot into facets
#' @param facet_lvls Levels to use for ordering plot facets
#' @param min_pct Minimum value to plot for feature
#' @param max_pct Maximum value to plot for feature
#' @param na_color Color to use for missing values
#' @param lm_line Add regression line to plot
#' @param cor_label Position of correlation coefficient label
#' @param label_size Size of correlation coefficient label
#' @param ... Additional parameters to pass to facet_wrap
#' @return ggplot object
plot_features <- function(sobj_in, x = "UMAP_1", y = "UMAP_2", feature, data_slot = "data",
                          pt_size = 0.25, plot_colors = NULL, feat_lvls = NULL, facet_col = NULL,
                          facet_lvls = NULL, min_pct = NULL, max_pct = NULL, na_color = "grey90",
                          lm_line = FALSE, cor_label = c(0.8, 0.9), label_size = 3.7, ...) {

  if (x == y) {
    stop("'x' and 'y' must be different.")
  }

  # Format imput data
  meta_df <- sobj_in

  if ("Seurat" %in% class(sobj_in)) {
    vars <- c(x, y, feature)

    if (!is.null(facet_col)) {
      vars <- c(vars, facet_col)
    }

    meta_df <- Seurat::FetchData(
      sobj_in,
      vars = unique(vars),
      slot = as.character(data_slot)
    )

    meta_df <- tibble::as_tibble(meta_df, rownames = ".cell_id")
  }

  if (!feature %in% colnames(meta_df)) {
    stop(paste(feature, "not found in object."))
  }

  # Rename features
  if (!is.null(names(feature))) {
    meta_df  <- dplyr::rename(meta_df, !!!syms(feature))
    feature <- names(feature)
  }

  if (!is.null(names(x))) {
    meta_df <- dplyr::rename(meta_df, !!!syms(x))
    x <- names(x)
  }

  if (!is.null(names(y))) {
    meta_df <- dplyr::rename(meta_df, !!!syms(y))
    y <- names(y)
  }

  # Adjust values based on min_pct and max_pct
  meta_df <- .set_lims(
    meta_df,
    ft = feature,
    mn = min_pct,
    mx = max_pct
  )

  # Set feature and facet order
  meta_df <- .set_lvls(meta_df, feature, feat_lvls)

  if (!is.null(facet_col) && length(facet_col) == 1) {
    meta_df <- .set_lvls(meta_df, facet_col, facet_lvls)
  }

  # Create scatter plot
  meta_df <- dplyr::arrange(meta_df, !!sym(feature))

  res <- ggplot2::ggplot(
    meta_df,
    ggplot2::aes(!!sym(x), !!sym(y), color = !!sym(feature))
  ) +
    ggplot2::geom_point(size = pt_size) +
    vdj_theme()

  # Add regression line
  if (lm_line) {
    res <- .add_lm(
      res,
      lab_pos = cor_label,
      lab_size = label_size
    )
  }

  # Set feature colors
  if (!is.null(plot_colors)) {
    if (is.numeric(meta_df[[feature]])) {
      res <- res +
        ggplot2::scale_color_gradientn(
          colors   = plot_colors,
          na.value = na_color
        )

    } else {
      res <- res +
        ggplot2::scale_color_manual(
          values   = plot_colors,
          na.value = na_color
        )
    }
  }

  # Split plot into facets
  if (!is.null(facet_col)) {
    if (length(facet_col) == 1) {
      res <- res +
        ggplot2::facet_wrap(~ !!sym(facet_col), ...)

    } else if (length(facet_col) == 2) {
      eq <- paste0(facet_col[1], " ~ ", facet_col[2])

      res <- res +
        ggplot2::facet_grid(stats::as.formula(eq), ...)
    }
  }

  res
}


#' Create bar graphs summarizing cell counts for clusters
#'
#' @param sobj_in Seurat object or data.frame containing data
#' @param x meta.data column containing clusters to use for creating bars
#' @param fill_col meta.data column to use for coloring bars
#' @param facet_col meta.data column containing groups to use for spliting plot
#' into facets
#' @param yaxis Units to plot on the y-axis, either "fraction" or "count"
#' @param plot_colors Character vector containing colors for plotting
#' @param plot_lvls Character vector containing levels for ordering
#' @param na_color Color to use for missing values
#' @param n_label Include label showing the number of cells represented by each
#' bar
#' @param label_aes Named list providing additional label aesthetics
#' (color, size, etc.)
#' @param facet_rows The number of facet rows. Use this argument if facet_col
#' is specified
#' @param facet_scales If facet_col is used, this argument passes a scales
#' specification to facet_wrap. Can be "fixed", "free", "free_x", or "free_y"
#' @param ... Additional arguments to pass to ggplot
#' @return ggplot object
#' @export
plot_cell_count <- function(sobj_in, x, fill_col = NULL, facet_col = NULL, yaxis = "fraction",
                            plot_colors = NULL, plot_lvls = NULL, na_color = "grey90", n_label = TRUE,
                            label_aes = list(), facet_rows = 1, facet_scales = "free_x", ...) {

  # Set y-axis unit
  y_types <- c("fraction", "count")

  if (!yaxis %in% y_types) {
    stop("yaxis must be either fraction or count.")
  }

  y_types <- purrr::set_names(c("fill", "identity"), y_types)

  bar_pos <- y_types[[yaxis]]

  # Format meta.data for plotting
  meta_df <- sobj_in

  if ("Seurat" %in% class(sobj_in)) {
    meta_df <- sobj_in@meta.data %>%
      as_tibble(rownames = ".cell_id")
  }

  meta_df <- .set_lvls(meta_df, x, plot_lvls)

  # Create bar graphs
  res <- meta_df %>%
    ggplot2::ggplot(ggplot2::aes(!!sym(x)))

  if (!is.null(fill_col)) {
    res <- ggplot2::ggplot(
      meta_df,
      ggplot2::aes(!!sym(x), fill = !!sym(fill_col))
    )
  }

  res <- res +
    ggplot2::geom_bar(position = bar_pos, ...) +
    ggplot2::labs(y = yaxis) +
    vdj_theme() +
    ggplot2::theme(
      axis.title.x = ggplot2::element_blank(),
      axis.text.x  = ggplot2::element_text(angle = 45, hjust = 1)
    )

  if (!is.null(facet_col)) {
    res <- res +
      ggplot2::facet_wrap(
        stats::as.formula(paste0("~ ", facet_col)),
        nrow   = facet_rows,
        scales = facet_scales
      )
  }

  # Add n labels
  if (n_label) {
    cell_counts <- dplyr::group_by(meta_df, !!sym(x))

    if (!is.null(facet_col)) {
      cell_counts <- group_by(cell_counts, !!sym(facet_col), .add = T)
    }

    cell_counts <- summarize(
      cell_counts,
      cell_count = paste0("n = ", dplyr::n_distinct(.data$.cell_id)),
      .groups = "drop"
    )

    res <- res +
      ggplot2::geom_text(
        ggplot2::aes(!!sym(x), 0, label = .data$cell_count),
        data          = cell_counts,
        inherit.aes   = F,
        check_overlap = T,
        size          = 3,
        angle         = 90,
        hjust         = -0.1
      )

    res <- .add_aes(res, label_aes, 2)
  }

  # Set plot colors
  if (!is.null(plot_colors)) {
    res <- res +
      ggplot2::scale_fill_manual(values = plot_colors, na.value = na_color)
  }

  res
}


#' Plot read support for chains
#'
#' @param sobj_in Seurat object
#' @param data_cols meta.data columns containing UMI and/or read counts
#' @param chain_col meta.data column containing chains. If chain_col is
#' provided, reads and UMIs will be plotted separately for each chain
#' @param cluster_col meta.data column containing cluster IDs. If cluster_col
#' is provided reads and UMIs will be plotted separately for each cluster.
#' @param type Type of plot to create, can be 'violin', 'histogram' or 'density'
#' @param plot_colors Character vector containing colors for plotting
#' @param plot_lvls Character vector containing levels for ordering
#' @param facet_rows The number of facet rows. Use this argument if both
#' chain_col and cluster_col are provided
#' @param facet_scales This argument passes a scales specification to
#' facet_wrap. Can be "fixed", "free", "free_x", or "free_y"
#' @param sep Separator to use for expanding data_cols
#' @param ... Additional arguments to pass to ggplot
#' @return ggplot object
#' @export
plot_reads <- function(sobj_in, data_cols = c("reads", "umis"), chain_col = NULL, cluster_col = NULL,
                       type = "violin", plot_colors = NULL, plot_lvls = NULL, facet_rows = 1,
                       facet_scales = "free_x", ..., sep = ";") {

  if (!type %in% c("violin", "histogram", "density")) {
    stop("type must be 'violin', 'histogram' or 'density'.")
  }

  # calculate average reads/umis each chain for each cell
  gg_df <- summarize_chains(
    sobj_in,
    data_cols    = data_cols,
    fn           = mean,
    chain_col    = chain_col,
    include_cols = cluster_col,
    sep          = sep
  )

  # Format data for plotting
  gg_df <- tidyr::pivot_longer(
    gg_df,
    cols      = all_of(data_cols),
    names_to  = "key",
    values_to = "counts"
  )

  if (!is.null(cluster_col)) {
    gg_df <- .set_lvls(gg_df, cluster_col, plot_lvls)
  }

  # Calculate median counts for each chain
  box_stats <- dplyr::group_by(
    gg_df,
    !!!syms(c(chain_col, cluster_col, "key"))
  )

  box_stats <- dplyr::summarize(
    box_stats,
    med     = stats::median(.data$counts),
    .groups = "drop"
  )

  # Set chain_col
  if (!is.null(chain_col)) {
    lvls <- dplyr::group_by(gg_df, !!sym(chain_col))

    lvls <- dplyr::summarise(
      lvls,
      n       = dplyr::n_distinct(.data$.cell_id),
      .groups = "drop"
    )

    lvls <- dplyr::arrange(lvls, desc(n))
    lvls <- pull(lvls, chain_col)

    gg_df <- .set_lvls(gg_df, chain_col, lvls)

    chain_col <- sym(chain_col)

  } else if (!is.null(cluster_col)) {
    chain_col <- sym(cluster_col)

  } else {
    chain_col <- "counts"
    box_stats <- dplyr::mutate(box_stats, !!sym(chain_col) := chain_col)
  }

  # Create violin plots
  if (type == "violin") {
    res <- ggplot2::ggplot(gg_df, ggplot2::aes(!!chain_col, .data$counts)) +
      ggplot2::geom_violin(
        ggplot2::aes(
          color = !!chain_col,
          fill  = !!chain_col,
          alpha = .data$key
        ),
        position = "identity",
        ...
      ) +
      ggplot2::geom_point(
        ggplot2::aes(!!sym(chain_col), .data$med),
        data        = box_stats,
        show.legend = FALSE
      ) +
      ggplot2::scale_y_log10(
        labels = scales::trans_format(
          "log10",
          scales::label_math()
        ),
        breaks = 10^(-10:10)
      ) +
      ggplot2::theme(axis.title.x = ggplot2::element_blank())

  # Create histogram
  } else {
    res <- ggplot2::ggplot(
      gg_df,
      ggplot2::aes(
        .data$counts,
        color = !!chain_col,
        fill  = !!chain_col,
        alpha = .data$key
      )
    ) +
      ggplot2::scale_x_log10(
        labels = scales::trans_format(
          "log10",
          scales::label_math()
        ),
        breaks = 10^(-10:10)
      )

    if (type == "histogram") {
      res <- res +
        ggplot2::geom_histogram(position = "identity", color = NA, ...)

    } else {
      res <- res +
        ggplot2::geom_density(...)
    }
  }

  res <- res +
    vdj_theme() +
    ggplot2::theme(legend.title = ggplot2::element_blank())

  if (!is.null(cluster_col) && cluster_col != chain_col) {
    res <- res +
      ggplot2::facet_wrap(
        stats::as.formula(paste0("~ ", cluster_col)),
        scales = facet_scales,
        nrow   = facet_rows
      )
  }

  # Adjust theme
  if (length(data_cols) > 1) {
    alphas <- 1 / length(data_cols)
    alphas <- seq(alphas, 1, alphas)

    res <- res +
      ggplot2::scale_alpha_manual(values = alphas) +
      ggplot2::guides(alpha = ggplot2::guide_legend(
        override.aes = list(fill = "black")
      ))
  }

  if (!is.null(plot_colors)) {
    res <- res +
      ggplot2::scale_color_manual(values = plot_colors) +
      ggplot2::scale_fill_manual(values = plot_colors)
  }

  res
}


#' Plot clonotype abundance
#'
#' @param sobj_in Seurat object containing V(D)J data
#' @param clonotype_col meta.data column containing clonotype IDs to use for
#' calculating clonotype abundance
#' @param cluster_col meta.data column containing cluster IDs to use for
#' grouping cells when calculating clonotype abundance
#' @param type Type of plot to create, can be 'bar' or 'line'
#' @param yaxis Units to plot on the y-axis, either "frequency" or "percent"
#' @param plot_colors Character vector containing colors for plotting
#' @param plot_lvls Character vector containing levels for ordering
#' @param label_col meta.data column containing labels to use for plot
#' @param n_clonotypes Number of clonotypes to label
#' @param color_col meta.data column to use for coloring bars
#' @param label_aes Named list providing additional label aesthetics (color,
#' size, etc.)
#' @param facet_rows The number of facet rows. Use this argument if type =
#' 'bar'
#' @param facet_scales If type = 'bar', this argument passes a scales
#' specification to facet_wrap, can be "fixed", "free", "free_x", or "free_y"
#' @param ... Additional arguments to pass to geom_line
#' @return ggplot object
#' @export
plot_abundance <- function(sobj_in, clonotype_col, cluster_col = NULL, type = "bar",
                           yaxis = "percent", plot_colors = NULL, plot_lvls = NULL, label_col = NULL,
                           n_clonotypes = 10, color_col = NULL, label_aes = list(), facet_rows = 1,
                           facet_scales = "free_x", ...) {

  if (!yaxis %in% c("frequency", "percent")) {
    stop("yaxis must be either 'frequency' or 'percent'.")
  }

  if (!type %in% c("bar", "line")) {
    stop("type must be either 'bar' or 'line'.")
  }

  if (type == "bar" && is.null(label_col)) {
    stop("Must include label_col when type = 'bar'.")
  }

  if (is.null(color_col)) {
    color_col <- cluster_col
  }

  # Calculate clonotype abundance
  # Return seurat since label_col is needed from meta.data
  sobj_in <- djvdj::calc_abundance(
    sobj_in,
    clonotype_col = clonotype_col,
    cluster_col   = cluster_col,
    prefix        = ".",
    return_seurat = TRUE
  )

  data_col <- ".clone_pct"

  if (yaxis == "frequency") {
    data_col <- ".clone_freq"
  }

  meta_df <- sobj_in@meta.data
  meta_df <- tibble::as_tibble(meta_df, rownames = ".cell_id")
  meta_df <- dplyr::filter(meta_df, !is.na(!!sym(clonotype_col)))

  abund_cols <- c(
    cluster_col, clonotype_col, label_col,
    data_col, color_col
  )

  meta_df <- dplyr::select(
    meta_df,
    all_of(abund_cols)
  )

  meta_df <- dplyr::distinct(meta_df)

  # Rank by abundance
  if (!is.null(cluster_col)) {
    meta_df <- .set_lvls(meta_df, cluster_col, plot_lvls)
    meta_df <- dplyr::group_by(meta_df, !!sym(cluster_col))
  }

  meta_df <- dplyr::mutate(
    meta_df,
    rank = dplyr::row_number(dplyr::desc(!!sym(data_col)))
  )

  # Identify top clonotypes
  top_genes <- dplyr::slice_min(
    meta_df,
    order_by  = rank,
    n         = n_clonotypes,
    with_ties = FALSE
  )

  meta_df <- dplyr::ungroup(meta_df)

  # Create bar graph
  if (type == "bar") {
    top_genes <- dplyr::arrange(top_genes, !!sym(data_col))

    top_genes <- .set_lvls(
      df_in = top_genes,
      clmn  = label_col,
      lvls  = dplyr::pull(top_genes, label_col)
    )

    res <- .create_bars(
      df_in = top_genes,
      x     = label_col,
      y     = data_col,
      y_ttl = yaxis,
      .fill = color_col,
      clrs  = plot_colors,
      ang   = 45,
      hjst  = 1
    )

    if (!is.null(cluster_col)){
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
  res <- ggplot2::ggplot(meta_df, ggplot2::aes(rank, !!sym(data_col))) +
    ggplot2::labs(y = yaxis) +
    vdj_theme()

  if (is.null(color_col)) {
    res <- res +
      ggplot2::geom_line(...)

  } else {
    res <- res +
      ggplot2::geom_line(ggplot2::aes(color = !!sym(color_col)), ...)
  }

  if (!is.null(plot_colors)) {
    res <- res +
      ggplot2::scale_color_manual(values = plot_colors)
  }

  # Add labels
  if (!is.null(label_col)) {
    res <- res +
      ggrepel::geom_text_repel(
        ggplot2::aes(label = !!sym(label_col)),
        data          = top_genes,
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


#' Plot repertoire diversity
#'
#' @param sobj_in Seurat object containing V(D)J data or data.frame containing
#' diversity values to plot. If a data.frame is used
#' @param clonotype_col meta.data column containing clonotype IDs to use for
#' calculating clonotype abundance
#' @param cluster_col meta.data column containing cluster IDs to use for
#' grouping cells when calculating clonotype abundance
#' @param method Method to use for calculating diversity. Can also pass a named
#' list to use multiple methods.
#' @param plot_colors Character vector containing colors for plotting
#' @param plot_lvls Character vector containing levels for ordering
#' @param facet_rows The number of facet rows. Use this argument if a list of
#' functions is passed to method
#' @param ... Additional arguments to pass to ggplot
#' @return ggplot object
#' @export
plot_diversity <- function(sobj_in, clonotype_col, cluster_col = NULL,
                           method = abdiv::shannon, plot_colors = NULL, plot_lvls = NULL,
                           facet_rows = 1, ...) {

  if (length(method) > 1 && is.null(names(method))) {
    stop("Must include names if using a list of methods.")
  }

  if (length(method) == 1 && is.null(names(method))) {
    nm <- as.character(substitute(method))
    nm <- dplyr::last(nm)

    method        <- list(method)
    names(method) <- nm
  }

  # Calculate diversity
  if ("Seurat" %in% class(sobj_in)) {
    sobj_in <- calc_diversity(
      sobj_in       = sobj_in,
      clonotype_col = clonotype_col,
      cluster_col   = cluster_col,
      method        = method,
      prefix        = "",
      return_seurat = FALSE
    )
  }

  # Format data for plotting
  include_strips <- TRUE

  if (is.null(cluster_col)) {
    meta_df <- tibble::tibble(
      diversity = unname(sobj_in),
      name      = names(sobj_in),
    )

    include_strips <- FALSE
    cluster_col <- "name"

  } else {
    meta_df <- tidyr::pivot_longer(
      sobj_in,
      cols = all_of(names(method)),
      values_to = "div"
    )

    meta_df <- .set_lvls(meta_df, cluster_col, plot_lvls)
  }

  # Set plot levels
  meta_df <- .set_lvls(meta_df, "name", unique(meta_df$name))

  # Create bar graphs
  res <- ggplot2::ggplot(meta_df, ggplot2::aes(
    !!sym(cluster_col),
    .data$div,
    fill = !!sym(cluster_col)
  )) +
    ggplot2::geom_col(...)

  # Create facets
  if (length(method) > 1) {
    res <- res +
      ggplot2::facet_wrap(~ name, nrow = facet_rows, scales = "free")
  }

  # Set plot colors
  if (!is.null(plot_colors)) {
    res <- res +
      ggplot2::scale_fill_manual(values = plot_colors)
  }

  # Set theme
  res <- res +
    vdj_theme() +
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


#' Plot repertoire overlap
#'
#' @param sobj_in Seurat object containing V(D)J data or matrix containing
#' similarity values to plot
#' @param clonotype_col meta.data column containing clonotype IDs to use for
#' calculating overlap
#' @param cluster_col meta.data column containing cluster IDs to use for
#' calculating overlap
#' @param method Method to use for calculating similarity between clusters
#' @param plot_colors Character vector containing colors for plotting
#' @param ... Additional arguments to pass to ggplot
#' @return ggplot object
#' @export
plot_similarity <- function(sobj_in, clonotype_col, cluster_col, method = abdiv::jaccard,
                            plot_colors = NULL, ...) {

  # Calculate similarity
  if ("Seurat" %in% class(sobj_in)) {
    if (is.null(clonotype_col) || is.null(cluster_col)) {
      stop("If a Seurat object is provided, clonotype_col and cluster_col must be specified.")
    }

    sobj_in <- calc_similarity(
      sobj_in       = sobj_in,
      clonotype_col = clonotype_col,
      cluster_col   = cluster_col,
      method        = method,
      prefix        = "",
      return_seurat = FALSE
    )
  }

  sim_col <- as.character(substitute(method))
  sim_col <- dplyr::last(sim_col)

  var_lvls <- unique(c(rownames(sobj_in), colnames(sobj_in)))
  var_lvls <- sort(var_lvls)

  gg_df <- tibble::as_tibble(sobj_in, rownames = "Var1")

  gg_df <- tidyr::pivot_longer(
    gg_df,
    cols      = -.data$Var1,
    names_to  = "Var2",
    values_to = sim_col
  )

  # Set Var levels
  gg_df <- dplyr::mutate(
    gg_df,
    Var1 = factor(.data$Var1, levels = rev(var_lvls)),
    Var2 = factor(.data$Var2, levels = var_lvls)
  )

  # Create heatmap
  res <- .create_heatmap(
    gg_df,
    x     = "Var1",
    y     = "Var2",
    .fill = sim_col,
    clrs  = plot_colors,
    ...
  )

  res
}


#' Plot V(D)J gene usage
#'
#' @param sobj_in Seurat object containing V(D)J data
#' @param gene_cols meta.data column containing genes used for each clonotype
#' @param cluster_col meta.data column containing cell clusters to use for
#' calculating gene usage
#' @param chain Chain to use for calculating gene usage, set to NULL to include
#' all chains
#' @param type Type of plot to create, can be either 'heatmap' or 'bar'
#' @param plot_colors Character vector containing colors for plotting
#' @param plot_genes Character vector containing genes to plot
#' @param n_genes Number of top genes to plot based on average usage
#' @param plot_lvls Levels to use for ordering clusters
#' @param yaxis Units to plot on the y-axis, either "frequency" or "percent"
#' @param chain_col meta.data column containing chains for each cell
#' @param sep Separator to use for expanding gene_cols
#' @param ... Additional arguments to pass to ggplot
#' @return ggplot object
#' @export
plot_usage <- function(sobj_in, gene_cols, cluster_col = NULL, chain = NULL, type = "bar",
                       plot_colors = NULL, plot_genes = NULL, n_genes = NULL, plot_lvls = NULL,
                       yaxis = "percent", chain_col = NULL, sep = ";", ...) {

  # Check inputs
  if (!type %in% c("heatmap", "bar")) {
    stop("type must be either 'heatmap' or 'bar'.")
  }

  if (length(gene_cols) > 2) {
    stop("Cannot specify more than two values for gene_cols.")
  }

  if (!yaxis %in% c("percent", "frequency")) {
    stop("yaxis must be either 'percent' or 'frequency'.")
  }

  # Set y-axis
  usage_col <- "pct"

  if (yaxis == "frequency") {
    usage_col <- "freq"
  }

  # Calculate gene usage
  gg_df <- calc_usage(
    sobj_in,
    gene_cols   = gene_cols,
    cluster_col = cluster_col,
    chain       = chain,
    chain_col   = chain_col,
    sep         = sep
  )

  gg_df <- dplyr::filter(gg_df, dplyr::across(
    dplyr::all_of(gene_cols),
    ~ .x != "None"
  ))

  # Order genes by average usage
  top_genes <- gg_df

  if (!is.null(cluster_col)) {
    top_genes <- dplyr::group_by(top_genes, !!sym(cluster_col))
  }

  if (!is.null(n_genes)) {
    top_genes <- dplyr::slice_max(
      top_genes,
      order_by  = !!sym(usage_col),
      n         = n_genes,
      with_ties = FALSE
    )
  }

  top_genes <- dplyr::group_by(top_genes, !!!syms(gene_cols))

  top_genes <- dplyr::summarize(
    top_genes,
    usage = mean(!!sym(usage_col)), .groups = "drop"
  )

  top_genes <- dplyr::arrange(top_genes, .data$usage)

  if (!is.null(n_genes)) {
    gene_lst <- unlist(top_genes[, gene_cols], use.names = FALSE)

    gg_df <- dplyr::filter(gg_df, dplyr::across(
      dplyr::all_of(gene_cols),
      ~ .x %in% gene_lst
    ))
  }

  if (length(gene_cols) == 1) {
    lvls  <- pull(top_genes, gene_cols)
    gg_df <- .set_lvls(gg_df, gene_cols, lvls)
  }

  # Order clusters
  gg_df <- .set_lvls(gg_df, cluster_col, plot_lvls)

  # Filter genes to plot
  if (!is.null(plot_genes)) {
    gg_df <- dplyr::filter(gg_df, dplyr::across(
      dplyr::all_of(gene_cols),
      ~ .x %in% plot_genes
    ))
  }

  # Create heatmap or bar graph for single gene
  if (length(gene_cols) == 1) {
    if (type == "bar") {
      gg_df <- dplyr::filter(gg_df, !!sym(usage_col) > 0)

      res <- .create_bars(
        gg_df,
        x     = gene_cols,
        y     = usage_col,
        .fill = cluster_col,
        clrs  = plot_colors,
        y_ttl = yaxis,
        ...
      )

      return(res)
    }

    res <- .create_heatmap(
      gg_df,
      x     = cluster_col,
      y     = gene_cols,
      .fill = usage_col,
      clrs  = plot_colors,
      ttl   = yaxis,
      ...
    )

    return(res)
  }

  # Create heatmap for two genes
  grps <- NULL

  if (!is.null(cluster_col)) {
    gg_df <- dplyr::group_by(gg_df, !!sym(cluster_col))

    grps <- pull(gg_df, cluster_col)
    grps <- sort(unique(grps))
  }

  res <- dplyr::group_split(gg_df)
  names(res) <- grps

  # Format data for plotting
  res <- purrr::map(res, ~ {
    res <- dplyr::select(.x, dplyr::all_of(c(gene_cols, usage_col)))

    res <- tidyr::pivot_wider(
      res,
      names_from  = gene_cols[2],
      values_from = all_of(usage_col)
    )

    res <- dplyr::mutate(res, dplyr::across(
      -dplyr::all_of(gene_cols[1]),
      ~ tidyr::replace_na(.x, 0)
    ))

    res <- tidyr::pivot_longer(
      res,
      cols      = -dplyr::all_of(gene_cols[1]),
      names_to  = gene_cols[2],
      values_to = usage_col
    )

    res <- dplyr::arrange(res, dplyr::desc(!!sym(usage_col)))

    v1_lvls <- rev(unique(pull(res, gene_cols[1])))
    v2_lvls <- rev(unique(pull(res, gene_cols[2])))

    res <- .set_lvls(res, gene_cols[1], v1_lvls)
    res <- .set_lvls(res, gene_cols[2], v2_lvls)

    res
  })

  # Create heatmaps
  res <- purrr::map(
    res,
    .create_heatmap,
    x     = gene_cols[1],
    y     = gene_cols[2],
    .fill = usage_col,
    clrs  = plot_colors,
    ttl   = yaxis,
    ...
  )

  if (length(res) == 1) {
    return(res[[1]])
  }

  res
}


#' Theme for djvdj plotting functions
#'
#' @param ttl_size Size of axis titles
#' @param txt_size Size of axis text
#' @param ln_size Size of axis lines
#' @param txt_col Color of axis text
#' @param ln_col Color of axis lines
#' @return ggplot theme
#' @export
vdj_theme <- function(ttl_size = 12, txt_size = 8, ln_size = 0.5, txt_col = "black",
                      ln_col = "grey85") {
  res <- ggplot2::theme(
    strip.background  = ggplot2::element_blank(),
    strip.text        = ggplot2::element_text(size = ttl_size),
    panel.background  = ggplot2::element_blank(),
    legend.background = ggplot2::element_blank(),
    legend.title      = ggplot2::element_text(size = ttl_size),
    legend.key        = ggplot2::element_blank(),
    legend.text       = ggplot2::element_text(size = txt_size, color = txt_col),
    axis.line         = ggplot2::element_line(size = ln_size, color = ln_col),
    axis.ticks        = ggplot2::element_line(size = ln_size, color = ln_col),
    axis.text         = ggplot2::element_text(size = txt_size, color = txt_col),
    axis.title        = ggplot2::element_text(size = ttl_size, color = txt_col)
  )

  res
}


#' Set min and max values for column
#'
#' @param df_in Input data.frame
#' @param ft Name of column containing feature values
#' @param mn Minimum percent rank
#' @param mx Maximum percent rank
#' @return data.frame with modified feature values
.set_lims <- function(df_in, ft, mn = NULL, mx = NULL) {

  if (is.null(mn) && is.null(mx)) {
    return(df_in)
  }

  ft <- sym(ft)

  res <- dplyr::mutate(
    df_in,
    pct = dplyr::percent_rank(!!ft)
  )

  if (!is.null(mn)) {
    res <- dplyr::mutate(
      res,
      !!ft := ifelse(.data$pct > mn, !!ft, NA),
      !!ft := ifelse(.data$pct <= mn, min(!!ft, na.rm = TRUE), !!ft)
    )
  }

  if (!is.null(mx)) {
    res <- dplyr::mutate(
      res,
      !!ft := ifelse(.data$pct < mx, !!ft, NA),
      !!ft := ifelse(.data$pct >= mx, max(!!ft, na.rm = TRUE), !!ft)
    )
  }

  res <- dplyr::select(res, -.data$pct)

  res
}


#' Add regression line to ggplot object
#'
#' @param gg_in ggplot object
#' @param lab_pos Position of correlation coefficient label. Set to NULL to
#' omit label
#' @param lab_size Size of label
#' @param ... Additional arguments to pass to geom_smooth
#' @return ggplot object with added regression line
.add_lm <- function(gg_in, lab_pos = NULL, lab_size = 3.5, ...) {

  # Add regression line
  res <- gg_in +
    ggplot2::geom_smooth(
      method   = "lm",
      formula  = y ~ x,
      se       = F,
      color    = "black",
      size     = 0.5,
      linetype = 2,
      ...
    )

  # Calculate correlation
  if (!is.null(lab_pos)) {
    gg_df <- gg_in$data

    x <- as_name(gg_in$mapping$x)
    y <- as_name(gg_in$mapping$y)

    gg_df <- dplyr::mutate(
      gg_df,
      r       = broom::tidy(stats::cor.test(!!sym(x), !!sym(y)))$estimate,
      r       = round(.data$r, digits = 2),
      pval    = broom::tidy(stats::cor.test(!!sym(x), !!sym(y)))$p.value,
      cor_lab = paste0("r = ", .data$r, ", p = ", format(.data$pval, digits = 2)),
      min_x   = min(!!sym(x)),
      max_x   = max(!!sym(x)),
      min_y   = min(!!sym(y)),
      max_y   = max(!!sym(y)),
      lab_x   = (.data$max_x - .data$min_x) * lab_pos[1] + .data$min_x,
      lab_y   = (.data$max_y - .data$min_y) * lab_pos[2] + .data$min_y
    )
  }

  # Add correlation coefficient label
  res <- res +
    ggplot2::geom_text(
      data          = gg_df,
      mapping       = ggplot2::aes(.data$lab_x, .data$lab_y, label = .data$cor_lab),
      color         = "black",
      size          = lab_size,
      check_overlap = T,
      show.legend   = F
    )

  res
}


#' Add list of aes params to ggplot object
#'
#' @param gg_in ggplot object
#' @param add_aes List of aes params to add or override
#' @param lay Layer number to modify
#' @return ggplot object
.add_aes <- function(gg_in, add_aes, lay) {

  # Need to use colour instead of color
  aes_names      <- names(add_aes)
  names(add_aes) <- replace(aes_names, aes_names == "color", "colour")

  # Add aes params
  curr_aes <- gg_in$layers[[lay]]$aes_params

  curr_aes[names(add_aes)]       <- add_aes
  gg_in$layers[[lay]]$aes_params <- curr_aes

  gg_in
}


#' Create ggplot heatmap
#'
#' @param df_in data.frame
#' @param x Variable to plot on the x-axis
#' @param y Variable to plot on the y-axis
#' @param .fill Variable to use for the fill color
#' @param clrs Vector of colors for plotting
#' @param na_color Color to use for missing values
#' @param ttl Legend title
#' @param ang Angle of x-axis text
#' @param hjst Horizontal justification for x-axis text
#' @param ... Additional arguments to pass to ggplot
#' @return ggplot object
.create_heatmap <- function(df_in, x, y, .fill, clrs = NULL, na_color = "white",
                            ttl = .fill, ang = 45, hjst = 1, ...) {

  if (is.null(clrs)) {
    clrs <- "#56B4E9"
  }

  if (length(clrs) == 1) {
    clrs <- c("grey90", clrs)
  }

  if (!is.null(x)) {
    res <- ggplot2::ggplot(
      df_in,
      ggplot2::aes(!!sym(x), !!sym(y), fill = !!sym(.fill))
    )

  } else {
    res <- ggplot2::ggplot(
      df_in,
      ggplot2::aes("sample", !!sym(y), fill = !!sym(.fill))
    )
  }

  res <- res +
    ggplot2::geom_tile(...) +
    ggplot2::guides(fill = ggplot2::guide_colorbar(title = ttl)) +
    ggplot2::scale_fill_gradientn(colors = clrs, na.value = na_color) +
    vdj_theme() +
    ggplot2::theme(
      axis.title  = ggplot2::element_blank(),
      axis.line   = ggplot2::element_blank(),
      axis.ticks  = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = ang, hjust = hjst)
    )

  res
}


#' Create ggplot bar graph
#'
#' @param df_in data.frame
#' @param x Variable to plot on x-axis
#' @param y Variable to plot on y-axis
#' @param .fill Variable to use for fill color
#' @param clrs Vector of colors for plotting
#' @param y_ttl Title for y-axis
#' @param ang Angle of x-axis text
#' @param hjst Horizontal justification for x-axis text
#' @param ... Additional arguments to pass to ggplot
#' @return ggplot object
.create_bars <- function(df_in, x, y, .fill, clrs = NULL, y_ttl = y, ang = 45,
                         hjst = 1, ...) {

  # Reverse bar order
  lvls  <- rev(levels(pull(df_in, x)))
  df_in <- .set_lvls(df_in, x, lvls)

  if (!is.null(.fill)) {
    res <- ggplot2::ggplot(
      df_in,
      ggplot2::aes(!!sym(x), !!sym(y), fill = !!sym(.fill))
    )

  } else {
    res <- ggplot2::ggplot(
      df_in,
      ggplot2::aes(!!sym(x), !!sym(y))
    )
  }

  res <- res +
    ggplot2::geom_col(..., position = "dodge") +
    ggplot2::labs(y = y_ttl) +
    vdj_theme() +
    ggplot2::theme(
      axis.title.x = ggplot2::element_blank(),
      axis.text.x  = ggplot2::element_text(angle = ang, hjust = hjst)
    )

  if (!is.null(clrs)) {
    res <- res +
      ggplot2::scale_fill_manual(values = clrs)
  }

  res
}


#' Set column levels
#'
#' @param df_in data.frame
#' @param clmn Column to modify
#' @param lvls Levels
#' @return data.frame
.set_lvls <- function(df_in, clmn, lvls) {

  if (!is.null(lvls) && !is.null(clmn)) {
    dat <- pull(df_in, clmn)

    if (!is.character(dat) && !is.factor(dat)) {
      warning("Plot levels were not modified, levels are only modified for characters and factors.")

      return(df_in)
    }

    if (!all(pull(df_in, clmn) %in% lvls)) {
      stop(paste0("Not all labels in ", clmn, " are included in plot_levels."))
    }

    df_in <- dplyr::mutate(
      df_in,
      !!sym(clmn) := factor(!!sym(clmn), levels = lvls)
    )
  }

  df_in
}

