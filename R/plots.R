#' Create 2D feature plot
#'
#' @export
plot_features <- function(input, ...) {

  UseMethod("plot_features", input)

}

#' @rdname plot_features
#' @param input Single cell object or data.frame containing V(D)J data. If a
#' data.frame is provided, the cell barcodes should be stored as row names.
#' @param x Variable to plot on x-axis
#' @param y Variable to plot on y-axis
#' @param feature Variable to use for coloring points
#' @param data_slot Slot in the Seurat object to pull data
#' @param plot_colors Vector of colors to use for plot
#' @param plot_lvls Levels to use for ordering feature
#' @param min_q Minimum quantile cutoff for color scale.
#' @param max_q Maximum quantile cutoff for color scale.
#' @param na_color Color to use for missing values
#' @param ... Additional arguments to pass to geom_point()
#' @return ggplot object
#' @export
plot_features.default <- function(input, x = "UMAP_1", y = "UMAP_2", feature, plot_colors = NULL,
                                  plot_lvls = NULL, min_q = NULL, max_q = NULL, na_color = "grey90", ...) {

  # Check arguments
  if (x == y) {
    stop("'x' and 'y' must be different")
  }

  plt_dat <- .get_meta(input)

  plt_vars  <- c(x, y, feature)
  miss_vars <- plt_vars[!plt_vars %in% colnames(plt_dat)]

  if (length(miss_vars) > 0) {
    stop(paste(miss_vars, collapse = ", "), " not found in object")
  }

  # Adjust values based on min_q and max_q
  if ((!is.null(min_q) || !is.null(max_q)) && is.numeric(plt_dat[[feature]])) {
    plt_dat <- .set_lims(
      plt_dat,
      ft = feature,
      mn = min_q,
      mx = max_q
    )
  }

  # Set feature and facet order
  plt_dat <- .set_lvls(plt_dat, feature, plot_lvls)

  # Create scatter plot
  # To add outline for each cluster create separate layers
  res <- arrange(plt_dat, !!sym(feature))

  res <- ggplot2::ggplot(
    res,
    ggplot2::aes(!!sym(x), !!sym(y), color = !!sym(feature))
  ) +
    ggplot2::geom_point(...)

  # Set feature colors
  if (!is.null(plot_colors)) {
    if (is.numeric(plt_dat[[feature]])) {
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

  # Set theme
  res <- res +
    djvdj_theme()

  if (!is.numeric(plt_dat[[feature]])) {
    res <- res +
      ggplot2::guides(
        color = ggplot2::guide_legend(override.aes = list(size = 3))
      )
  }

  res
}

#' @rdname plot_features
#' @importFrom Seurat FetchData
#' @export
plot_features.Seurat <- function(input, x = "UMAP_1", y = "UMAP_2", feature, data_slot = "data",
                                 plot_colors = NULL, plot_lvls = NULL, min_q = NULL, max_q = NULL,
                                 na_color = "grey90", ...) {

  # Fetch variables and add to meta.data
  # want input data to include meta.data and any features from FetchData
  plt_vars <- c(x, y, feature)

  suppressWarnings({
    plt_dat <- Seurat::FetchData(
      input,
      vars = unique(plt_vars),
      slot = as.character(data_slot)
    )
  })

  # Format input data
  # keep rownames since default method will create rowname column
  plt_dat <- .merge_meta(input, plt_dat)
  plt_dat <- .get_meta(plt_dat)

  plot_features(
    input       = plt_dat,
    x           = x,
    y           = y,
    feature     = feature,
    plot_colors = plot_colors,
    plot_lvls   = plot_lvls,
    min_q       = min_q,
    max_q       = max_q,
    na_color    = na_color,
    ...
  )
}


#' Plot read support for chains
#'
#' @param input Single cell object or data.frame containing V(D)J data. If a
#' data.frame is provided, the cell barcodes should be stored as row names.
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
#' facet_wrap. Can be "fixed", "free", "free_x", or "free_y".
#' @param sep Separator used for storing per cell V(D)J data
#' @param ... Additional arguments to pass to ggplot
#' @return ggplot object
#' @export
plot_reads <- function(input, data_cols = c("reads", "umis"), chain_col = "chains", cluster_col = NULL,
                       type = "violin", plot_colors = NULL, plot_lvls = NULL, facet_rows = 1,
                       facet_scales = "free_x", sep = ";", ...) {

  if (!type %in% c("violin", "histogram", "density")) {
    stop("type must be 'violin', 'histogram' or 'density'.")
  }

  # Calculate mean reads/umis for each chain for each cell
  plt_dat <- summarize_chains(
    input,
    data_cols    = data_cols,
    fn           = mean,
    chain_col    = chain_col,
    include_cols = cluster_col,
    sep          = sep
  )

  # Format data for plotting
  plt_dat <- tidyr::pivot_longer(
    plt_dat,
    cols      = all_of(data_cols),
    names_to  = "key",
    values_to = "counts"
  )

  if (!is.null(cluster_col)) {
    plt_dat <- .set_lvls(plt_dat, cluster_col, plot_lvls)
  }

  # Calculate median counts for each chain
  box_stats <- dplyr::group_by(
    plt_dat,
    !!!syms(c(chain_col, cluster_col, "key"))
  )

  box_stats <- dplyr::summarize(
    box_stats,
    med     = stats::median(.data$counts),
    .groups = "drop"
  )

  # Set chain_col
  if (!is.null(chain_col)) {
    lvls <- dplyr::group_by(plt_dat, !!sym(chain_col))

    lvls <- dplyr::summarise(
      lvls,
      n       = dplyr::n_distinct(.data$.cell_id),
      .groups = "drop"
    )

    lvls <- dplyr::arrange(lvls, desc(n))
    lvls <- pull(lvls, chain_col)

    plt_dat <- .set_lvls(plt_dat, chain_col, lvls)

    chain_col <- sym(chain_col)

  } else if (!is.null(cluster_col)) {
    chain_col <- sym(cluster_col)

  } else {
    chain_col <- "counts"
    box_stats <- dplyr::mutate(box_stats, !!sym(chain_col) := chain_col)
  }

  # Create violin plots
  if (type == "violin") {
    res <- ggplot2::ggplot(plt_dat, ggplot2::aes(!!chain_col, .data$counts)) +
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
      plt_dat,
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
    djvdj_theme() +
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
#' @param input Single cell object or data.frame containing V(D)J data. If a
#' data.frame is provided, the cell barcodes should be stored as row names.
#' @param cluster_col meta.data column containing cluster IDs to use for
#' grouping cells when calculating clonotype abundance
#' @param clonotype_col meta.data column containing clonotype IDs to use for
#' calculating clonotype abundance
#' @param type Type of plot to create, can be 'bar' or 'line'
#' @param yaxis Units to plot on the y-axis, either 'frequency' or 'percent'
#' @param plot_colors Character vector containing colors for plotting
#' @param plot_lvls Character vector containing levels for ordering
#' @param label_col meta.data column to use for labeling clonotypes. This is
#' useful if clonotype_col contains names that are too long to include on the
#' plot.
#' @param n_clonotypes Number of clonotypes to plot. If type is set to 'line',
#' this will specify the number of clonotypes to label.
#' @param color_col meta.data column to use for coloring bars
#' @param label_aes Named list providing additional label aesthetics (color,
#' size, etc.)
#' @param facet_rows The number of facet rows. Use this argument if type is set
#' to 'bar'
#' @param facet_scales If type is 'bar', this argument passes a scales
#' specification to facet_wrap, can be 'fixed', 'free', 'free_x', or 'free_y'
#' @param ... Additional arguments to pass to ggplot
#' @return ggplot object
#' @export
plot_abundance <- function(input, cluster_col = NULL, clonotype_col = "clonotype_id", type = "bar",
                           yaxis = "percent", plot_colors = NULL, plot_lvls = NULL, label_col = clonotype_col,
                           n_clonotypes = 10, color_col = cluster_col, label_aes = list(), facet_rows = 1,
                           facet_scales = "free_x", ...) {

  if (!yaxis %in% c("frequency", "percent")) {
    stop("yaxis must be either 'frequency' or 'percent'.")
  }

  if (!type %in% c("bar", "line")) {
    stop("type must be either 'bar' or 'line'.")
  }

  if (type == "bar" && n_clonotypes <= 0) {
    stop("If type is set to 'bar', n_clonotypes must be >0.")
  }

  # Calculate clonotype abundance
  plt_dat <- calc_abundance(
    input         = input,
    cluster_col   = cluster_col,
    clonotype_col = clonotype_col,
    prefix        = ".",
    return_df     = TRUE
  )

  dat_col <- ".clone_pct"

  if (yaxis == "frequency") {
    dat_col <- ".clone_freq"
  }

  plt_dat <- tibble::as_tibble(plt_dat, rownames = ".cell_id")
  plt_dat <- dplyr::filter(plt_dat, !is.na(!!sym(clonotype_col)))

  abund_cols <- c(
    cluster_col, clonotype_col,
    label_col,   dat_col,
    color_col
  )

  plt_dat <- dplyr::distinct(plt_dat, !!!syms(abund_cols))

  # Collapse values in label_col so there is a single label for each value in
  # clonotype_col
  if (label_col != clonotype_col) {
    plt_dat <- dplyr::group_by(plt_dat, !!!syms(c(clonotype_col, cluster_col)))

    plt_dat <- dplyr::mutate(
      plt_dat,
      !!sym(label_col) := paste0(!!sym(label_col), collapse = ";")
    )

    plt_dat <- dplyr::distinct(plt_dat)
    plt_dat <- dplyr::ungroup(plt_dat)
  }

  # Add and format label column for plotting
  len <- 25

  plt_dat <- dplyr::rowwise(plt_dat)

  plt_dat <- dplyr::mutate(
    plt_dat,
    .x   = !!sym(clonotype_col),
    .len = length(unlist(strsplit(!!sym(label_col), ""))),
    .lab = strtrim(!!sym(label_col), len),
    .lab = paste0(.data$.lab, ifelse(.data$.len > len, "...", ""))
  )

  plt_dat <- dplyr::ungroup(plt_dat)

  # Rank by abundance
  if (!is.null(cluster_col)) {
    plt_dat <- .set_lvls(plt_dat, cluster_col, plot_lvls)
    plt_dat <- dplyr::group_by(plt_dat, !!sym(cluster_col))

    plt_dat <- dplyr::mutate(
      plt_dat,
      .x = paste0(!!sym(cluster_col), "_", .data$.x)
    )
  }

  plt_dat <- dplyr::mutate(
    plt_dat,
    rank = dplyr::row_number(dplyr::desc(!!sym(dat_col)))
  )

  # Identify top clonotypes
  top_clones <- dplyr::slice_min(
    plt_dat,
    order_by  = rank,
    n         = n_clonotypes,
    with_ties = FALSE
  )

  plt_dat    <- dplyr::ungroup(plt_dat)
  top_clones <- dplyr::ungroup(top_clones)

  # Create bar graph
  if (type == "bar") {
    plt_labs <- set_names(top_clones$.lab, top_clones$.x)

    top_clones <- dplyr::arrange(top_clones, desc(!!sym(dat_col)))

    lvls <- rev(unique(top_clones$.x))

    top_clones <- .set_lvls(
      df_in = top_clones,
      clmn  = ".x",
      lvls  = lvls
    )

    res <- .create_bars(
      df_in = top_clones,
      x     = ".x",
      y     = dat_col,
      y_ttl = yaxis,
      .fill = color_col,
      clrs  = plot_colors,
      ang   = 45,
      hjst  = 1,
      ...
    )

    res <- res +
      ggplot2::scale_x_discrete(labels = plt_labs)

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
  res <- ggplot2::ggplot(plt_dat, ggplot2::aes(rank, !!sym(dat_col))) +
    ggplot2::labs(y = yaxis) +
    djvdj_theme()

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
  if (n_clonotypes > 0) {
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


#' Plot repertoire diversity
#'
#' @param input Single cell object or data.frame containing V(D)J data. If a
#' data.frame is provided, the cell barcodes should be stored as row names.
#' @param cluster_col meta.data column containing cluster IDs to use for
#' grouping cells when calculating clonotype abundance
#' @param method Method to use for calculating diversity. Can also pass a named
#' list to use multiple methods.
#' @param clonotype_col meta.data column containing clonotype IDs to use for
#' calculating clonotype abundance
#' @param plot_colors Character vector containing colors for plotting
#' @param plot_lvls Character vector containing levels for ordering
#' @param facet_rows The number of facet rows, use this argument if a list of
#' functions is passed to method
#' @param ... Additional arguments to pass to ggplot
#' @return ggplot object
#' @export
plot_diversity <- function(input, cluster_col = NULL, method = abdiv::simpson, clonotype_col = "clonotype_id",
                           plot_colors = NULL, plot_lvls = NULL, facet_rows = 1, ...) {

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
  plt_dat <- calc_diversity(
    input         = input,
    cluster_col   = cluster_col,
    method        = method,
    clonotype_col = clonotype_col,
    prefix        = "",
    return_df     = TRUE
  )

  plt_dat <- dplyr::filter(plt_dat, !is.na(!!sym(clonotype_col)))
  plt_dat <- dplyr::distinct(plt_dat, !!!syms(c(cluster_col, names(method))))

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


#' Plot repertoire overlap
#'
#' @param input Single cell object or data.frame containing V(D)J data. If a
#' data.frame is provided, the cell barcodes should be stored as row names.
#' @param cluster_col meta.data column containing cluster IDs to use for
#' calculating overlap
#' @param method Method to use for calculating similarity between clusters
#' @param clonotype_col meta.data column containing clonotype IDs to use for
#' calculating overlap
#' @param plot_colors Character vector containing colors for plotting
#' @param ... Additional arguments to pass to ggplot
#' @return ggplot object
#' @export
plot_similarity <- function(input, cluster_col, method = abdiv::jaccard, clonotype_col = "clonotype_id",
                            plot_colors = NULL, ...) {

  # Calculate similarity
  plt_dat <- calc_similarity(
    input         = input,
    cluster_col   = cluster_col,
    method        = method,
    clonotype_col = clonotype_col,
    prefix        = "",
    return_mat    = TRUE
  )

  sim_col <- as.character(substitute(method))
  sim_col <- dplyr::last(sim_col)

  var_lvls <- unique(c(rownames(plt_dat), colnames(plt_dat)))
  var_lvls <- sort(var_lvls)

  plt_dat <- tibble::as_tibble(plt_dat, rownames = "Var1")

  plt_dat <- tidyr::pivot_longer(
    plt_dat,
    cols      = -.data$Var1,
    names_to  = "Var2",
    values_to = sim_col
  )

  # Set Var levels
  plt_dat <- dplyr::mutate(
    plt_dat,
    Var1 = factor(.data$Var1, levels = rev(var_lvls)),
    Var2 = factor(.data$Var2, levels = var_lvls)
  )

  # Create heatmap
  res <- .create_heatmap(
    plt_dat,
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
#' @param input Object containing V(D)J data. If a data.frame is provided, the
#' cell barcodes should be stored as row names.
#' @param gene_cols meta.data column containing genes for each clonotype,
#' provide a vector with two column names to plot paired usage of genes
#' @param cluster_col meta.data column containing cell clusters to use for
#' calculating gene usage
#' @param chain Chain to use for calculating gene usage, set to NULL to include
#' all chains
#' @param type Type of plot to create, can be either 'heatmap' or 'bar', a
#' bar plot can only be created when a single gene_col is provided
#' @param plot_colors Character vector containing colors for plotting
#' @param plot_genes Character vector containing genes to plot
#' @param n_genes Number of top genes to plot based on usage. If cluster_col is
#' provided, top genes will be identified for each cluster.
#' @param plot_lvls Levels to use for ordering clusters
#' @param yaxis Units to plot on the y-axis, either 'frequency' or 'percent'
#' @param chain_col meta.data column containing chains for each cell
#' @param sep Separator used for storing per cell V(D)J data
#' @param ... Additional arguments to pass to ggplot
#' @return ggplot object
#' @export
plot_usage <- function(input, gene_cols, cluster_col = NULL, chain = NULL, type = "bar",
                       plot_colors = NULL, plot_genes = NULL, n_genes = 50, plot_lvls = NULL,
                       yaxis = "percent", chain_col = "chains", sep = ";", ...) {

  # Check inputs
  if (!type %in% c("heatmap", "bar")) {
    stop("type must be either 'heatmap' or 'bar'")
  }

  if (length(gene_cols) > 2) {
    stop("Cannot specify more than two values for gene_cols")
  }

  if (!yaxis %in% c("percent", "frequency")) {
    stop("yaxis must be either 'percent' or 'frequency'")
  }

  # Set y-axis
  usage_col <- "pct"

  if (yaxis == "frequency") {
    usage_col <- "freq"
  }

  # Calculate gene usage
  plt_dat <- calc_usage(
    input       = input,
    gene_cols   = gene_cols,
    cluster_col = cluster_col,
    chain       = chain,
    chain_col   = chain_col,
    sep         = sep
  )

  plt_dat <- dplyr::filter(plt_dat, dplyr::across(
    dplyr::all_of(gene_cols),
    ~ .x != "None"
  ))

  # If plot_genes provided check plt_dat
  top_genes <- unique(plot_genes)

  if (!is.null(top_genes)) {
    absent <- unlist(plt_dat[, gene_cols], use.names = FALSE)
    absent <- top_genes[!top_genes %in% absent]

    if (identical(absent, top_genes)) {
      stop("None of the provided genes were found")

    } else if (!is_empty(absent)) {
      absent <- paste0(absent, collapse = ", ")

      warning("Some genes not found: ", absent)
    }

  # If plot_genes not provided, identify top n_genes for each cluster based on
  # usage
  } else {
    top_genes <- plt_dat

    if (!is.null(cluster_col)) {
      top_genes <- dplyr::group_by(top_genes, !!sym(cluster_col))
    }

    # If two gene_cols are provided this will select the top gene pairs
    top_genes <- dplyr::slice_max(
      top_genes,
      order_by  = !!sym(usage_col),
      n         = n_genes,
      with_ties = FALSE
    )

    # Do not use dplyr::pull since user can pass multiple gene_cols
    top_genes <- unlist(top_genes[, gene_cols], use.names = FALSE)
    top_genes <- unique(top_genes)
  }

  # Filter for top genes
  # if two gene_cols are provided, all gene pairs containing at least one gene
  # from the top gene pairs will be included
  plt_dat <- dplyr::filter(plt_dat, dplyr::across(
    dplyr::all_of(gene_cols),
    ~ .x %in% top_genes
  ))

  # Order top genes based on max usage for clusters
  if (length(gene_cols) == 1) {
    lvls <- dplyr::group_by(plt_dat, !!!syms(gene_cols))

    lvls <- dplyr::summarize(
      lvls,
      max_usage = max(!!sym(usage_col)),
      .groups = "drop"
    )

    lvls <- dplyr::arrange(lvls, .data$max_usage)
    lvls <- pull(lvls, gene_cols)

    plt_dat <- .set_lvls(plt_dat, gene_cols, lvls)
  }

  # Order clusters based on plot_lvls
  plt_dat <- .set_lvls(plt_dat, cluster_col, plot_lvls)

  # Create heatmap or bar graph for single gene
  if (length(gene_cols) == 1) {
    if (type == "bar") {
      plt_dat <- dplyr::filter(plt_dat, !!sym(usage_col) > 0)

      res <- .create_bars(
        plt_dat,
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
      plt_dat,
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
    plt_dat <- dplyr::group_by(plt_dat, !!sym(cluster_col))

    grps <- pull(plt_dat, cluster_col)
    grps <- sort(unique(grps))
  }

  res <- dplyr::group_split(plt_dat)
  names(res) <- grps

  # Format data for plotting
  res <- map(res, ~ {
    res <- dplyr::select(.x, dplyr::all_of(c(gene_cols, usage_col)))

    res <- tidyr::pivot_wider(
      res,
      names_from  = gene_cols[2],
      values_from = all_of(usage_col)
    )

    res <- dplyr::mutate(
      res,
      dplyr::across(-dplyr::all_of(gene_cols[1]), ~ tidyr::replace_na(.x, 0))
    )

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
  res <- map(
    res,
    .create_heatmap,
    x     = gene_cols[1],
    y     = gene_cols[2],
    .fill = usage_col,
    clrs  = plot_colors,
    ttl   = yaxis,
    ...
  )

  if (!is.null(names(res))) {
    res <- imap(res, ~ .x + ggplot2::ggtitle(.y))
  }

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
djvdj_theme <- function(ttl_size = 12, txt_size = 8, ln_size = 0.5, txt_col = "black",
                        ln_col = "grey85") {
  res <- ggplot2::theme(
    strip.background  = ggplot2::element_blank(),
    strip.text        = ggplot2::element_text(size = ttl_size),
    panel.background  = ggplot2::element_blank(),
    legend.background = ggplot2::element_blank(),
    legend.title      = ggplot2::element_text(size = ttl_size),
    legend.key        = ggplot2::element_blank(),
    legend.text       = ggplot2::element_text(size = txt_size, color = txt_col),
    axis.line         = ggplot2::element_line(size = ln_size,  color = ln_col),
    axis.ticks        = ggplot2::element_line(size = ln_size,  color = ln_col),
    axis.text         = ggplot2::element_text(size = txt_size, color = txt_col),
    axis.title        = ggplot2::element_text(size = ttl_size, color = txt_col)
  )

  res
}


#' Set min and max values for column
#'
#' @param df_in Input data.frame
#' @param ft Name of column containing feature values
#' @param mn Minimum quantile
#' @param mx Maximum quantile
#' @return data.frame with modified feature values
#' @noRd
.set_lims <- function(df_in, ft, mn = NULL, mx = NULL) {

  if (!is.numeric(df_in[[ft]])) {
    stop(ft, " is not numeric")
  }

  if (is.null(mn) && is.null(mx)) {
    stop("mn or mx must be provided")
  }

  ft <- sym(ft)

  res <- mutate(df_in, pct = dplyr::percent_rank(!!ft))

  if (!is.null(mn)) {
    mn <- mn / 100

    res <- dplyr::mutate(
      res,
      !!ft := ifelse(.data$pct > mn, !!ft, NA),
      !!ft := ifelse(.data$pct <= mn, min(!!ft, na.rm = TRUE), !!ft)
    )
  }

  if (!is.null(mx)) {
    mx <- mx / 100

    res <- dplyr::mutate(
      res,
      !!ft := ifelse(.data$pct < mx, !!ft, NA),
      !!ft := ifelse(.data$pct >= mx, max(!!ft, na.rm = TRUE), !!ft)
    )
  }

  res <- dplyr::select(res, -.data$pct)

  res
}


#' Add list of aes params to ggplot object
#'
#' @param gg_in ggplot object
#' @param add_aes List of aes params to add or override
#' @param lay Layer number to modify
#' @return ggplot object
#' @noRd
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
#' @noRd
.create_heatmap <- function(df_in, x, y, .fill, clrs = NULL, na_color = "white",
                            ttl = .fill, ang = 45, hjst = 1, ...) {

  if (is.null(clrs)) {
    clrs <- "#6A51A3"
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
    djvdj_theme() +
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
#' @noRd
.create_bars <- function(df_in, x, y, .fill, clrs = NULL, y_ttl = y, ang = 45,
                         hjst = 1, ...) {

  # Bar position to use for plot
  bar_pos <- ggplot2::position_dodge(preserve = "single")

  # Reverse bar order
  lvls  <- rev(levels(pull(df_in, x)))
  df_in <- .set_lvls(df_in, x, lvls)

  # Color fill by variable
  if (!is.null(.fill)) {
    res <- ggplot2::ggplot(
      df_in,
      ggplot2::aes(!!sym(x), !!sym(y), fill = !!sym(.fill))
    ) +
    ggplot2::geom_col(..., position = bar_pos)

    if (!is.null(clrs)) {
      res <- res +
        ggplot2::scale_fill_manual(values = clrs)
    }

  # Use single fill color
  } else {
    res <- ggplot2::ggplot(
      df_in,
      ggplot2::aes(!!sym(x), !!sym(y))
    )

    if (!is.null(clrs)) {
      res <- res +
        ggplot2::geom_col(..., fill = clrs, position = bar_pos)

    } else {
      res <- res +
        ggplot2::geom_col(..., position = bar_pos)
    }
  }

  res <- res +
    ggplot2::labs(y = y_ttl) +
    djvdj_theme() +
    ggplot2::theme(
      axis.title.x = ggplot2::element_blank(),
      axis.text.x  = ggplot2::element_text(angle = ang, hjust = hjst)
    )

  res
}


#' Set column levels
#'
#' @param df_in data.frame
#' @param clmn Column to modify
#' @param lvls Levels
#' @return data.frame
#' @noRd
.set_lvls <- function(df_in, clmn, lvls) {

  if (!is.null(lvls) && !is.null(clmn)) {
    dat <- pull(df_in, clmn)

    if (!is.character(dat) && !is.factor(dat)) {
      warning("Plot levels were not modified, levels are only modified for characters and factors")

      return(df_in)
    }

    if (!all(pull(df_in, clmn) %in% lvls)) {
      stop("Not all labels in ", clmn, " are included in plot_levels")
    }

    df_in <- dplyr::mutate(
      df_in,
      !!sym(clmn) := factor(!!sym(clmn), levels = unique(lvls))
    )
  }

  df_in
}


#' Plot CDR3 lengths
#'
#' @param input Single cell object or data.frame containing V(D)J data. If a
#' data.frame is provided, the cell barcodes should be stored as row names.
#' @param length_col meta.data column containing CDR3 lengths to plot
#' @param cluster_col meta.data column containing cluster IDs to use for
#' grouping cells when plotting CDR3 lengths
#' @param chain Chain to use for plotting CDR3 lengths, set to NULL to include
#' all chains
#' @param type Type of plot to create, can be 'histogram' or 'violin'
#' @param yaxis Units to use for y-axis when plotting histogram. Use
#' 'frequency' to show the raw number of CDR3 sequences or 'percent' to show
#' the percentage of total CDR3 sequences.
#' @param plot_colors Character vector containing colors for plotting
#' @param plot_lvls Character vector containing levels for ordering
#' @param chain_col meta.data column containing chains for each cell
#' @param sep Separator used for storing per cell V(D)J data
#' @param ... Additional arguments to pass to ggplot
#' @return ggplot object
#' @export
#' @importFrom stats median
plot_cdr3_length <- function(input, length_col = "cdr3_length", cluster_col = NULL, chain = NULL, type = "histogram",
                             yaxis = "frequency", plot_colors = NULL, plot_lvls = NULL, chain_col = "chains",
                             sep = ";", ...) {

  # Check inputs
  if (!type %in% c("histogram", "density", "violin")) {
    stop("type must be either 'histogram', 'density', or 'violin'")
  }

  if (length(length_col) != 1) {
    stop("Must specify only one value for length_col")
  }

  if (!yaxis %in% c("percent", "frequency")) {
    stop("yaxis must be either 'percent' or 'frequency'")
  }

  # Format input data
  sep_cols <- c(length_col, chain_col)
  vdj_cols <- c(".cell_id", cluster_col, sep_cols)

  meta <- .get_meta(input)
  meta <- dplyr::select(meta, all_of(vdj_cols))

  meta <- dplyr::filter(meta, !is.na(!!sym(length_col)))

  # Expand meta.data based on sep
  coerce_cols <- purrr::set_names(c("numeric", "character"), sep_cols)

  plt_dat <- .split_vdj(
    meta,
    sep         = sep,
    sep_cols    = sep_cols,
    coerce_cols = coerce_cols,
    expand      = TRUE
  )

  # Filter for chain
  if (!is.null(chain)) {
    plt_dat <- dplyr::filter(plt_dat, !!sym(chain_col) %in% chain)
  }

  # Order clusters based on plot_lvls
  plt_dat <- .set_lvls(plt_dat, cluster_col, plot_lvls)

  # Create violin plot
  if (type == "violin") {
    if (!is.null(cluster_col)) {
      plt_aes <- ggplot2::aes(
        !!sym(cluster_col),
        !!sym(length_col),
        color = !!sym(cluster_col),
        fill  = !!sym(cluster_col)
      )

    } else {
      plt_aes <- ggplot2::aes(length_col, !!sym(length_col))
    }

    res <- ggplot2::ggplot(plt_dat, plt_aes) +
      ggplot2::geom_violin(...) +
      ggplot2::stat_summary(geom = "point", fun = stats::median, color = "black")

    return(res)
  }

  # Create histogram
  y_lab <- "number of"

  if (!is.null(cluster_col)) {
    plt_aes <- ggplot2::aes(
      !!sym(length_col),
      color = !!sym(cluster_col),
      fill  = !!sym(cluster_col)
    )

    if (yaxis == "percent" && type != "density") {
      y_lab <- "percentage of"

      plt_aes <- ggplot2::aes(
        !!sym(length_col),
        .data$..count.. / sum(.data$..count..) * 100,
        color = !!sym(cluster_col),
        fill  = !!sym(cluster_col)
      )
    }

  } else {
    plt_aes <- ggplot2::aes(!!sym(length_col))

    if (yaxis == "percent" && type != "density") {
      y_lab <- "percentage of"

      plt_aes <- ggplot2::aes(
        !!sym(length_col),
        .data$..count.. / sum(.data$..count..) * 100
      )
    }
  }

  res <- ggplot2::ggplot(plt_dat, plt_aes)

  # Set plot colors
  if (!is.null(plot_colors)) {
    res <- res +
      ggplot2::scale_color_manual(values = plot_colors) +
      ggplot2::scale_fill_manual(values = plot_colors)
  }

  if (type == "density") {
    res <- res +
      ggplot2::geom_density(..., alpha = 0.5)

    return(res)
  }

  res <- res +
    ggplot2::geom_histogram(...) +
    ggplot2::labs(y = paste(y_lab, "CDR3 sequences"))

  res
}

