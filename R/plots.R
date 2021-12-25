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
#' @param method Function to use for calculating diversity. A named list of
#' functions can be passed to plot multiple diversity metrics,
#' e.g. list(simpson = abdiv::simpson, shannon = abdiv::shannon)
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


#' Create ggplot boxplot
#'
#' @param df_in data.frame
#' @param x Variable to plot on x-axis
#' @param y Variable to plot on y-axis
#' @param .color Variable to use for line color
#' @param .fill Variable to use for fill color
#' @param clrs Vector of colors for plotting
#' @param type Type of plot to create, either 'boxplot' or 'violin'
#' @param log_trans If TRUE, log10 transform the y-axis
#' @param ... Additional arguments to pass to ggplot
#' @return ggplot object
#' @noRd
.create_boxes <- function(df_in, x = NULL, y, .color = NULL, .fill = NULL, clrs = NULL,
                          type = "boxplot", log_trans = FALSE, ...) {

  # Check input
  typs <- c("boxplot", "violin")

  if (!type %in% typs) {
    stop("type must be one of ", paste0(typs, collapse = ", "))
  }

  # Set aesthetics
  plt_aes <- ggplot2::aes(y, !!sym(y))

  if (!is.null(x)) {
    plt_aes$x <- sym(x)
  }

  if (!is.null(.fill)) {
    plt_aes$fill <- sym(.fill)
  }

  if (!is.null(.color)) {
    plt_aes$colour <- sym(.color)
  }

  # Adjust theme
  res <- ggplot2::ggplot(df_in, plt_aes) +
    djvdj_theme() +
    ggplot2::theme(
      legend.position = "none",
      axis.title.x    = ggplot2::element_blank()
    )

  if (!is.null(clrs)) {
    res <- res +
      ggplot2::scale_color_manual(values = clrs) +
      ggplot2::scale_fill_manual(values = clrs)
  }

  # Log10 tranform y-axis
  if (log_trans) {
    res <- res +
      ggplot2::scale_y_log10()
  }

  # Return violins
  if (type == "violin") {
    res <- res +
      ggplot2::geom_violin(...) +
      ggplot2::stat_summary(geom = "point", fun = stats::median, color = "black")

    return(res)
  }

  # Return boxes
  res <- res +
    ggplot2::geom_boxplot(...)

  res
}


#' Create ggplot histogram
#'
#' @param df_in data.frame
#' @param x Variable to plot on x-axis
#' @param .color Variable to use for line color
#' @param .fill Variable to use for fill color
#' @param clrs Vector of colors for plotting
#' @param type Type of plot to create, either 'histogram' or 'density'
#' @param yaxis Units to use for y-axis when plotting histogram. Use
#' 'frequency' to show the raw number of values or 'percent' to show
#' the percentage of total values.
#' @param log_trans If TRUE, log10 transform the x-axis
#' @param ... Additional arguments to pass to ggplot
#' @return ggplot object
#' @noRd
.create_hist <- function(df_in, x, .color = NULL, .fill = NULL, clrs = NULL, type = "histogram",
                         yaxis = "frequency", log_trans = FALSE, ...) {

  # Check inputs
  typs <- c("histogram", "density")

  if (!type %in% typs) {
    stop("type must be one of ", paste0(typs, collapse = ", "))
  }

  axs <- c("frequency", "percent")

  if (!yaxis %in% axs) {
    stop("yaxis must be one of ", paste0(axs, collapse = ", "))
  }

  # Set aesthetics
  plt_aes <- ggplot2::aes()

  # Only plot percent for histogram
  if (yaxis == "percent" && type == "histogram") {
    plt_aes <- ggplot2::aes(y = .data$..count.. / sum(.data$..count..) * 100)
  }

  plt_aes$x <- sym(x)

  if (!is.null(.fill)) {
    plt_aes$fill <- sym(.fill)
  }

  if (!is.null(.color)) {
    plt_aes$colour <- sym(.color)
  }

  # Adjust theme
  res <- ggplot2::ggplot(df_in, plt_aes) +
    djvdj_theme() +
    ggplot2::theme(legend.position = "top")

  if (!is.null(clrs)) {
    res <- res +
      ggplot2::scale_color_manual(values = clrs) +
      ggplot2::scale_fill_manual(values = clrs)
  }

  # Log10 transform axis
  if (log_trans) {
    res <- res +
      ggplot2::scale_x_log10()
  }

  # Return density plot
  if (type == "density") {
    res <- res +
      ggplot2::geom_density(...)

    return(res)
  }

  # Return histogram
  res <- res +
    ggplot2::geom_histogram(...) +
    ggplot2::labs(y = yaxis)

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


#' Plot continuous V(D)J data
#'
#' @param input Single cell object or data.frame containing V(D)J data. If a
#' data.frame is provided, cell barcodes should be stored as row names.
#' @param data_cols meta.data column(s) containing continuous V(D)J data to
#' plot
#' @param cluster_col meta.data column containing cluster IDs to use for
#' grouping cells when plotting CDR3 lengths
#' @param chain Chain(s) to use for plotting, set to NULL to include all chains
#' @param type Type of plot to create, can be 'histogram', 'density',
#' 'boxplot', or 'violin'
#' @param yaxis Units to use for y-axis when plotting histogram. Use
#' 'frequency' to show the raw number of values or 'percent' to show
#' the percentage of total values.
#' @param plot_colors Character vector containing colors for plotting
#' @param plot_lvls Character vector containing levels for ordering
#' @param chain_col meta.data column containing chains for each cell
#' @param sep Separator used for storing per cell V(D)J data
#' @param ... Additional arguments to pass to ggplot geom
#' @return ggplot object
#' @name plot_vdj
NULL

#' @rdname plot_vdj
#' @param per_cell Should values be plotted per cell, i.e. each data point
#' would represent one cell. If TRUE, values will be summarized for each cell
#' using summary_fn. If FALSE, values will be plotted per chain.
#' @param summary_fn Function to use for summarizing values for each cell, e.g.
#' specifying stats::median will plot the median value per cell
#' @param alpha Color transparency
#' @param log_trans If TRUE, axis will be log10 transformed
#' @export
plot_vdj <- function(input, data_cols, per_cell = FALSE, summary_fn = mean, cluster_col = NULL, chain = NULL,
                     type = "histogram", yaxis = "frequency", plot_colors = NULL, plot_lvls = NULL, alpha = 0.5,
                     log_trans = FALSE, chain_col = "chains", sep = ";", ...) {

  fn <- function(clmn) {
    .plot_vdj(
      input,
      data_col    = clmn,
      per_cell    = per_cell,
      cluster_col = cluster_col,
      chain       = chain,
      type        = type,
      yaxis       = yaxis,
      plot_colors = plot_colors,
      plot_lvls   = plot_lvls,
      log_trans   = log_trans,
      chain_col   = chain_col,
      sep         = sep,
      ...
    )
  }

  res <- purrr::map(data_cols, fn)

  # Combine plots
  res <- purrr::reduce(res, `+`)

  res
}

#' @rdname plot_vdj
#' @noRd
.plot_vdj <- function(input, data_cols, per_cell = FALSE, summary_fn = mean, cluster_col = NULL, chain = NULL,
                      type = "boxplot", yaxis = "frequency", plot_colors = NULL, plot_lvls = NULL, alpha = 0.5,
                      log_trans = FALSE, chain_col = "chains", sep = ";", ...) {

  # Check inputs
  if (length(data_col) != 1) {
    stop("Must specify only one value for data_col")
  }

  # Format input data
  if (per_cell) {
    plt_dat <- summarize_vdj(
      input,
      vdj_cols  = data_col,
      fn        = summary_fn,
      chain     = chain,
      chain_col = chain_col,
      sep       = sep,
      return_df = TRUE
    )

  } else {
    fetch_cols <- data_col

    if (!is.null(chain)) {
      fetch_cols <- c(chain_col, fetch_cols)
    }

    plt_dat <- fetch_vdj(
      input,
      vdj_cols      = fetch_cols,
      clonotype_col = clonotype_col,
      unnest        = TRUE
    )

    if (!is.null(chain)) {
      plt_dat <- dplyr::filter(plt_dat, !!sym(chain_col) %in% chain)
    }
  }

  plt_dat <- dplyr::filter(plt_dat, !is.na(!!sym(data_col)))

  # Order clusters based on plot_lvls
  plt_dat <- .set_lvls(plt_dat, cluster_col, plot_lvls)

  # Create violin plot
  if (type %in% c("histogram", "density")) {
    res <- .create_hist(
      plt_dat,
      x         = data_col,
      .color    = cluster_col,
      .fill     = cluster_col,
      clrs      = plot_colors,
      type      = type,
      yaxis     = yaxis,
      log_trans = log_trans,
      alpha     = alpha,
      ...
    )

    return(res)
  }

  res <- .create_boxes(
    plt_dat,
    x         = cluster_col,
    y         = data_col,
    .color    = cluster_col,
    .fill     = cluster_col,
    clrs      = plot_colors,
    type      = type,
    log_trans = log_trans,
    alpha     = alpha,
    ...
  )

  res
}

#' @rdname plot_vdj
#' @export
plot_vdj_reads <- function(input, data_cols = c("reads", "umis"), cluster_col = NULL, chain = NULL,
                           type = "violin", yaxis = "frequency", plot_colors = NULL, plot_lvls = NULL,
                           chain_col = "chains", sep = ";", ...) {

  res <- plot_vdj(
    input,
    data_cols   = data_cols,
    per_cell    = FALSE,
    cluster_col = cluster_col,
    chain       = chain,
    type        = type,
    yaxis       = yaxis,
    plot_colors = plot_colors,
    plot_lvls   = plot_lvls,
    log_trans   = TRUE,
    chain_col   = chain_col,
    sep         = sep,
    ...
  )

  res
}

#' @rdname plot_vdj
#' @export
plot_cdr3_length <- function(input, data_cols = "cdr3_length", cluster_col = NULL, chain = NULL,
                             type = "histogram", yaxis = "frequency", plot_colors = NULL, plot_lvls = NULL,
                             chain_col = "chains", sep = ";", ...) {

  res <- plot_vdj(
    input,
    data_cols   = data_cols,
    per_cell    = FALSE,
    cluster_col = cluster_col,
    chain       = chain,
    type        = type,
    yaxis       = yaxis,
    plot_colors = plot_colors,
    plot_lvls   = plot_lvls,
    log_trans   = TRUE,
    chain_col   = chain_col,
    sep         = sep,
    ...
  )

  res
}

#' @rdname plot_vdj
#' @export
plot_vdj_indels <- function(input, data_cols = c("n_insertion", "n_deletion", "n_mismatch"), cluster_col = NULL, chain = NULL,
                            type = "violin", yaxis = "frequency", plot_colors = NULL, plot_lvls = NULL,
                            chain_col = "chains", sep = ";", ...) {

  res <- plot_vdj(
    input,
    data_cols   = data_cols,
    per_cell    = FALSE,
    cluster_col = cluster_col,
    chain       = chain,
    type        = type,
    yaxis       = yaxis,
    plot_colors = plot_colors,
    plot_lvls   = plot_lvls,
    log_trans   = FALSE,
    chain_col   = chain_col,
    sep         = sep,
    ...
  )

  res
}

