#' Set min and max values for column
#'
#' @param df_in Input data.frame
#' @param feat_col Name of column containing feature values
#' @param lim The value cutoff
#' @param op The operator to use for comparing values with lim
#' (either "<" or ">")
#' @return data.frame with modified feature values
set_val_limits <- function(df_in, feat_col, lim, op) {

  if (!op %in% c("<", ">")) {
    stop("op must be either \"<\" or \">\"")
  }

  func <- "min"

  if (op == "<") {
    func <- "max"
  }

  res <- dplyr::mutate(
    df_in,
    pct_rank = dplyr::percent_rank(!!dplyr::sym(feat_col)),

    lim = eval(parse(text = paste0(
      "ifelse(pct_rank ", op, " lim, ", feat_col, ", NA)"
    ))),

    lim = eval(parse(text = paste0(
      func, "(lim, na.rm = T)"
    ))),

    !!dplyr::sym(feat_col) := eval(parse(text = paste0(
      "dplyr::if_else(", feat_col, op, " lim, lim, ", feat_col, ")"
    )))
  )

  res <- dplyr::select(res, -.data$pct_rank, -lim)

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
add_lm_line <- function(gg_in, lab_pos = NULL, lab_size = 3.5, ...) {

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
    gg_data <- gg_in$data

    x <- as.character(gg_in$mapping$x)[2]
    y <- as.character(gg_in$mapping$y)[2]

    gg_data <- dplyr::mutate(
      gg_data,
      r       = broom::tidy(stats::cor.test(!!dplyr::sym(x), !!dplyr::sym(y)))$estimate,
      r       = round(.data$r, digits = 2),
      pval    = broom::tidy(stats::cor.test(!!dplyr::sym(x), !!dplyr::sym(y)))$p.value,
      cor_lab = stringr::str_c("r = ", .data$r, ", p = ", format(.data$pval, digits = 2)),
      min_x   = min(!!dplyr::sym(x)),
      max_x   = max(!!dplyr::sym(x)),
      min_y   = min(!!dplyr::sym(y)),
      max_y   = max(!!dplyr::sym(y)),
      lab_x   = (.data$max_x - .data$min_x) * lab_pos[1] + .data$min_x,
      lab_y   = (.data$max_y - .data$min_y) * lab_pos[2] + .data$min_y
    )
  }

  # Add correlation coefficient label
  res <- res +
    ggplot2::geom_text(
      data          = gg_data,
      mapping       = ggplot2::aes(.data$lab_x, .data$lab_y, label = .data$cor_lab),
      color         = "black",
      size          = lab_size,
      check_overlap = T,
      show.legend   = F
    )

  res
}


#' Create two dimensional scatter plot
#'
#' @param obj_in Seurat object or data.frame containing data for plotting
#' @param x Variable to plot on x-axis
#' @param y Variable to plot on y-axis
#' @param feature Variable to use for coloring points
#' @param data_slot Slot in the Seurat object to pull data from
#' @param pt_size Point size
#' @param plot_cols Vector of colors to use for plot
#' @param feat_levels Levels to use for ordering feature
#' @param split_id Variable to use for splitting plot into facets
#' @param split_levels Levels to use for ordering plot facets
#' @param min_pct Minimum value to plot for feature
#' @param max_pct Maximum value to plot for feature
#' @param lm_line Add regression line to plot
#' @param cor_label Position of correlation coefficient label
#' @param label_size Size of correlation coefficient label
#' @param ... Additional parameters to pass to facet_wrap
#' @return ggplot object
#' @export
plot_features <- function(obj_in, x = "UMAP_1", y = "UMAP_2", feature,
                          data_slot = "data", pt_size = 0.25, plot_cols = NULL,
                          feat_levels = NULL, split_id = NULL,
                          split_levels = NULL, min_pct = NULL, max_pct = NULL,
                          lm_line = F, cor_label = c(0.8, 0.9),
                          label_size = 3.7, ...) {

  # Format imput data
  meta_df <- obj_in

  if ("Seurat" %in% class(obj_in)) {
    vars <- c(x, y, feature)

    if (!is.null(split_id)) {
      vars <- c(vars, split_id)
    }

    meta_df <- Seurat::FetchData(obj_in, vars = unique(vars), slot = data_slot)
    meta_df <- tibble::as_tibble(meta_df, rownames = ".cell_id")

    if (!feature %in% colnames(meta_df)) {
      stop(paste(feature, "not found in object."))
    }
  }

  # Rename features
  if (!is.null(names(feature))) {
    meta_df  <- dplyr::rename(meta_df, !!!dplyr::syms(feature))
    feature <- names(feature)
  }

  if (!is.null(names(x))) {
    meta_df <- dplyr::rename(meta_df, !!!dplyr::syms(x))
    x <- names(x)
  }

  if (!is.null(names(y))) {
    meta_df <- dplyr::rename(meta_df, !!!dplyr::syms(y))
    y <- names(y)
  }

  # Adjust values based on min_pct and max_pct
  if (!is.null(min_pct)) {
    meta_df <- set_val_limits(
      meta_df,
      feat_col = feature,
      lim      = min_pct,
      op       = "<"
    )
  }

  if (!is.null(max_pct)) {
    meta_df <- set_val_limits(
      meta_df,
      feat_col = feature,
      lim      = max_pct,
      op       = ">"
    )
  }

  # Set feature order
  if (!is.null(feat_levels)) {
    meta_df <- dplyr::mutate(
      meta_df,
      !!dplyr::sym(feature) := factor(
        !!dplyr::sym(feature),
        levels = feat_levels
      )
    )
  }

  # Set facet order
  if (!is.null(split_id) && length(split_id) == 1 && !is.null(split_levels)) {
    meta_df <- dplyr::mutate(
      meta_df,
      !!dplyr::sym(split_id) := factor(
        !!dplyr::sym(split_id),
        levels = split_levels
      )
    )
  }

  # Create scatter plot
  meta_df <- dplyr::arrange(meta_df, !!dplyr::sym(feature))

  res <- ggplot2::ggplot(meta_df, ggplot2::aes(
    !!dplyr::sym(x),
    !!dplyr::sym(y),
    color = !!dplyr::sym(feature)
  )) +
    ggplot2::geom_point(size = pt_size)

  # Add regression line
  if (lm_line) {
    res <- add_lm_line(
      res,
      lab_pos = cor_label,
      lab_size = label_size
    )
  }

  # Set feature colors
  if (!is.null(plot_cols)) {
    if (is.numeric(meta_df[[feature]])) {
      res <- res +
        ggplot2::scale_color_gradientn(colors = plot_cols)

    } else {
      res <- res +
        ggplot2::scale_color_manual(values = plot_cols)
    }
  }

  # Split plot into facets
  if (!is.null(split_id)) {
    if (length(split_id) == 1) {
      res <- res +
        ggplot2::facet_wrap(~ !!dplyr::sym(split_id), ...)

    } else if (length(split_id) == 2) {
      eq <- stringr::str_c(split_id[1], " ~ ", split_id[2])

      res <- res +
        ggplot2::facet_grid(stats::as.formula(eq), ...)
    }
  }

  res
}


#' Plot clonotype abundance
#'
#' @param sobj_in Seurat object containing VDJ data
#' @param clonotype_col meta.data column containing clonotype IDs to use for
#' calculating clonotype abundance
#' @param cluster_col meta.data column containing cluster IDs to use for
#' grouping cells when calculating clonotype abundance
#' @param yaxis Units to plot on the y-axis, either "frequency" or "percent"
#' @param abundance_col meta.data column containing pre-calculated abundances
#' @param plot_colors Character vector containing colors to use for plotting
#' @param plot_levels Character vector containing levels to use for ordering
#' @param ... Additional arguments to pass to geom_line
#' @return Seurat object with calculated clonotype abundance
#' @export
plot_abundance <- function(sobj_in, clonotype_col = "clonotype_id", cluster_col = NULL,
                           yaxis = "percent", abundance_col = NULL, plot_colors = NULL, plot_levels = NULL,
                           ...) {

  # Calculate clonotype abundance
  if (is.null(abundance_col)) {
    sobj_in <- djvdj::calc_abundance(
      sobj_in,
      clonotype_col = clonotype_col,
      cluster_col   = cluster_col,
      prefix        = "."
    )

    abundance_col <- ".clone_pct"

    if (yaxis == "frequency") {
      abundance_col <- ".clone_freq"
    }
  }

  meta_df <- sobj_in@meta.data
  meta_df <- tibble::as_tibble(meta_df, rownames = ".cell_id")
  meta_df <- dplyr::filter(meta_df, !is.na(!!dplyr::sym(clonotype_col)))

  # Rank by abundance
  if (!is.null(cluster_col)) {
    if (!is.null(plot_levels)) {
      meta_df <- dplyr::mutate(
        meta_df,
        !!dplyr::sym(cluster_col) := factor(
          !!dplyr::sym(cluster_col),
          levels = plot_levels
        )
      )
    }

    meta_df <- dplyr::group_by(meta_df, !!dplyr::sym(cluster_col))
  }

  meta_df <- dplyr::mutate(
    meta_df,
    rank = dplyr::row_number(dplyr::desc(!!dplyr::sym(abundance_col)))
  )

  # Plot abundance vs rank
  gg <- ggplot2::ggplot(meta_df, ggplot2::aes(rank, !!dplyr::sym(abundance_col))) +
    ggplot2::labs(y = yaxis)

  if (is.null(cluster_col)) {
    gg <- gg +
      ggplot2::geom_line(...)

  } else {
    gg <- gg +
      ggplot2::geom_line(ggplot2::aes(color = !!dplyr::sym(cluster_col)), ...)
  }

  if (!is.null(plot_colors)) {
    gg <- gg +
      ggplot2::scale_color_manual(values = plot_colors)
  }

  gg
}


#' Plot VDJ gene usage
#'
#' @param sobj_in Seurat object containing VDJ data
#' @param gene_col meta.data column containing genes used for each clonotype
#' @param cluster_col meta.data column containing cell clusters to use for
#' calculating gene usage
#' @param chain Chain to use for calculating gene usage. Set to NULL to include
#' all chains.
#' @param plot_colors Character vector containing colors to use for plotting
#' @param plot_genes Character vector of genes to plot
#' @param n_genes Number of top genes to plot based on average usage
#' @param clust_levels Levels to use for ordering clusters
#' @param yaxis Units to plot on the y-axis, either "frequency" or "percent"
#' @param vline_color Color of vertical line to separate clusters
#' @param ... Additional arguments to pass to geom_tile
#' @param chain_col meta.data column containing chains for each cell
#' @param sep Separator to use for expanding gene_col
#' @return ggplot object
#' @export
plot_usage <- function(sobj_in, gene_col, cluster_col = NULL, chain = NULL, plot_colors = NULL,
                       plot_genes = NULL, n_genes = NULL, clust_levels = NULL, yaxis = "percent",
                       vline_color = NULL, ..., chain_col = "chains", sep = ";") {

  # Calculate gene usage
  gg_data <- calc_usage(
    sobj_in,
    gene_col    = gene_col,
    cluster_col = cluster_col,
    chain       = chain,
    chain_col   = chain_col,
    sep         = sep
  )

  gg_data <- dplyr::filter(gg_data, !!dplyr::sym(gene_col) != "None")

  usage_col <- "pct"

  if (yaxis == "frequency") {
    usage_col <- "freq"
  }

  # Order genes by average usage
  top_genes <- dplyr::group_by(gg_data, !!dplyr::sym(gene_col))

  top_genes <- dplyr::summarize(
    top_genes,
    usage   = mean(!!dplyr::sym(usage_col)),
    .groups = "drop"
  )

  top_genes <- dplyr::arrange(top_genes, .data$usage)

  gg_data <- dplyr::mutate(
    gg_data,
    !!dplyr::sym(gene_col) := factor(
      !!dplyr::sym(gene_col),
      levels = dplyr::pull(top_genes, gene_col)
    )
  )

  # Order clusters
  if (!is.null(clust_levels)) {
    gg_data <- dplyr::mutate(
      gg_data,
      !!dplyr::sym(cluster_col) := factor(
        !!dplyr::sym(cluster_col),
        levels = clust_levels
      )
    )
  }

  # Filter genes to plot
  if (!is.null(plot_genes)) {
    gg_data <- dplyr::filter(gg_data, !!dplyr::sym(gene_col) %in% plot_genes)
  }

  # Select top genes to plot
  if (!is.null(n_genes)) {
    top_genes <- dplyr::slice_max(top_genes, .data$usage, n = n_genes)
    top_genes <- pull(top_genes, gene_col)

    gg_data <- dplyr::filter(gg_data, !!dplyr::sym(gene_col) %in% top_genes)
  }

  # Create heatmap
  res <- ggplot2::ggplot(gg_data, ggplot2::aes(
    !!dplyr::sym(cluster_col),
    !!dplyr::sym(gene_col),
    fill = !!dplyr::sym(usage_col)
  )) +
    ggplot2::geom_tile(...) +
    ggplot2::guides(fill = ggplot2::guide_colorbar(title = yaxis)) +
    ggplot2::theme(
      axis.title = ggplot2::element_blank(),
      axis.line  = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank()
    )

  # Add dividing line
  if (!is.null(cluster_col) && !is.null(vline_color)) {
    n_clusts <- dplyr::n_distinct(gg_data[, cluster_col])
    ival     <- 1
    xint     <- seq(ival + 0.5, n_clusts - 0.5, ival)

    res <- res +
      ggplot2::geom_vline(xintercept = xint, color = vline_color)
  }

  # Set colors
  if (!is.null(plot_colors)) {
    res <- res +
      ggplot2::scale_fill_gradientn(colors = plot_colors, name = yaxis)
  }

  res
}


#' Plot repertoire overlap
#'
#' @param obj_in Seurat object containing VDJ data or matrix
#' @param clonotype_col meta.data column containing clonotype IDs to use for
#' calculating overlap
#' @param cluster_col meta.data column containing cluster IDs to use for
#' calculating overlap
#' @param plot_colors Character vector containing colors to use for plotting
#' @param ... Additional arguments to pass to geom_tile
#' @return ggplot object
#' @export
plot_overlap <- function(obj_in, clonotype_col = NULL, cluster_col = NULL,
                         plot_colors = NULL, ...) {

  if ("Seurat" %in% class(obj_in)) {
    if (is.null(clonotype_col) || is.null(cluster_col)) {
      stop("if a Seurat object is provided, clonotype_col and cluster_col must be specified")
    }

    obj_in <- calc_overlap(
      sobj_in       = obj_in,
      clonotype_col = clonotype_col,
      cluster_col   = cluster_col,
      prefix        = "",
      return_seurat = FALSE
    )
  }

  if (!identical(rownames(obj_in), colnames(obj_in))) {
    stop("matrix must have the same column and row names")
  }

  gg_data <- tibble::as_tibble(obj_in, rownames = "Var1")

  gg_data <- tidyr::pivot_longer(
    gg_data,
    cols      = -Var1,
    names_to  = "Var2",
    values_to = "jaccard"
  )

  gg_data <- dplyr::mutate(gg_data, jaccard = ifelse(Var1 == Var2, NA, jaccard))
  gg_data <- dplyr::rowwise(gg_data)

  gg_data <- dplyr::mutate(
    gg_data,
    key = stringr::str_c(sort(c(Var1, Var2)), collapse = "_")
  )

  # Add NAs so each comparison is only included once
  gg_data <- dplyr::group_by(gg_data, key)
  gg_data <- dplyr::mutate(gg_data, jaccard = ifelse(dplyr::row_number() == 2, NA, jaccard))
  gg_data <- dplyr::ungroup(gg_data)

  # Set Var levels
  var_levels <- unique(gg_data$Var1)

  gg_data <- dplyr::mutate(
    gg_data,
    Var1 = factor(Var1, levels = var_levels),
    Var2 = factor(Var2, levels = rev(var_levels))
  )

  gg_data <- dplyr::filter(
    gg_data,
    Var1 != levels(Var2)[1],
    Var2 != levels(Var1)[1]
  )

  # Create heatmap
  res <- ggplot2::ggplot(gg_data, aes(Var1, Var2, fill = jaccard)) +
    ggplot2::geom_tile(...) +
    ggplot2::theme(
      axis.title  = element_blank(),
      axis.line   = element_blank(),
      axis.ticks  = element_blank()
    )

  if (!is.null(plot_colors)) {
    res <- res +
      ggplot2::scale_fill_gradientn(colors = plot_colors, na.value = "white")
  }

  res
}







