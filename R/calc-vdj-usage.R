#' Calculate V(D)J gene segment usage
#'
#' @param input Object containing V(D)J data. If a data.frame is provided, the
#' cell barcodes should be stored as row names.
#' @param gene_cols meta.data column containing V(D)J genes identified for each
#' clonotype. If multiple columns are provided, paired usage of genes will be
#' calculated.
#' @param cluster_col meta.data column containing cell clusters to use when
#' calculating gene usage
#' @param chain Chain to use for calculating gene usage. Set to NULL to include
#' all chains.
#' @param chain_col meta.data column containing chains for each cell
#' @param sep Separator used for storing per cell V(D)J data
#' @return data.frame containing gene usage summary
#'
#' @examples
#' # Calculate V(D)J segment usage for all cells
#' calc_vdj_usage(
#'   vdj_so,
#'   gene_cols = "v_gene"
#' )
#'
#' # Calculate gene usage separately for cell clusters
#' calc_vdj_usage(
#'   vdj_so,
#'   gene_cols = "v_gene",
#'   cluster_col = "orig.ident"
#' )
#'
#' # Calculate gene usage for a specific chain
#' calc_vdj_usage(
#'   vdj_so,
#'   gene_cols = "v_gene",
#'   chain = c("IGH", "IGK")
#' )
#'
#' # Calculate paired usage of V(D)J segments
#' calc_vdj_usage(
#'   vdj_so,
#'   gene_cols = c("v_gene", "j_gene"),
#' )
#'
#' # Using calc_vdj_usage outside of Seurat
#' # SingleCellExperiment objects or data.frames containing V(D)J data are also
#' # compatible. If a data.frame is provided, cell barcodes should be stored as
#' # row names.
#' calc_vdj_usage(vdj_sce, gene_cols = "v_gene")
#'
#' df <- vdj_so@meta.data
#'
#' calc_vdj_usage(df, gene_cols = "v_gene")
#'
#' @export
calc_vdj_usage <- function(input, gene_cols, cluster_col = NULL, chain = NULL,
                           chain_col = "chains", sep = ";") {

  # Format input data
  sep_cols <- gene_cols

  if (!is.null(chain)) {
    sep_cols <- c(sep_cols, chain_col)
  }

  vdj_cols <- c(".cell_id", cluster_col, sep_cols)

  meta <- .get_meta(input)
  meta <- dplyr::select(meta, all_of(vdj_cols))

  meta <- dplyr::filter(meta, across(
    all_of(gene_cols),
    ~ !is.na(.x)
  ))

  res <- dplyr::mutate(meta, across(
    all_of(sep_cols),
    ~ strsplit(as.character(.x), sep)
  ))

  # Filter chains
  if (!is.null(chain)) {
    res <- .filter_chains(
      res,
      vdj_cols  = gene_cols,
      chain     = chain,
      chain_col = chain_col,
      col_names = "{.col}",
      empty_val = "None"
    )

    res <- dplyr::select(res, -all_of(chain_col))
  }

  res <- tidyr::unnest(res, cols = all_of(gene_cols))
  res <- dplyr::distinct(res)

  # Count genes used
  res <- dplyr::group_by(res, !!!syms(gene_cols))

  if (!is.null(cluster_col)) {
    res <- dplyr::group_by(res, !!sym(cluster_col), .add = TRUE)
  }

  res <- dplyr::summarize(
    res,
    n_cells = dplyr::n_distinct(meta$.cell_id),
    freq    = dplyr::n_distinct(.data$.cell_id),
    .groups = "drop"
  )

  # Calculate percentage used
  if (!is.null(cluster_col)) {
    clsts       <- pull(meta, cluster_col)
    clst_counts <- table(clsts)
    clsts       <- unique(clsts)

    res <- tidyr::pivot_wider(
      res,
      names_from  = all_of(cluster_col),
      values_from = .data$freq,
      values_fill = 0
    )

    res <- tidyr::pivot_longer(
      res,
      cols      = all_of(clsts),
      names_to  = cluster_col,
      values_to = "freq"
    )

    res <- dplyr::mutate(
      res,
      n_cells = as.numeric(clst_counts[!!sym(cluster_col)])
    )

    res <- dplyr::relocate(res, .data$n_cells, .before = .data$freq)
  }

  res <- dplyr::mutate(res, pct = (.data$freq / .data$n_cells) * 100)

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
#' @param type Type of plot to create, can be either 'heatmap' or 'bar'. If
#' multiple columns are provided to gene_cols, a heatmap is generated.
#' @param plot_colors Character vector containing colors to use for plot. If a
#' bar graph is created this will specify how to color cell clusters. For a
#' heatmap, these colors will be used to generate the color gradient.
#' @param vdj_genes V(D)J genes to plot, if NULL the top genes will be shown
#' @param n_genes Number of top genes to plot based on usage. If cluster_col is
#' provided, top genes will be identified for each cluster.
#' @param plot_lvls Levels to use for ordering clusters
#' @param yaxis Units to plot on the y-axis, either 'frequency' or 'percent'
#' @param chain_col meta.data column containing chains for each cell
#' @param sep Separator used for storing per-chain V(D)J data for each cell
#' @param ... Additional arguments to pass to ggplot2, e.g. color, fill, size,
#' linetype, etc.
#' @return ggplot object
#'
#' @examples
#' # Plot V(D)J segment usage for all cells
#' plot_vdj_usage(
#'   vdj_so,
#'   gene_cols = "v_gene"
#' )
#'
#' # Plot gene usage separately for cell clusters
#' plot_vdj_usage(
#'   vdj_so,
#'   gene_cols = "v_gene",
#'   cluster_col = "orig.ident"
#' )
#'
#' # Plot gene usage for a specific chain
#' plot_vdj_usage(
#'   vdj_so,
#'   gene_cols = "v_gene",
#'   chain = c("IGH", "IGK")
#' )
#'
#' # Plot gene usage for a specific chain
#' plot_vdj_usage(
#'   vdj_so,
#'   gene_cols = "v_gene",
#'   chain = c("IGH", "IGK")
#' )
#'
#' # Create a heatmap
#' plot_vdj_usage(
#'   vdj_so,
#'   gene_cols = "v_gene",
#'   type = "heatmap"
#' )
#'
#' # Plot paired usage of V(D)J segments
#' plot_vdj_usage(
#'   vdj_so,
#'   gene_cols = c("v_gene", "j_gene"),
#' )
#'
#' # Specify colors to use for each cell cluster
#' plot_vdj_usage(
#'   vdj_so,
#'   gene_cols = "v_gene",
#'   cluster_col = "orig.ident",
#'   plot_colors = c(avid_2 = "blue", avid_1 = "green")
#' )
#'
#' # Specify order to use for plotting cell clusters
#' plot_vdj_usage(
#'   vdj_so,
#'   gene_cols = "v_gene",
#'   cluster_col = "orig.ident",
#'   plot_lvls = c("avid_2", "avid_1")
#' )
#'
#' # Specify certain V(D)J genes to include in plot
#' plot_vdj_usage(
#'   vdj_so,
#'   gene_cols = "v_gene",
#'   vdj_genes = c("IGKV5-43", "IGLV1", "IGHV1-64", "IGHV1-11")
#' )
#'
#' # Specify the number of top V(D)J genes to include in plot
#' plot_vdj_usage(
#'   vdj_so,
#'   gene_cols = "v_gene",
#'   n_genes = 10
#' )
#'
#' # Plot the frequency of each V(D)J segment instead of percent
#' plot_vdj_usage(
#'   vdj_so,
#'   gene_cols = "v_gene",
#'   yaxis = "frequency"
#' )
#'
#' # Using plot_vdj_usage outside of Seurat
#' # SingleCellExperiment objects or data.frames containing V(D)J data are also
#' # compatible. If a data.frame is provided, cell barcodes should be stored as
#' # row names.
#' plot_vdj_usage(vdj_sce, gene_cols = "v_gene")
#'
#' df <- vdj_so@meta.data
#'
#' plot_vdj_usage(df, gene_cols = "v_gene")
#'
#' @export
plot_vdj_usage <- function(input, gene_cols, cluster_col = NULL, chain = NULL, type = "bar",
                           plot_colors = NULL, vdj_genes = NULL, n_genes = 50, plot_lvls = NULL,
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

  if (identical(yaxis, "frequency")) {
    usage_col <- "freq"
  }

  # Calculate gene usage
  plt_dat <- calc_vdj_usage(
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

  # If vdj_genes provided check plt_dat
  top_genes <- unique(vdj_genes)

  if (!is.null(top_genes)) {
    absent <- unlist(plt_dat[, gene_cols], use.names = FALSE)
    absent <- top_genes[!top_genes %in% absent]

    if (identical(absent, top_genes)) {
      stop("None of the provided genes were found")

    } else if (!purrr::is_empty(absent)) {
      absent <- paste0(absent, collapse = ", ")

      warning("Some genes not found: ", absent)
    }

    # If vdj_genes not provided, identify top n_genes for each cluster based on
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
    if (identical(type, "bar")) {
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

  if (!is.null(names(res))) {
    res <- purrr::imap(res, ~ .x + ggplot2::labs(title = .y))
  }

  if (length(res) == 1) {
    return(res[[1]])
  }

  res
}

