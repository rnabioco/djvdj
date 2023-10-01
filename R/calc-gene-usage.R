#' Calculate V(D)J segment usage
#'
#' Quantify the usage of different V(D)J segments for each cell cluster. The
#' usage of two V(D)J segments can also be calculated for a single chain.
#' For example, calc_gene_usage() can calculate the frequency that different
#' heavy chain V and J segments appear together.
#'
#' @param input Object containing V(D)J data. If a data.frame is provided, the
#' cell barcodes should be stored as row names.
#' @param data_cols meta.data column(s) containing V(D)J genes identified for
#' each clonotype. If multiple columns are provided, paired usage of genes will
#' be calculated.
#' @param cluster_col meta.data column containing cell clusters to use when
#' calculating gene usage
#' @param chain Chain(s) to use for calculating gene usage. Set to `NULL` to
#' include all chains.
#' @param chain_col meta.data column containing chains for each cell
#' @param prefix Prefix to add to new columns
#' @param return_df Return results as a data.frame. If `FALSE`, results will be
#' added to the input object.
#' @param sep Separator used for storing per cell V(D)J data
#' @return data.frame containing gene usage summary
#' @seealso [plot_gene_usage()], [calc_gene_pairs()], [plot_gene_pairs()]
#'
#' @examples
#' # Calculate V(D)J segment usage for all cells
#' calc_gene_usage(
#'   vdj_sce,
#'   data_cols = "v_gene"
#' )
#'
#' # Calculate gene usage separately for cell clusters
#' calc_gene_usage(
#'   vdj_sce,
#'   data_cols = "v_gene",
#'   cluster_col = "orig.ident"
#' )
#'
#' # Calculate gene usage for a specific chain(s)
#' calc_gene_usage(
#'   vdj_sce,
#'   data_cols = "v_gene",
#'   chain = c("IGK", "IGL")
#' )
#'
#' # Calculate paired usage of V(D)J segments
#' calc_gene_usage(
#'   vdj_sce,
#'   data_cols = c("v_gene", "j_gene"),
#' )
#'
#' @export
calc_gene_usage <- function(input, data_cols, cluster_col = NULL, chain = NULL,
                            chain_col = global$chain_col,
                            prefix = paste0(data_cols[1], "_"),
                            return_df = FALSE, sep = global$sep) {

  # Check that columns are present in object
  .check_obj_cols(
    input, data_cols, cluster_col, chain = chain, chain_col = chain_col
  )

  # Check input classes
  .check_args(cluster_col = list(allow_null = TRUE, len_one = FALSE))

  # Format input data
  sep_cols <- data_cols

  if (!is.null(chain)) sep_cols <- c(sep_cols, chain_col)

  vdj_cols <- c(global$cell_col, cluster_col, sep_cols)

  meta <- .get_meta(input)
  vdj  <- dplyr::select(meta, all_of(vdj_cols))

  # Calculate frequency
  if (length(data_cols) > 1) return_df <- TRUE
  if (return_df)             prefix <- ""

  vdj <- .calc_freq(
    df_in       = vdj,
    data_cols   = data_cols,
    cluster_col = cluster_col,
    chain       = chain,
    chain_col   = chain_col,
    per_chain   = TRUE,
    return_df   = return_df,
    prefix      = prefix,
    sep         = sep
  )

  if (return_df) return(vdj)

  # Format results
  res <- dplyr::left_join(meta, vdj, by = global$cell_col)
  res <- .add_meta(input, meta = res)

  res
}

#' Calculate paired usage of V(D)J segments across chains
#'
#' Quantify the paired usage of V(D)J segments across two chains. For example,
#' calc_gene_pairs() can calculate the frequency that different TRA and TRB V
#' segments appear together.
#'
#' @param input Object containing V(D)J data. If a data.frame is provided, the
#' cell barcodes should be stored as row names.
#' @param data_col meta.data column containing V(D)J genes identified for each
#' clonotype.
#' @param chains Chains to use for calculating usage of different gene pairs.
#' This should be a character vector containing the two chains to use for
#' calculations, e.g. `c("IGH", "IGK")`.
#' @param cluster_col meta.data column containing cell clusters to use when
#' calculating gene usage
#' @param chain_col meta.data column containing chains for each cell
#' @param sep Separator used for storing per cell V(D)J data
#' @return data.frame containing gene pair summary
#' @seealso [plot_gene_pairs()], [calc_gene_usage()], [plot_gene_usage()]
#'
#' @examples
#' # Calculate the frequency of different V genes for IGH and IGK chains
#' calc_gene_pairs(
#'   vdj_sce,
#'   data_col    = "v_gene",
#'   chains      = c("IGH", "IGK"),
#'   cluster_col = "orig.ident"
#' )
#'
#' @export
calc_gene_pairs <- function(input, data_col, chains, cluster_col = NULL,
                            chain_col = global$chain_col, sep = global$sep) {

  # Check that columns are present in object
  .check_obj_cols(
    input, data_col, cluster_col, chain = chains, chain_col = chain_col
  )

  # Check input classes
  .check_args(cluster_col = list(allow_null = TRUE, len_one = FALSE))

  # Format input data
  vdj_cols <- c(global$cell_col, data_col, cluster_col, chain_col)

  res <- fetch_vdj(
    input,
    data_cols    = c(data_col, chain_col),
    filter_cells = TRUE,
    unnest       = FALSE,
    per_chain    = TRUE,
    sep          = sep
  )

  res  <- dplyr::select(res, all_of(vdj_cols))

  res <- .filter_chains(
    res,
    data_cols  = c(data_col, chain_col),
    chain      = chains,
    chain_col  = chain_col,
    col_names  = "{.col}",
    allow_dups = FALSE
  )

  res <- dplyr::filter(res, !is.na(!!sym(data_col)))

  res <- tidyr::pivot_wider(
    res,
    names_from  = !!sym(chain_col),
    values_from = !!sym(data_col),
    values_fill = "None"
  )

  # res <- dplyr::filter(res, dplyr::if_all(dplyr::all_of(chains), ~ !is.na(.x)))

  # Calculate frequency
  res <- .calc_freq(
    df_in       = res,
    data_cols   = chains,
    cluster_col = cluster_col,
    per_chain   = FALSE,
    return_df   = TRUE,
    sep         = sep
  )

  res
}


#' Plot V(D)J segment usage
#'
#' Plot the usage of different V(D)J segments for each cell cluster. The
#' usage of two V(D)J segments can also be plotted for a single chain.
#'
#' @param input Object containing V(D)J data. If a data.frame is provided, the
#' cell barcodes should be stored as row names.
#' @param data_cols meta.data column containing genes for each clonotype,
#' provide a vector with two column names to plot paired usage of genes
#' @param cluster_col meta.data column containing cell clusters to use for
#' calculating gene usage
#' @param group_col meta.data column to use for grouping cluster IDs present in
#' cluster_col. This is useful when there are multiple replicates or patients
#' for each treatment condition.
#' @param genes An integer specifying the number of genes to plot, or
#' a vector giving the names of genes to include.
#' @param method Method to use for plotting, possible values are:
#'
#' - 'bar', create a bargraph, this is the default when a single column is
#' passed to the data_cols argument
#' - 'boxplot', create boxplots, this can only be used when group_col is
#' provided
#' - 'heatmap', create a heatmap, this is the default when two columns are
#' passed to the data_cols argument
#' - 'circos', create a circos plot, this requires two columns to be provided
#' to the data_cols argument
#'
#' @param units Units to plot on the y-axis, either 'frequency' or 'percent'
#' @param return_list Should a list of plots be returned, if FALSE plots will be
#' combined and arranged into panels
#'
#'
#' @param plot_colors Character vector containing colors to use for plot. If a
#' bar graph is created this will specify how to color cell clusters. For a
#' heatmap, these colors will be used to generate the color gradient.
#' @param plot_lvls Levels to use for ordering clusters
#' @param trans Transformation to use when plotting segment usage, e.g.
#' 'log10'. By default values are not transformed, refer to
#' [ggplot2::continuous_scale()] for more options.
#' @param rotate_labels Should labels on circos plot be rotated to reduce
#' overlapping text
#' @param panel_nrow The number of rows to use for arranging plots when
#' return_list is FALSE
#' @param show_points If `TRUE` data points will be shown on boxplots, the point
#' size can be adjusted using the `point.size` parameter
#' @param show_zeros If `TRUE` cell labels that are missing from a cluster will
#' still be shown on the plot
#' @param n_label Location on plot where n label should be added, this is only
#' applicable when `method` is 'bar' and can be any combination of the
#' following:
#'
#' - 'corner', display the total number of cells plotted in the top right
#'   corner, the position of the label can be modified by passing `x` and `y`
#'   specifications with the `label_params` argument
#' - 'legend', display the number of cells plotted for each group shown in the
#'   plot legend
#' - 'none', do not display the number of cells plotted
#'
#' @param p_label Specification indicating how p-values should be labeled on
#' plot, this can one of the following:
#'
#' - 'none', do not display p-values
#' - 'all', show p-values for all groups
#' - A named vector providing p-value cutoffs and labels to display,
#'   e.g. `c('*' = 0.05, '**' = 0.01, '***' = 0.001)`. The keyword 'value' can
#'   be used to display the p-value for those less than a certain cutoff,
#'   e.g. `c(value = 0.05, ns = Inf)` will show significant p-values, all others
#'   will be labeled 'ns'.
#'
#' @param p_method Method to use for calculating p-values. By default when
#' comparing two groups a t-test will be performed, when comparing more than
#' two groups the Kruskal-Wallis test will be used. p-values are adjusted for
#' multiple testing using Bonferroni correction. Possible methods include:
#'
#' - 't', two sample t-test performed with `stats::t.test()`
#' - 'wilcox', Wilcoxon rank sum test performed with `stats::wilcox.test()`
#' - 'kruskal', Kruskal-Wallis test performed with `stats::kruskal.test()`
#'
#' @param p_file File path to save table containing p-values for each
#' comparison.
#' @param label_params Named list providing additional parameters to modify
#' n label aesthetics, e.g. list(size = 4, color = "red")
#' @param ... Additional arguments to pass to plotting function,
#' [ggplot2::geom_col()] for bargraph, [ggplot2::geom_tile()] for heatmap,
#' [circlize::chordDiagram()] for circos plot
#' @param chain Chain to use for calculating gene usage, set to NULL to include
#' all chains
#' @param chain_col meta.data column containing chains for each cell
#' @param sep Separator used for storing per-chain V(D)J data for each cell
#' @seealso [calc_gene_usage()], [calc_gene_pairs()], [plot_gene_pairs()]
#' @return ggplot object
#'
#' @examples
#' # Plot V(D)J segment usage for all cells
#' plot_gene_usage(
#'   vdj_sce,
#'   data_cols = "v_gene"
#' )
#'
#' # Plot gene usage separately for cell clusters
#' plot_gene_usage(
#'   vdj_sce,
#'   data_cols = "v_gene",
#'   cluster_col = "orig.ident"
#' )
#'
#' # Plot gene usage for a specific chain
#' plot_gene_usage(
#'   vdj_sce,
#'   data_cols = "v_gene",
#'   chain = c("IGH", "IGK")
#' )
#'
#' # Plot paired usage of V(D)J segments
#' plot_gene_usage(
#'   vdj_sce,
#'   data_cols = c("v_gene", "j_gene"),
#'   type = "circos"
#' )
#'
#' # Specify colors to use for each cell cluster
#' plot_gene_usage(
#'   vdj_sce,
#'   data_cols = "v_gene",
#'   cluster_col = "orig.ident",
#'   plot_colors = c(avid_2 = "blue", avid_1 = "green")
#' )
#'
#' # Specify order to use for plotting cell clusters
#' plot_gene_usage(
#'   vdj_sce,
#'   data_cols = "v_gene",
#'   cluster_col = "orig.ident",
#'   plot_lvls = c("avid_2", "avid_1")
#' )
#'
#' # Specify certain V(D)J genes to include in plot
#' plot_gene_usage(
#'   vdj_sce,
#'   data_cols = "v_gene",
#'   vdj_genes = c("IGKV5-43", "IGLV1", "IGHV1-64")
#' )
#'
#' # Specify the number of top V(D)J genes to include in plot
#' plot_gene_usage(
#'   vdj_sce,
#'   data_cols = "v_gene",
#'   genes = 10
#' )
#'
#' @export
plot_gene_usage <- function(input, data_cols, cluster_col = NULL,
                            group_col = NULL, method = NULL,
                            units = "percent", genes = 20,
                            return_list = FALSE,
                            plot_colors = NULL, plot_lvls = NULL,
                            trans = "identity", rotate_labels = FALSE,
                            panel_nrow = NULL, show_points = TRUE,
                            show_zeros = TRUE,
                            n_label = NULL, p_label = c(value = 0.05),
                            p_method = NULL, p_file = NULL,
                            label_params = list(),
                            ...,
                            chain = NULL,
                            chain_col = global$chain_col,
                            sep = global$sep) {

  # Check that columns are present in object
  .check_obj_cols(
    input,
    data_cols, cluster_col, group_col, chain = chain, chain_col = chain_col
  )

  # Check input classes
  .check_args(method = list(allow_null = TRUE))

  # Check input values
  paired <- length(data_cols) == 2
  method <- method %||% ifelse(paired, "heatmap", "bar")

  if (identical(method, "circos")) units  <- "frequency"

  .check_usage_args(method, data_cols, group_col, units, paired)
  .check_group_cols(cluster_col, group_col, input)

  # Set y-axis
  usage_col <- switch(units, frequency = "freq", percent   = "pct")

  # Calculate gene usage
  plt_dat <- calc_gene_usage(
    input       = input,
    data_cols   = data_cols,
    cluster_col = c(cluster_col, group_col),
    chain       = chain,
    chain_col   = chain_col,
    sep         = sep,
    return_df   = TRUE
  )

  plt_dat <- dplyr::filter(plt_dat, dplyr::if_all(
    dplyr::all_of(data_cols), ~ .x != "None"
  ))

  # If vdj_genes provided check plt_dat
  if (is.character(genes)) {
    absent <- unlist(plt_dat[, data_cols], use.names = FALSE)
    absent <- genes[!genes %in% absent]

    if (identical(absent, genes)) {
      cli::cli_abort("None of the provided genes were found")

    } else if (!purrr::is_empty(absent)) {
      cli::cli_warn("The following genes were not found: {absent}")
    }

    # if vdj_genes are provided and two data_cols are provided, gene pairs
    # containing at least one of the specified genes will be included
    plt_dat <- dplyr::filter(plt_dat, dplyr::if_any(
      dplyr::all_of(data_cols), ~ .x %in% genes
    ))

  # If vdj_genes not provided, identify top genes for each cluster based on
  # usage
  } else {
    plt_dat <- .filter_top_genes(
      plt_dat,
      dat_col  = usage_col,
      gn_col   = data_cols,
      clst_col = cluster_col,
      n        = genes
    )
  }

  # Plotting arguments
  gg_args <- list(
    df_in    = plt_dat,
    gn_col   = data_cols,
    dat_col  = usage_col,
    clst_col = cluster_col,
    method   = method,
    clrs     = plot_colors,
    lvls     = plot_lvls,
    trans    = trans,
    ttl      = .get_axis_label(units),
    n_row    = panel_nrow,
    ...
  )

  # Create plot for paired usage
  # heatmap and circos return invisibly
  if (paired) {
    gg_args$rotate_labels <- rotate_labels
    gg_args$return_list   <- return_list

    res <- .lift(.plot_paired_usage)(gg_args)

    return(invisible())
  }

  # Create plot for single usage
  gg_args$clst_col     <- cluster_col
  gg_args$grp_col      <- group_col
  gg_args$show_points  <- show_points
  gg_args$show_zeros   <- show_zeros
  gg_args$n_label      <- n_label
  gg_args$p_label      <- p_label
  gg_args$p_method     <- p_method
  gg_args$p_file       <- p_file
  gg_args$label_params <- label_params

  res <- .lift(.plot_single_usage)(gg_args)

  res
}

#' Plot paired usage of V(D)J segments across chains
#'
#' Plot the paired usage of V(D)J segments across two chains. For example,
#' plot_gene_pairs() can be used to plot the frequency that different TRA and
#' TRB V segments appear together.
#'
#' @param input Object containing V(D)J data. If a data.frame is provided, the
#' cell barcodes should be stored as row names.
#' @param data_col meta.data column containing V(D)J genes identified for each
#' clonotype.
#' @param chains Chains to use for calculating usage of different gene pairs.
#' This should be a character vector containing the two chains to use for
#' calculations, e.g. `c("IGH", "IGK")`.
#' @param cluster_col meta.data column containing cell clusters to use when
#' calculating gene usage
#' @param genes An integer specifying the number of genes to plot, or
#' a vector giving the names of genes to include.
#' @param chain_col meta.data column containing chains for each cell
#' @param method Method to use for plotting, possible values are:
#'
#' - 'heatmap', create a heatmap, this is the default when two columns are
#' passed to the data_cols argument
#' - 'circos', create a circos plot, this requires two columns to be provided
#' to the data_cols argument
#'
#' @param units Units to show on scale, either 'frequency' or 'percent'
#' @param return_list Should a list of plots be returned, if FALSE plots will be
#' combined and arranged into panels
#' @param sep Separator used for storing per-chain V(D)J data for each cell
#'
#'
#' @param plot_colors Character vector containing colors to use for plot. If a
#' bar graph is created this will specify how to color cell clusters. For a
#' heatmap, these colors will be used to generate the color gradient.
#' @param plot_lvls Levels to use for ordering clusters
#' @param trans Transformation to use when plotting segment usage, e.g.
#' 'log10'. By default values are not transformed, refer to
#' [ggplot2::continuous_scale()] for more options.
#' @param rotate_labels Should labels on circos plot be rotated to reduce
#' overlapping text
#' @param panel_nrow The number of rows to use for arranging plots when
#' return_list is FALSE
#' @param ... Additional arguments to pass to plotting function,
#' [ggplot2::geom_tile()] for heatmap, [circlize::chordDiagram()] for circos
#' plot
#' @seealso [calc_gene_pairs()], [calc_gene_usage()], [plot_gene_usage()]
#' @return ggplot object
#'
#' @examples
#' # Plot the frequency of different V genes for IGH and IGK chains
#' plot_gene_pairs(
#'   vdj_sce,
#'   data_col    = "v_gene",
#'   chains      = c("IGH", "IGK"),
#'   cluster_col = "orig.ident"
#' )
#'
#' @export
plot_gene_pairs <- function(input, data_col, chains, cluster_col = NULL,
                            genes = 20,
                            chain_col = global$chain_col,
                            method = "heatmap",
                            units = "percent",
                            return_list = FALSE,
                            sep = global$sep,
                            plot_colors = NULL, plot_lvls = NULL,
                            trans = "identity",
                            rotate_labels = FALSE,
                            panel_nrow = NULL, ...) {

  # Check that columns are present in object
  .check_obj_cols(
    input,
    data_col, cluster_col, chain = chains, chain_col = chain_col
  )

  # Check input classes
  .check_args(method = list(allow_null = TRUE))

  # Check input values
  if (identical(method, "circos")) units  <- "frequency"

  .check_usage_args(method, data_col, grp_col = NULL, units, paired = TRUE)

  # Set y-axis
  usage_col <- switch(units, frequency = "freq", percent   = "pct")

  # Calculate gene usage
  plt_dat <- calc_gene_pairs(
    input       = input,
    data_col    = data_col,
    chains      = chains,
    cluster_col = cluster_col,
    chain_col   = chain_col,
    sep         = sep
  )

  # If vdj_genes provided check plt_dat
  if (is.character(genes)) {
    absent <- unlist(plt_dat[, chains], use.names = FALSE)
    absent <- genes[!genes %in% absent]

    if (identical(absent, genes)) {
      cli::cli_abort("None of the provided genes were found")

    } else if (!purrr::is_empty(absent)) {
      cli::cli_warn("The following genes were not found: {absent}")
    }

    # if vdj_genes are provided and two data_cols are provided, gene pairs
    # containing at least one of the specified genes will be included
    plt_dat <- dplyr::filter(plt_dat, dplyr::if_any(
      dplyr::all_of(chains), ~ .x %in% genes
    ))

  # If vdj_genes not provided, identify top genes for each cluster based on
  # usage
  } else {
    plt_dat <- .filter_top_genes(
      plt_dat,
      dat_col  = usage_col,
      gn_col   = chains,
      clst_col = cluster_col,
      n        = genes
    )
  }

  # Plotting arguments
  gg_args <- list(
    df_in    = plt_dat,
    gn_col   = chains,
    dat_col  = usage_col,
    clst_col = cluster_col,
    method   = method,
    clrs     = plot_colors,
    lvls     = plot_lvls,
    trans    = trans,
    ttl      = .get_axis_label(units),
    n_row    = panel_nrow,
    ...
  )

  # Create plot for paired usage
  # heatmap and circos return invisibly
  gg_args$rotate_labels <- rotate_labels
  gg_args$return_list   <- return_list

  .lift(.plot_paired_usage)(gg_args)
}

#' Plot usage for single gene column
#'
#' @param df_in data.frame
#' @param gn_col Column containing genes
#' @param dat_col Column containing gene usage metric
#' @param clst_col Column containing clusters
#' @param grp_col Column with IDs to use for grouping clusters
#' @param method Method to use for generating plot
#' @param clrs Plot colors
#' @param lvls Levels for ordering clusters
#' @param show_points If `TRUE` data points will be shown on boxplots, the point
#' size can be adjusted using the `point.size` parameter
#' @param show_zeros If `TRUE` cell labels that are missing from a cluster will
#' still be shown on the plot
#' @param n_label n label specification
#' @param p_label p label specification
#' @param p_method Method to calculate p-values
#' @param p_file File to save p-values
#' @param label_params List of parameters to modify label aesthetics
#' @param trans Method to use for transforming data
#' @param ttl Title for y-axis or legend depending on type of graph
#' @param order Should genes be ordered based on usage
#' @param n_row Number of rows to use for arranging plots
#' @param ... Additional parameters to pass to plotting function
#' @noRd
.plot_single_usage <- function(df_in, gn_col, dat_col, method, clst_col = NULL,
                               grp_col = NULL, lvls = NULL, show_points = TRUE,
                               show_zeros = TRUE, n_label = NULL,
                               p_label = c(value = 0.05), p_method = NULL,
                               p_file = NULL, label_params = list(), n_row = 1,
                               ttl = dat_col, order = TRUE, ...) {

  # Order clusters based on plot_lvls
  lvl_col <- grp_col %||% clst_col
  df_in   <- .set_lvls(df_in, lvl_col, lvls)

  # Order plot levels
  if (order) {
    decr <- (!is.null(grp_col) && !identical(method, "heatmap")) ||
      identical(method, "bar")

    g_col <- sym(gn_col)

    df_in <- dplyr::mutate(
      df_in,
      !!g_col := .order_lvls(!!g_col, !!sym(dat_col), max, decreasing = decr)
    )
  }

  # Set possible n label parameters
  psbl_labs <- switch(
    method,
    heatmap = c("none", "axis"),
    bar     = c("none", "corner")
  )

  if (identical(method, "bar") && !is.null(clst_col)) {
    psbl_labs <- c(psbl_labs, "legend")
  }

  n_label <- n_label %||% psbl_labs
  n_label <- n_label[n_label %in% psbl_labs]

  # Set n label data
  n_lab_dat <- list()

  grp_cols  <- c(.n = "n_cells", clst_col, grp_col)

  n_lab_dat$corner <- dplyr::distinct(df_in, !!!syms(grp_cols))

  n_lab_dat$legend <- n_lab_dat$axis <- n_lab_dat$corner

  if (
    (!is.null(grp_col) || !is.null(clst_col)) &&
    (identical(method, "bar") || identical(method, "boxplot"))
  ) {
    n_lab_dat$axis <- dplyr::rename(df_in, .n = "freq")
  }

  # Set common arguments
  gg_args <- list(
    df_in        = df_in,
    n_label      = n_label,
    label_data   = n_lab_dat,
    label_params = label_params,
    n_fn         = sum,
    ...
  )

  if (identical(method, "bar") || !is.null(grp_col)) {
    gg_args$x <- gn_col
    gg_args$y <- dat_col
  }

  # Create grouped boxplot
  if (!is.null(grp_col) && !identical(method, "heatmap")) {
    .check_possible_values(p_method = c("t", "wilcox", "kruskal"))

    gg_args$grp <- gg_args$.color <- gg_args$.fill <- grp_col

    gg_args$clst        <- clst_col
    gg_args$method      <- method
    gg_args$y_ttl       <- ttl
    gg_args$show_points <- show_points
    gg_args$show_zeros  <- show_zeros
    gg_args$p_label     <- p_label
    gg_args$p_method    <- p_method
    gg_args$p_file      <- p_file

    res <- .lift(.create_grouped_plot)(gg_args)

    return(res)

  # Create bargraph
  } else if (identical(method, "bar")) {

    # Remove zeros from data
    if (!show_zeros) {
      df_in <- dplyr::filter(df_in, !!sym(dat_col) > 0)
    }

    gg_args$df_in <- df_in
    gg_args$.fill <- gg_args$.color <- clst_col
    gg_args$y_ttl <- ttl

    res <- .lift(.create_bars)(gg_args)

    return(res)
  }

  # Create heatmap
  # heatmap can be generated for both grouped and ungrouped plots
  gg_args$x       <- clst_col
  gg_args$y       <- gn_col
  gg_args$grp     <- grp_col
  gg_args$.fill   <- dat_col
  gg_args$lgd_ttl <- ttl
  gg_args$nrow    <- n_row
  gg_args$scales  <- "free_x"

  res <- .lift(.create_gg_heatmap)(gg_args)

  res
}

#' Plot paired usage for two gene columns
#'
#' @param df_in data.frame
#' @param gn_col Columns containing genes
#' @param dat_col Column containing gene usage metric
#' @param clst_col Column containing clusters
#' @param method Method to use for generating plot
#' @param clrs Plot colors
#' @param lvls Levels for ordering clusters
#' @param trans Method to use for transforming data
#' @param ttl Legend title
#' @param order Should genes be ordered based on usage
#' @param return_list Should a list of plot be returned, if FALSE plots will be
#' combined
#' @param n_row Number of rows to use for arranging plots
#' @param ... Additional parameters to pass to plotting function
#' @noRd
.plot_paired_usage <- function(df_in, gn_col, dat_col, clst_col = NULL,
                               method, clrs = NULL, lvls = NULL,
                               trans = "identity", ttl = dat_col, order = TRUE,
                               rotate_labels = FALSE, return_list = FALSE,
                               n_row = NULL, ...) {

  # Split input data into clusters and order based on lvls
  res <- list(df_in)

  if (!is.null(clst_col)) {
    res <- .set_lvls(df_in, clst_col, lvls)
    res <- split(res, df_in[[clst_col]])

    if (!is.null(lvls)) res <- res[unique(lvls)]
  }

  # Format data for plotting
  # sort level order, ascending order for y
  res <- purrr::map(res, ~ {
    .x <- dplyr::select(.x, dplyr::all_of(c(gn_col, dat_col)))
    d  <- tidyr::expand(.x, !!!syms(gn_col))
    d  <- dplyr::left_join(d, .x, by = gn_col)

    d <- dplyr::mutate(
      d,
      !!sym(dat_col) := tidyr::replace_na(!!sym(dat_col), 0)
    )

    if (order) {
      d <- dplyr::mutate(d, dplyr::across(all_of(gn_col), ~ {
        decr <- dplyr::cur_column() == gn_col[1]

        .order_lvls(.x, !!sym(dat_col), max, decreasing = decr)
      }))
    }

    d
  })

  plt_args <- list(clrs = clrs, ...)
  add_ttl  <- !is.null(names(res))

  # Set number of rows/columns for plots
  nplts <- length(res)
  sq    <- sqrt(nplts)
  nr    <- floor(sq)

  n_row <- n_row %||% nr
  n_col <- ceiling(nplts / n_row)

  # Create heatmap
  if (identical(method, "heatmap")) {
    plt_args$x       <- gn_col[1]
    plt_args$y       <- gn_col[2]
    plt_args$.fill   <- dat_col
    plt_args$lgd_ttl <- ttl
    plt_args$trans   <- trans

    res <- purrr::imap(res, ~ {
      plt_args$df_in <- .x

      if (add_ttl) plt_args$plt_ttl <- .y

      .lift(.create_gg_heatmap)(plt_args)
    })

    if (length(res) == 1) return(res[[1]])
    if (return_list)      return(res)

    # Arrange multiple plots
    wd <- 1 / n_col
    ht <- 1 / n_row

    graphics::plot.new()

    purrr::iwalk(unname(res), ~ {
      x <- (.y - 1) %% n_col + 1
      x <- x / n_col - (wd / 2)
      y <- n_row - (ceiling(.y / n_col))
      y <- y / n_row + (ht / 2)

      plt_vp <- grid::viewport(
        height = grid::unit(ht, "npc"),
        width  = grid::unit(wd, "npc"),
        x      = grid::unit(x, "npc"),
        y      = grid::unit(y, "npc")
      )

      print(.x, vp = plt_vp)
    })

    return(invisible())
  }

  # Create circos plot
  # use lift to pass dots into walk
  plt_args$symmetric <- FALSE
  plt_args$rotate_labels <- rotate_labels

  if (!return_list) graphics::par(mfrow = c(n_row, n_col))

  purrr::iwalk(res, ~ {
    d <- tidyr::pivot_wider(
      .x,
      names_from  = gn_col[1],
      values_from = dat_col
    )

    d <- tibble::column_to_rownames(d, gn_col[2])
    d <- as.matrix(d)

    plt_args$mat_in <- d

    if (add_ttl) plt_args$plt_ttl <- .y

    .lift(.create_circos)(plt_args)
  })

  return(invisible())
}

#' Filter genes based on top usage
#'
#' @param df_in data.frame
#' @param dat_col Column containing gene usage metric
#' @param gn_col Columns containing genes
#' @param clst_col Column containing clusters
#' @param n Number of top genes to pull
#' @noRd
.filter_top_genes <- function(df_in, dat_col, gn_col, clst_col, n) {

  gn_args <- list(
    df_in   = df_in,
    dat_col = dat_col,
    gn_col  = gn_col,
    n       = n
  )

  # Filter for single gene column
  if (length(gn_col) == 1) {
    top <- .lift(.get_top_genes)(gn_args)
    res <- dplyr::filter(df_in, !!sym(gn_col) %in% top)

    return(res)
  }

  # Filter for paired gene columns
  # will return n top genes for column 1, for column 2 gene to be included
  # must be a top gene or have the highest usage
  gn_args$clst_col <- clst_col

  top_genes <- purrr::map(gn_col, ~ {
    gn_args$gn_col <- .x
    .lift(.get_top_genes)(gn_args)
  })

  gn_1  <- sym(gn_col[1])
  gn_2  <- sym(gn_col[2])
  top_1 <- top_genes[[1]]
  top_2 <- top_genes[[2]]

  if (is.null(clst_col)) {
    res <- dplyr::ungroup(df_in)
    res <- dplyr::filter(res, !!gn_1 %in% top_1)

  } else {
    res <- dplyr::group_by(df_in, !!sym(clst_col))
    res <- dplyr::filter(
      res,
      !!gn_1 %in% top_1[[dplyr::cur_group()[[clst_col]]]]
    )
  }

  res <- dplyr::group_by(res, !!gn_1, .add = TRUE)
  res <- dplyr::mutate(res, rnk = row_number(desc(!!sym(dat_col))))

  if (is.null(clst_col)) {
    res <- dplyr::filter(res, !!gn_2 %in% top_2 | .data$rnk == 1)

  } else {
    res <- dplyr::filter(
      res,
      !!gn_2 %in% top_2[[dplyr::cur_group()[[clst_col]]]] | .data$rnk == 1
    )
  }

  res <- dplyr::ungroup(res)

  res
}

#' Get top genes based on usage
#'
#' @param df_in data.frame
#' @param dat_col Column containing gene usage metric
#' @param gn_col Columns containing genes
#' @param clst_col Column containing clusters
#' @param n Number of top genes to pull
#' @noRd
.get_top_genes <- function(df_in, dat_col, gn_col, clst_col = NULL, n) {

  grps <- c(gn_col, clst_col)
  res  <- dplyr::group_by(df_in, !!!syms(grps))

  res <- dplyr::summarize(
    res,
    !!sym(dat_col) := max(!!sym(dat_col)),
    .groups = "drop"
  )

  res <- dplyr::arrange(res, dplyr::desc(!!sym(dat_col)))

  if (!is.null(clst_col)) res <- dplyr::group_by(res, !!sym(clst_col))

  res <- dplyr::slice_max(
    res,
    order_by  = !!sym(dat_col),
    n         = n,
    with_ties = FALSE
  )

  res <- dplyr::ungroup(res)

  # Pull gene lists
  .pull_gns <- function(df_in, clmns) {
    res <- unlist(df_in[, clmns], use.names = FALSE)
    res <- unique(res)
    res
  }

  if (!is.null(clst_col)) {
    res <- split(res, res[[clst_col]])
    res <- purrr::map(res, .pull_gns, gn_col)

  } else {
    res <- .pull_gns(res, gn_col)
  }

  res
}

#' Check plot_gene_usage arguments
#'
#' @param typ type
#' @param gn_cols gene_cols
#' @param grp_col group_col
#' @param axis units
#' @param paired paired
#' @noRd
.check_usage_args <- function(method, gn_cols, grp_col, units, paired) {

  .check_possible_values(
    method = c("heatmap", "bar", "boxplot", "circos"),
    units  = c("percent", "frequency")
  )

  mets <- c("heatmap", "circos")

  if (paired && !method %in% mets) {
    cli::cli_abort(
      "`method` must be {.or {typs}} when two columns are passed to `data_cols`"
    )
  }

  if (identical(method, "circos") && !paired) {
    cli::cli_abort(
      "A circos plot can only be generated when two columns
       are passed to `data_cols`"
    )
  }

  if (length(gn_cols) > 2) {
    cli::cli_abort("Cannot specify more than two values for `data_cols`")
  }

  if (paired && !is.null(grp_col)) {
    cli::cli_abort(
      "`group_col` can only be used when a single column
       is passed to `data_cols`"
    )
  }
}

