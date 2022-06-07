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
#' calc_gene_usage(
#'   vdj_so,
#'   gene_cols = "v_gene"
#' )
#'
#' # Calculate gene usage separately for cell clusters
#' calc_gene_usage(
#'   vdj_sce,
#'   gene_cols = "v_gene",
#'   cluster_col = "orig.ident"
#' )
#'
#' # Calculate gene usage for a specific chain(s)
#' calc_gene_usage(
#'   vdj_so,
#'   gene_cols = "v_gene",
#'   chain = c("IGK", "IGL")
#' )
#'
#' # Calculate paired usage of V(D)J segments
#' calc_gene_usage(
#'   vdj_sce,
#'   gene_cols = c("v_gene", "j_gene"),
#' )
#'
#' @export
calc_gene_usage <- function(input, gene_cols, cluster_col = NULL, chain = NULL,
                            chain_col = "chains", sep = ";") {

  # Format input data
  sep_cols <- gene_cols

  if (!is.null(chain)) {
    sep_cols <- c(sep_cols, chain_col)
  }

  vdj_cols <- c(CELL_COL, cluster_col, sep_cols)



  # >>> SHOULD THIS BE REPLACED WITH FETCH_VDJ CALL???
  meta <- .get_meta(input)
  meta <- dplyr::select(meta, all_of(vdj_cols))

  meta <- dplyr::filter(meta, dplyr::if_all(
    all_of(gene_cols),
    ~ !is.na(.x)
  ))

  res <- dplyr::mutate(meta, across(
    all_of(sep_cols),
    ~ strsplit(as.character(.x), sep)
  ))
  # <<< SHOULD THIS BE REPLACED WITH FETCH_VDJ CALL???



  # Combine cluster_col columns for calculations
  if (!is.null(cluster_col)) {
    clst_nm <- paste0(cluster_col, collapse = "_")

    res <- tidyr::unite(
      res,
      !!sym(clst_nm),
      !!!syms(cluster_col),
      remove = FALSE
    )

    clst_key <- dplyr::select(res, -all_of(gene_cols))
  }

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
    res <- dplyr::group_by(res, !!!syms(c(clst_nm, cluster_col)), .add = TRUE)
  }

  res <- dplyr::summarize(
    res,
    n_cells = dplyr::n_distinct(meta[[CELL_COL]]),
    freq    = dplyr::n_distinct(!!sym(CELL_COL)),
    .groups = "drop"
  )

  # Calculate percentage used
  if (!is.null(cluster_col)) {
    clsts       <- pull(clst_key, clst_nm)
    clst_counts <- table(clsts)
    clsts       <- unique(clsts)

    res <- dplyr::select(res, -all_of(cluster_col), !!sym(clst_nm))

    res <- tidyr::pivot_wider(
      res,
      names_from  = clst_nm,
      values_from = .data$freq,
      values_fill = 0
    )

    res <- tidyr::pivot_longer(
      res,
      cols      = all_of(clsts),
      names_to  = clst_nm,
      values_to = "freq"
    )

    res <- dplyr::mutate(
      res,
      n_cells = as.numeric(clst_counts[!!sym(clst_nm)])
    )

    # Add group columns
    clst_key <- dplyr::distinct(clst_key, !!!syms(c(clst_nm, cluster_col)))

    res <- dplyr::left_join(res, clst_key, by = clst_nm)

    res <- dplyr::select(
      res,
      all_of(c(gene_cols, cluster_col, "n_cells", "freq"))
    )
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
#' @param group_col meta.data column to use for grouping cluster IDs present in
#' cluster_col. This is useful when there are multiple replicates or patients
#' for each treatment condition.
#' @param chain Chain to use for calculating gene usage, set to NULL to include
#' all chains
#' @param chain_col meta.data column containing chains for each cell
#' @param type Type of plot to create, possible values are:
#'
#' - 'bar', bargraphs can only be plotted when a single column is passed to the
#' gene_cols argument
#' - 'heatmap', heatmap implemented with ComplexHeatmap::Heatmap()
#' - 'circos', circos plot implemented with circlize::chordDiagram().
#' to create a circos plot, two columns must be passed to the gene_cols
#' argument
#'
#' @param plot_colors Character vector containing colors to use for plot. If a
#' bar graph is created this will specify how to color cell clusters. For a
#' heatmap, these colors will be used to generate the color gradient.
#' @param vdj_genes V(D)J genes to plot, if NULL the top genes will be shown
#' @param n_genes Number of top genes to plot based on usage. If cluster_col is
#' provided, top genes will be identified for each cluster.
#' @param plot_lvls Levels to use for ordering clusters
#' @param yaxis Units to plot on the y-axis, either 'frequency' or 'percent'
#' @param sep Separator used for storing per-chain V(D)J data for each cell
#' @param ... Additional arguments to pass to ggplot2, e.g. color, fill, size,
#' linetype, etc.
#' @return ggplot object
#'
#' @examples
#' # Plot V(D)J segment usage for all cells
#' plot_gene_usage(
#'   vdj_so,
#'   gene_cols = "v_gene"
#' )
#'
#' # Plot gene usage separately for cell clusters
#' plot_gene_usage(
#'   vdj_sce,
#'   gene_cols = "v_gene",
#'   cluster_col = "orig.ident"
#' )
#'
#' # Plot gene usage for a specific chain
#' plot_gene_usage(
#'   vdj_so,
#'   gene_cols = "v_gene",
#'   chain = c("IGH", "IGK")
#' )
#'
#' # Plot gene usage for a specific chain
#' plot_gene_usage(
#'   vdj_sce,
#'   gene_cols = "v_gene",
#'   chain = c("IGH", "IGK")
#' )
#'
#' # Create a heatmap
#' plot_gene_usage(
#'   vdj_so,
#'   gene_cols = "v_gene",
#'   type = "heatmap"
#' )
#'
#' # Plot paired usage of V(D)J segments
#' plot_gene_usage(
#'   vdj_sce,
#'   gene_cols = c("v_gene", "j_gene"),
#' )
#'
#' # Specify colors to use for each cell cluster
#' plot_gene_usage(
#'   vdj_so,
#'   gene_cols = "v_gene",
#'   cluster_col = "orig.ident",
#'   plot_colors = c(avid_2 = "blue", avid_1 = "green")
#' )
#'
#' # Specify order to use for plotting cell clusters
#' plot_gene_usage(
#'   vdj_sce,
#'   gene_cols = "v_gene",
#'   cluster_col = "orig.ident",
#'   plot_lvls = c("avid_2", "avid_1")
#' )
#'
#' # Specify certain V(D)J genes to include in plot
#' plot_gene_usage(
#'   vdj_so,
#'   gene_cols = "v_gene",
#'   vdj_genes = c("IGKV5-43", "IGLV1", "IGHV1-64")
#' )
#'
#' # Specify the number of top V(D)J genes to include in plot
#' plot_gene_usage(
#'   vdj_sce,
#'   gene_cols = "v_gene",
#'   n_genes = 10
#' )
#'
#' # Plot the frequency of each V(D)J segment instead of percent
#' plot_gene_usage(
#'   vdj_so,
#'   gene_cols = "v_gene",
#'   yaxis = "frequency"
#' )
#'
#' @export
plot_gene_usage <- function(input, gene_cols, cluster_col = NULL,
                            group_col = NULL, chain = NULL, type = "bar",
                            plot_colors = NULL, vdj_genes = NULL, n_genes = 20,
                            plot_lvls = names(plot_colors), yaxis = "percent",
                            chain_col = "chains", sep = ";", ...) {

  # Check inputs
  ## NEED TO CHECK WHETHER GENE_COLS INPUT IS COMPATIBLE WITH TYPE ##
  # heat: 1 or 2
  # bar: 1
  # circos: 2
  typs <- c("heatmap", "bar", "circos")

  if (!type %in% typs) {
    stop("type must be one of the following: ", paste(typs, collapse = ", "))
  }

  if (length(gene_cols) > 2) {
    stop("Cannot specify more than two values for gene_cols")
  }

  if (!yaxis %in% c("percent", "frequency")) {
    stop("yaxis must be either 'percent' or 'frequency'")
  }

  .chk_group_cols(cluster_col, group_col)

  if (length(gene_cols) > 1 && !is.null(group_col)) {
    warning(
      "The group_col argument can only be used when a single column ",
      "is passed to the gene_cols argument."
    )

    group_col <- NULL
  }

  if (identical(type, "circos")) yaxis <- "frequency"

  # Set y-axis
  usage_col <- ifelse(identical(yaxis, "frequency"), "freq", "pct")

  # Calculate gene usage
  plt_dat <- calc_gene_usage(
    input       = input,
    gene_cols   = gene_cols,
    cluster_col = c(cluster_col, group_col),
    chain       = chain,
    chain_col   = chain_col,
    sep         = sep
  )

  plt_dat <- dplyr::filter(plt_dat, dplyr::if_all(
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

    # if vdj_genes are provided and two gene_cols are provided, gene pairs
    # containing at least one of the specified genes will be included
    plt_dat <- dplyr::filter(plt_dat, dplyr::if_any(
      dplyr::all_of(gene_cols),
      ~ .x %in% top_genes
    ))

  # If vdj_genes not provided, identify top n_genes for each cluster based on
  # usage
  } else {
    top_genes <- plt_dat
    top_genes <- dplyr::group_by(top_genes, !!!syms(c(cluster_col, group_col)))

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

    # Filter for top gene pairs
    # filter after identifying top genes, instead of at the same time, since
    # want the same top pairs to be plotted for all samples
    plt_dat <- dplyr::filter(plt_dat, dplyr::if_all(
      dplyr::all_of(gene_cols),
      ~ .x %in% top_genes
    ))
  }

  # Order clusters based on plot_lvls
  lvl_col <- dplyr::if_else(is.null(group_col), cluster_col, group_col)
  plt_dat <- .set_lvls(plt_dat, lvl_col, plot_lvls)

  # Plot usage for single gene column
  if (length(gene_cols) == 1) {
    res <- .plot_single_usage(
      plt_dat,
      gn_col   = gene_cols,
      dat_col  = usage_col,
      clst_col = cluster_col,
      grp_col  = group_col,
      typ      = type,
      clrs     = plot_colors,
      ax_ttl   = yaxis,
      ...
    )

    return(res)
  }

  # Plot paired usage for two gene columns
  res <- .plot_paired_usage(
    plt_dat,
    gn_cols  = gene_cols,
    dat_col  = usage_col,
    clst_col = cluster_col,
    typ      = type,
    clrs     = plot_colors,
    ax_ttl   = yaxis,
    ...
  )

  res
}

#' Plot usage for single gene columns
#'
#' @param df_in data.frame
#' @param gn_col Column containing genes
#' @param dat_col Column containing gene usage metric
#' @param clst_col Column containing clusters
#' @param grp_col Column with IDs to use for grouping clusters
#' @param typ Plot type
#' @param clrs Plot colors
#' @param ax_ttl Y-axis title
#' @param order Should genes be ordered based on usage
#' @param ... Additional parameters to pass to plotting function
#' @noRd
.plot_single_usage <- function(df_in, gn_col, dat_col, type, clst_col = NULL,
                               grp_col = NULL, typ, clrs = NULL, ax_ttl,
                               order = TRUE, ...) {

  # Order plot levels
  if (order) {
    df_in <- .order_lvls(df_in, gn_col, dat_col, max, decreasing = TRUE)
  }

  # Set common arguments
  plt_args <- list(clrs = clrs, ...)

  if (identical(typ, "bar") || !is.null(grp_col)) {
    plt_args$x <- gn_col
    plt_args$y <- dat_col
  }

  # Create grouped boxplot
  if (!is.null(grp_col)) {
    plt_args$df_in         <- df_in
    plt_args$.color        <- plt_args$.fill <- grp_col
    plt_args$outlier.color <- plt_args$outlier.color %||% NA
    plt_args$alpha         <- plt_args$alpha %||% 0.5

    res <- purrr::lift_dl(.create_boxes)(plt_args)

    res <- res +
      ggplot2::geom_jitter(
        position = ggplot2::position_jitterdodge(jitter.width = 0.05)
      ) +
      labs(y = ax_ttl) +
      theme(
        legend.position = "right",
        axis.text.x     = ggplot2::element_text(angle = 45, hjust = 1)
      )

    return(res)
  }

  # Create bargraph
  if (identical(typ, "bar")) {
    df_in <- dplyr::filter(df_in, !!sym(dat_col) > 0)

    plt_args$df_in <- df_in
    plt_args$.fill <- clst_col
    plt_args$y_ttl <- ax_ttl

    res <- purrr::lift_dl(.create_bars)(plt_args)

    return(res)
  }

  # Create heatmap
  plt_args$df_in <- df_in
  plt_args$x     <- clst_col
  plt_args$y     <- gn_col
  plt_args$.fill <- dat_col
  plt_args$ttl   <- ax_ttl

  res <- purrr::lift_dl(.create_gg_heatmap)(plt_args)

  res
}

#' Plot paired usage for two gene columns
#'
#' @param df_in data.frame
#' @param gn_cols Columns containing genes
#' @param dat_col Column containing gene usage metric
#' @param clst_col Column containing clusters
#' @param typ Plot type
#' @param clrs Plot colors
#' @param ax_ttl Y-axis title
#' @param order Should genes be ordered based on usage
#' @param ... Additional parameters to pass to plotting function
#' @noRd
.plot_paired_usage <- function(df_in, gn_cols, dat_col, clst_col = NULL,
                               typ, clrs = NULL, ax_ttl, order = TRUE, ...) {

  # Create plot when two columns are passed to gene_col
  # group_split will order list based on the sorted column
  grps <- NULL

  if (!is.null(clst_col)) {
    df_in <- dplyr::group_by(df_in, !!sym(clst_col))

    grps <- df_in[[clst_col]]
    grps <- sort(unique(grps))
  }

  res <- dplyr::group_split(df_in)
  names(res) <- grps

  # Format data for plotting
  # sort level order
  res <- map(res, ~ {
    .x <- dplyr::select(.x, dplyr::all_of(c(gn_cols, dat_col)))
    d  <- tidyr::expand(.x, !!!syms(gn_cols))
    d  <- dplyr::left_join(d, .x, by = gn_cols)

    d <- dplyr::mutate(
      d,
      !!sym(dat_col) := tidyr::replace_na(!!sym(dat_col), 0)
    )

    if (order) {
      d <- dplyr::mutate(
        d,
        dplyr::across(all_of(gn_cols), ~ {
          decr <- dplyr::cur_column() == gn_cols[1]

          .order_lvls(.x, freq, max, decreasing = decr)
        })
      )
    }

    d
  })

  plt_args <- list(clrs = clrs, ...)

  # Create circos plot
  # use purrr::lift_dl to pass dots into walk
  if (identical(typ, "circos")) {
    plt_args$symmetric <- FALSE

    iwalk(res, ~ {
      d <- tidyr::pivot_wider(
        .x,
        names_from  = gn_cols[2],
        values_from = dat_col
      )

      d <- tibble::column_to_rownames(d, gn_cols[1])
      d <- as.matrix(d)

      plt_args$mat_in <- d
      plt_args$ttl <- .y

      purrr::lift_dl(.create_circos)(plt_args)
    })

    return(invisible())
  }

  # Create heatmap
  res <- purrr::map(
    res,
    .create_gg_heatmap,
    x     = gn_cols[1],
    y     = gn_cols[2],
    .fill = dat_col,
    clrs  = clrs,
    ttl   = ax_ttl,
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







