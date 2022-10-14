#' Calculate V(D)J segment usage
#'
#' @param input Object containing V(D)J data. If a data.frame is provided, the
#' cell barcodes should be stored as row names.
#' @param data_cols meta.data column containing V(D)J genes identified for each
#' clonotype. If multiple columns are provided, paired usage of genes will be
#' calculated.
#' @param cluster_col meta.data column containing cell clusters to use when
#' calculating gene usage
#' @param chain Chain to use for calculating gene usage. Set to NULL to include
#' all chains.
#' @param chain_col meta.data column containing chains for each cell
#' @param sep Separator used for storing per cell V(D)J data
#' @return data.frame containing gene usage summary
#' @seealso [plot_gene_usage()]
#'
#' @examples
#' # Calculate V(D)J segment usage for all cells
#' calc_gene_usage(
#'   vdj_so,
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
#'   vdj_so,
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
                            chain_col = "chains", sep = ";") {

  # Format input data
  sep_cols <- data_cols

  if (!is.null(chain)) sep_cols <- c(sep_cols, chain_col)

  vdj_cols <- c(CELL_COL, cluster_col, sep_cols)



  # >>> SHOULD THIS BE REPLACED WITH FETCH_VDJ CALL???
  meta <- .get_meta(input)
  meta <- dplyr::select(meta, all_of(vdj_cols))

  meta <- dplyr::filter(meta, dplyr::if_all(
    all_of(data_cols),
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

    clst_key <- dplyr::select(res, -all_of(data_cols))
  }

  # Filter chains
  if (!is.null(chain)) {
    res <- .filter_chains(
      res,
      data_cols = data_cols,
      chain     = chain,
      chain_col = chain_col,
      col_names = "{.col}",
      empty_val = "None"
    )

    res <- dplyr::select(res, -all_of(chain_col))
  }

  res <- tidyr::unnest(res, cols = all_of(data_cols))
  res <- dplyr::distinct(res)

  # Count genes used
  res <- dplyr::group_by(res, !!!syms(data_cols))

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
      all_of(c(data_cols, cluster_col, "n_cells", "freq"))
    )
  }

  res <- dplyr::mutate(res, pct = (.data$freq / .data$n_cells) * 100)
  res <- dplyr::arrange(res, desc(.data$pct))

  res
}


#' Plot V(D)J segment usage
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
#' @param chain Chain to use for calculating gene usage, set to NULL to include
#' all chains
#' @param chain_col meta.data column containing chains for each cell
#' @param method Method to use for plotting, possible values are:
#'
#' - 'bar', create a bargraph, this is the default when a single column is
#' passed to the data_cols argument
#' - 'heatmap', create a heatmap, this is the default when two columns are
#' passed to the data_cols argument
#' - 'circos', create a circos plot, this requires two columns to be provided
#' to the data_cols argument
#'
#' @param plot_colors Character vector containing colors to use for plot. If a
#' bar graph is created this will specify how to color cell clusters. For a
#' heatmap, these colors will be used to generate the color gradient.
#' @param vdj_genes V(D)J genes to plot, if NULL the top genes will be shown
#' @param n_genes Number of top genes to plot based on usage. If cluster_col is
#' provided, top genes will be identified for each cluster.
#' @param plot_lvls Levels to use for ordering clusters
#' @param units Units to plot on the y-axis, either 'frequency' or 'percent'
#' @param sep Separator used for storing per-chain V(D)J data for each cell
#' @param ... Additional arguments to pass to plotting function,
#' [ggplot2::geom_col()] for bargraph, [ggplot2::geom_tile()] for heatmap,
#' [circlize::chordDiagram()] for circos plot
#' @seealso [calc_gene_usage()]
#' @return ggplot object
#'
#' @examples
#' # Plot V(D)J segment usage for all cells
#' plot_gene_usage(
#'   vdj_so,
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
#'   vdj_so,
#'   data_cols = "v_gene",
#'   chain = c("IGH", "IGK")
#' )
#'
#' # Plot gene usage for a specific chain
#' plot_gene_usage(
#'   vdj_sce,
#'   data_cols = "v_gene",
#'   chain = c("IGH", "IGK")
#' )
#'
#' # Create a heatmap
#' plot_gene_usage(
#'   vdj_so,
#'   data_cols = "v_gene",
#'   type = "heatmap"
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
#'   vdj_so,
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
#'   vdj_so,
#'   data_cols = "v_gene",
#'   vdj_genes = c("IGKV5-43", "IGLV1", "IGHV1-64")
#' )
#'
#' # Specify the number of top V(D)J genes to include in plot
#' plot_gene_usage(
#'   vdj_sce,
#'   data_cols = "v_gene",
#'   n_genes = 10
#' )
#'
#' # Plot the frequency of each V(D)J segment instead of percent
#' plot_gene_usage(
#'   vdj_so,
#'   data_cols = "v_gene",
#'   units = "frequency"
#' )
#'
#' @export
plot_gene_usage <- function(input, data_cols, cluster_col = NULL,
                            group_col = NULL, chain = NULL, method = NULL,
                            plot_colors = NULL, vdj_genes = NULL, n_genes = 20,
                            plot_lvls = names(plot_colors), units = "percent",
                            chain_col = "chains", sep = ";", ...) {

  # Check inputs
  paired <- length(data_cols) == 2

  if (is.null(method))             method  <- ifelse(paired, "heatmap", "bar")
  if (identical(method, "circos")) units <- "frequency"

  .chk_usage_args(method, data_cols, group_col, units, paired)
  .chk_group_cols(cluster_col, group_col)

  # Set y-axis
  usage_col <- ifelse(identical(units, "frequency"), "freq", "pct")

  # Calculate gene usage
  plt_dat <- calc_gene_usage(
    input       = input,
    data_cols   = data_cols,
    cluster_col = c(cluster_col, group_col),
    chain       = chain,
    chain_col   = chain_col,
    sep         = sep
  )

  plt_dat <- dplyr::filter(plt_dat, dplyr::if_all(
    dplyr::all_of(data_cols),
    ~ .x != "None"
  ))

  # If vdj_genes provided check plt_dat
  top_genes <- unique(vdj_genes)

  if (!is.null(top_genes)) {
    absent <- unlist(plt_dat[, data_cols], use.names = FALSE)
    absent <- top_genes[!top_genes %in% absent]

    if (identical(absent, top_genes)) {
      stop("None of the provided genes were found")

    } else if (!purrr::is_empty(absent)) {
      absent <- paste0(absent, collapse = ", ")

      warning("Some genes not found: ", absent)
    }

    # if vdj_genes are provided and two data_cols are provided, gene pairs
    # containing at least one of the specified genes will be included
    plt_dat <- dplyr::filter(plt_dat, dplyr::if_any(
      dplyr::all_of(data_cols),
      ~ .x %in% top_genes
    ))

  # If vdj_genes not provided, identify top n_genes for each cluster based on
  # usage
  } else {
    plt_dat <- .filter_top_genes(
      plt_dat,
      dat_col  = usage_col,
      gn_col   = data_cols,
      clst_col = cluster_col,
      n        = n_genes
    )
  }

  # Order clusters based on plot_lvls
  lvl_col <- dplyr::if_else(is.null(group_col), cluster_col, group_col)
  plt_dat <- .set_lvls(plt_dat, lvl_col, plot_lvls)

  # Plot gene usage
  plt_args <- list(
    df_in    = plt_dat,
    gn_col   = data_cols,
    dat_col  = usage_col,
    clst_col = cluster_col,
    typ      = method,
    clrs     = plot_colors,
    ax_ttl   = units,
    ...
  )

  usage_fn <- .plot_paired_usage

  if (!paired) {
    usage_fn         <- .plot_single_usage
    plt_args$grp_col <- group_col
  }

  purrr::lift_dl(usage_fn)(plt_args)
}

#' Plot usage for single gene column
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
.plot_single_usage <- function(df_in, gn_col, dat_col, typ, clst_col = NULL,
                               grp_col = NULL, clrs = NULL, ax_ttl,
                               order = TRUE, ...) {

  # Order plot levels
  if (order) {
    decr  <- !is.null(grp_col)
    g_col <- sym(gn_col)

    df_in <- dplyr::mutate(
      df_in,
      !!g_col := .order_lvls(!!g_col, !!sym(dat_col), max, decreasing = decr)
    )
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
    df_in          <- dplyr::filter(df_in, !!sym(dat_col) > 0)
    plt_args$df_in <- df_in
    plt_args$.fill <- clst_col
    plt_args$y_ttl <- ax_ttl

    res <- purrr::lift_dl(.create_bars)(plt_args)

    return(res)
  }

  # Create heatmap
  plt_args$df_in   <- df_in
  plt_args$x       <- clst_col
  plt_args$y       <- gn_col
  plt_args$.fill   <- dat_col
  plt_args$lgd_ttl <- ax_ttl

  res <- purrr::lift_dl(.create_gg_heatmap)(plt_args)

  res
}

#' Plot paired usage for two gene columns
#'
#' @param df_in data.frame
#' @param gn_col Columns containing genes
#' @param dat_col Column containing gene usage metric
#' @param clst_col Column containing clusters
#' @param typ Plot type
#' @param clrs Plot colors
#' @param ax_ttl Y-axis title
#' @param order Should genes be ordered based on usage
#' @param ... Additional parameters to pass to plotting function
#' @noRd
.plot_paired_usage <- function(df_in, gn_col, dat_col, clst_col = NULL,
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
  # sort level order, ascending order for y
  res <- map(res, ~ {
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

  # Create heatmap
  if (identical(typ, "heatmap")) {
    plt_args$x       <- gn_col[1]
    plt_args$y       <- gn_col[2]
    plt_args$.fill   <- dat_col
    plt_args$lgd_ttl <- ax_ttl

    res <- purrr::imap(res, ~ {
      plt_args$df_in <- .x

      if (add_ttl) plt_args$plt_ttl <- .y

      purrr::lift_dl(.create_gg_heatmap)(plt_args)
    })

    if (length(res) == 1) return(res[[1]])

    return(res)
  }

  # Create circos plot
  # use purrr::lift_dl to pass dots into walk
  plt_args$symmetric <- FALSE

  iwalk(res, ~ {
    d <- tidyr::pivot_wider(
      .x,
      names_from  = gn_col[1],
      values_from = dat_col
    )

    d <- tibble::column_to_rownames(d, gn_col[2])
    d <- as.matrix(d)

    plt_args$mat_in <- d

    if (add_ttl) plt_args$plt_ttl <- .y

    purrr::lift_dl(.create_circos)(plt_args)
  })
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
    top <- purrr::lift_dl(.get_top_genes)(gn_args)
    res <- dplyr::filter(df_in, !!sym(gn_col) %in% top)

    return(res)
  }

  # Filter for paired gene columns
  # will return n top genes for column 1, for column 2 gene to be included
  # must be a top gene or have the highest usage
  gn_args$clst_col <- clst_col

  top_genes <- purrr::map(gn_col, ~ {
    gn_args$gn_col <- .x
    purrr::lift_dl(.get_top_genes)(gn_args)
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

  } else res <- .pull_gns(res, gn_col)

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
.chk_usage_args <- function(typ, gn_cols, grp_col, axis, paired) {

  typs <- c("heatmap", "bar", "circos")

  if (!typ %in% typs) {
    stop("type must be one of the following: ", paste(typs, collapse = ", "))
  }

  if (paired && !typ %in% c("heatmap", "circos")) {
    stop(
      "type must be 'heatmap' or 'circos' when two columns ",
      "are passed to the data_cols argument."
    )
  }

  if (identical(typ, "circos") && !paired) {
    stop(
      "A circos plot can only be generated when two columns ",
      "are passed to the data_cols argument."
    )
  }

  if (length(gn_cols) > 2) {
    stop("Cannot specify more than two values for data_cols")
  }

  if (!axis %in% c("percent", "frequency")) {
    stop("units must be either 'percent' or 'frequency'")
  }

  if (paired && !is.null(grp_col)) {
    stop(
      "The group_col argument can only be used when a single column ",
      "is passed to the data_cols argument."
    )
  }
}






### EXTRA TOP GENE FILTERING CODE ###
#
# grp_clmns <- gene_cols[1]
#
# if (paired) grp_clmns <- c(grp_clmns, cluster_col)
#
# top_genes <- dplyr::group_by(plt_dat, !!!syms(grp_clmns))
#
# top_genes <- dplyr::summarize(
#   top_genes,
#   !!sym(usage_col) := max(!!sym(usage_col)),
#   .groups = "drop"
# )
#
# top_genes <- dplyr::arrange(top_genes, dplyr::desc(!!sym(usage_col)))
#
# if (paired) {
#   top_genes <- dplyr::group_by(top_genes, !!sym(cluster_col))
#
#   top_genes <- dplyr::slice_max(
#     top_genes,
#     order_by  = !!sym(usage_col),
#     n         = n_genes,
#     with_ties = FALSE
#   )
#
#   top_genes <- dplyr::ungroup(top_genes)
#
#   top_genes <- split(top_genes, top_genes[[cluster_col]])
#   top_genes <- purrr::map(top_genes, dplyr::pull, gene_cols[1])
#
#   plt_dat <- dplyr::group_by(plt_dat, !!sym(cluster_col))
#
#   plt_dat <- dplyr::filter(
#     plt_dat,
#     !!sym(gene_cols[1]) %in% top_genes[[unlist(dplyr::cur_group())]]
#   )
#
#   plt_dat <- dplyr::ungroup(plt_dat)
#
# } else {
#   top_genes <- dplyr::slice_max(
#     top_genes,
#     order_by  = !!sym(usage_col),
#     n         = n_genes,
#     with_ties = FALSE
#   )
#
#   top_genes <- top_genes[[gene_cols[1]]]
#
#   plt_dat <- dplyr::filter(plt_dat, !!sym(gene_cols[1]) %in% top_genes)
# }
#
# # Filter for top gene pairs
# # filter after identifying top genes, instead of at the same time, since
# # want the same top pairs to be plotted for all samples
# gn_args$clst_col <- cluster_col
#
# top_genes <- purrr::lift_dl(.get_top_genes)(gn_args)
#
# plt_dat <- dplyr::group_by(plt_dat, !!sym(cluster_col))
#
# plt_dat <- dplyr::filter(plt_dat, dplyr::if_all(
#   dplyr::all_of(gene_cols),
#   ~ .x %in% top_genes[[unlist(dplyr::cur_group())]]
# ))
#
# plt_dat <- dplyr::ungroup(plt_dat)
#
# ORIGINAL VERSION
# top_genes <- plt_dat
# top_genes <- dplyr::group_by(top_genes, !!!syms(c(cluster_col, group_col)))
#
# # If two gene_cols are provided this will select the top gene pairs
# top_genes <- dplyr::slice_max(
#   top_genes,
#   order_by  = !!sym(usage_col),
#   n         = n_genes,
#   with_ties = FALSE
# )
#
# # Do not use dplyr::pull since user can pass multiple gene_cols
# top_genes <- unlist(top_genes[, gene_cols], use.names = FALSE)
# top_genes <- unique(top_genes)
#
# # Filter for top gene pairs
# # filter after identifying top genes, instead of at the same time, since
# # want the same top pairs to be plotted for all samples
# plt_dat <- dplyr::filter(plt_dat, dplyr::if_all(
#   dplyr::all_of(gene_cols),
#   ~ .x %in% top_genes
# ))






