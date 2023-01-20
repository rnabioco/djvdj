#' Calculate V(D)J segment usage
#'
#' @param input Object containing V(D)J data. If a data.frame is provided, the
#' cell barcodes should be stored as row names.
#' @param data_cols meta.data column containing V(D)J genes identified for each
#' clonotype. If multiple columns are provided, paired usage of genes will be
#' calculated.
#' @param cluster_col meta.data column containing cell clusters to use when
#' calculating gene usage
#' @param chain Chain(s) to use for calculating gene usage. Set to `NULL` to
#' include all chains.
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
                            chain_col = global$chain_col,
                            prefix = paste0(data_cols[1], "_"), return_df = FALSE,
                            sep = global$sep) {

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

#' old version
#' @export
old_calc_gene_usage <- function(input, data_cols, cluster_col = NULL, chain = NULL,
                            chain_col = global$chain_col, sep = global$sep) {

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



  # >>> SHOULD THIS BE REPLACED WITH FETCH_VDJ CALL???
  # Currently fetch_vdj() requires clonotype_col for filtering
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

  # Filter for unique rows, since any duplicated rows are cells with the same
  # gene for multiple chains
  # Do this since we want to calculate the fraction of cells with a given gene,
  # NOT the fraction of chains
  res <- tidyr::unnest(res, cols = all_of(data_cols))
  res <- dplyr::distinct(res)

  # Count genes used
  res <- dplyr::group_by(res, !!!syms(data_cols))

  if (!is.null(cluster_col)) {
    res <- dplyr::group_by(res, !!!syms(c(clst_nm, cluster_col)), .add = TRUE)
  }

  res <- dplyr::summarize(
    res,
    n_cells = dplyr::n_distinct(meta[[global$cell_col]]),
    freq    = dplyr::n_distinct(!!sym(global$cell_col)),
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
      names_from  = all_of(clst_nm),
      values_from = "freq",
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
    if (length(cluster_col) > 1) {
      clst_key <- dplyr::distinct(clst_key, !!!syms(c(clst_nm, cluster_col)))

      res <- dplyr::left_join(res, clst_key, by = clst_nm)
    }

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
#' @param genes An integer specifying the number of genes to plot, or
#' a vector giving the names of genes to include.
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
#' @param plot_lvls Levels to use for ordering clusters
#' @param trans Transformation to use when plotting segment usage, e.g.
#' 'log10'. By default values are not transformed, refer to
#' [ggplot2::continuous_scale()] for more options.
#' @param units Units to plot on the y-axis, either 'frequency' or 'percent'
#' @param rotate_labels Should labels on circos plot be rotated to reduce
#' overlapping text
#' @param panel_nrow The number of rows to use for arranging plots when
#' return_list is FALSE
#' @param show_points If `TRUE` data points will be shown on boxplots, the point
#' size can be adjusted using the `point.size` parameter
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
#' @param label_params Named list providing additional parameters to modify
#' n label aesthetics, e.g. list(size = 4, color = "red")
#' @param return_list Should a list of plots be returned, if FALSE plots will be
#' combined and arranged into panels
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
#'   genes = 10
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
                            group_col = NULL, genes = 20, method = NULL,
                            plot_colors = NULL, plot_lvls = NULL, chain = NULL,
                            chain_col = global$chain_col, trans = "identity",
                            units = "percent", rotate_labels = FALSE,
                            panel_nrow = NULL, show_points = TRUE,
                            n_label = NULL, label_params = list(),
                            return_list = FALSE, sep = global$sep, ...) {

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

    res <- lift(.plot_paired_usage)(gg_args)

    return(invisible())
  }

  # Create plot for single usage
  gg_args$grp_col      <- group_col
  gg_args$show_points  <- show_points
  gg_args$n_label      <- n_label
  gg_args$label_params <- label_params

  res <- lift(.plot_single_usage)(gg_args)

  res
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
#' @param trans Method to use for transforming data
#' @param ttl Title for y-axis or legend depending on type of graph
#' @param order Should genes be ordered based on usage
#' @param n_row Number of rows to use for arranging plots
#' @param ... Additional parameters to pass to plotting function
#' @noRd
.plot_single_usage <- function(df_in, gn_col, dat_col, method, clst_col = NULL,
                               grp_col = NULL, lvls = NULL, show_points = TRUE,
                               n_label = NULL, n_row = 1, ttl = dat_col,
                               order = TRUE, ...) {

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

  # Set n label parameters
  psbl_labs <- switch(
    method,
    heatmap = c("none", "axis"),
    bar     = c("none", "corner", "legend", "axis")
  )

  n_label   <- n_label %||% psbl_labs
  n_label   <- n_label[n_label %in% psbl_labs]

  # Set n label data
  n_lab_dat <- list()

  grp_cols  <- c(.n = "n_cells", clst_col, grp_col)

  n_lab_dat$corner <- dplyr::distinct(df_in, !!!syms(grp_cols))

  n_lab_dat$legend <- n_lab_dat$axis <- n_lab_dat$corner

  if ((!is.null(grp_col) || !is.null(clst_col)) && identical(method, "bar")) {
    n_lab_dat$axis <- dplyr::rename(df_in, .n = freq)
  }

  # Set common arguments
  gg_args <- list(n_label = n_label, label_data = n_lab_dat, n_fn = sum, ...)

  if (identical(method, "bar") || !is.null(grp_col)) {
    gg_args$x <- gn_col
    gg_args$y <- dat_col
  }

  # Create grouped boxplot
  if (identical(method, "bar") && !is.null(grp_col)) {
    gg_args$df_in       <- df_in
    gg_args$.color      <- gg_args$.fill <- grp_col
    gg_args$alpha       <- gg_args$alpha %||% 0.5
    gg_args$show_points <- show_points
    gg_args$y_ttl       <- ttl

    res <- lift(.create_boxes)(gg_args)

    res <- res +
      theme(legend.position = "right")

    return(res)

  # Create bargraph
  } else if (identical(method, "bar")) {

    # Remove zeros from data
    df_in <- dplyr::filter(df_in, !!sym(dat_col) > 0)

    gg_args$df_in <- df_in
    gg_args$.fill <- clst_col
    gg_args$y_ttl <- ttl

    res <- lift(.create_bars)(gg_args)

    return(res)
  }

  # Create heatmap
  gg_args$df_in   <- df_in
  gg_args$x       <- clst_col
  gg_args$y       <- gn_col
  gg_args$grp     <- grp_col
  gg_args$.fill   <- dat_col
  gg_args$lgd_ttl <- ttl
  gg_args$nrow    <- n_row
  gg_args$scales  <- "free_x"

  res <- lift(.create_gg_heatmap)(gg_args)

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
#' @importFrom grid unit viewport
#' @importFrom graphics par plot.new
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

      lift(.create_gg_heatmap)(plt_args)
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

    lift(.create_circos)(plt_args)
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
    top <- lift(.get_top_genes)(gn_args)
    res <- dplyr::filter(df_in, !!sym(gn_col) %in% top)

    return(res)
  }

  # Filter for paired gene columns
  # will return n top genes for column 1, for column 2 gene to be included
  # must be a top gene or have the highest usage
  gn_args$clst_col <- clst_col

  top_genes <- purrr::map(gn_col, ~ {
    gn_args$gn_col <- .x
    lift(.get_top_genes)(gn_args)
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
    method = c("heatmap", "bar", "circos"),
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

