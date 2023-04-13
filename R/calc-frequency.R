#' Calculate frequency of cell groups present in object
#'
#' Calculate the frequency of each cell label present in the provided meta.data
#' column. This is useful for comparing the proportion of cells belonging to
#' different samples, cell types, clonotypes, isotypes, etc.
#'
#' @param input Single cell object or data.frame containing V(D)J data. If a
#' data.frame is provided, the cell barcodes should be stored as row names.
#' @param data_col meta.data column containing cell labels to use for
#' calculating frequency. To calculate clonotype frequencies, provide the column
#' containing clonotype IDs, to calculate isotype frequencies provide the column
#' containing cell isotypes. By default the clonotype_id is used for
#' calculations.
#' @param cluster_col meta.data column containing cluster IDs to use for
#' grouping cells when calculating clonotype abundance
#' @param prefix Prefix to add to new columns
#' @param return_df Return results as a data.frame. If set to `FALSE`, results
#' will be added to the input object.
#'
#'
#' @param per_chain If `TRUE` the frequency of each per-chain value will be
#' calculated. If `FALSE` per-chain data will not be parsed and the values
#' present in `data_col` will be used as is.
#' @param chain Chain(s) to use for calculating frequency. Set to `NULL` to
#' include all chains.
#' @param chain_col meta.data column(s) containing chains for each cell
#' @param sep Separator used for storing per-chain V(D)J data for each cell
#' @return Single cell object or data.frame with clonotype frequencies
#' @seealso [plot_frequency()], [plot_clone_frequency()]
#'
#' @examples
#' # Calculate clonotype abundance using all cells
#' res <- calc_frequency(
#'   vdj_so,
#'   data_col = "clonotype_id"
#' )
#'
#' head(slot(res, "meta.data"), 1)
#'
#' # Group cells based on meta.data column before calculating abundance
#' res <- calc_frequency(
#'   vdj_sce,
#'   data_col = "clonotype_id",
#'   cluster_col = "orig.ident"
#' )
#'
#' head(slot(res, "colData"), 1)
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
#' head(slot(res, "meta.data"), 1)
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
calc_frequency <- function(input, data_col, cluster_col = NULL,
                           prefix = paste0(data_col, "_"), return_df = FALSE,

                           per_chain = FALSE,
                           chain = NULL,
                           chain_col = global$chain_col,
                           sep = global$sep) {

  # Check that columns are present in object
  .check_obj_cols(
    input, data_col, cluster_col, chain = chain, chain_col = chain_col
  )

  # Check arguments
  .check_args()

  # Format input data
  vdj_cols <- c(global$cell_col, data_col, cluster_col)

  if (!is.null(chain)) vdj_cols <- c(vdj_cols, chain_col)

  meta <- .get_meta(input)
  vdj  <- dplyr::select(meta, all_of(vdj_cols))

  # Calculate frequency
  vdj <- .calc_freq(
    df_in       = vdj,
    data_cols   = data_col,
    cluster_col = cluster_col,
    chain       = chain,
    chain_col   = chain_col,
    per_chain   = per_chain,
    prefix      = prefix,
    sep         = sep
  )

  freq_clmn <- paste0(prefix, "freq")
  grp_clmn  <- paste0(prefix, "grp")

  if (!per_chain) {
    vdj <- dplyr::mutate(
      vdj, !!sym(grp_clmn) := .calc_freq_grp(!!sym(freq_clmn))
    )
  }

  # Format results
  res <- dplyr::left_join(meta, vdj, by = global$cell_col)

  if (return_df) input <- meta

  res <- .add_meta(input, meta = res)

  res
}

#' Calculate frequency of a cell label
#'
#' @param df_in Input data.frame
#' @param data_cols Column containing data for calculating abundance
#' (e.g. clonotype IDs)
#' @param cluster_col Column(s) containing cluster IDs to use for grouping cells
#' @param chain Chain(s) to use for calculating frequency. Set to `NULL` to
#' include all chains.
#' @param chain_col meta.data column(s) containing chains for each cell
#' @param prefix Prefix to add to output columns
#' @param per_chain Calculate frequency for per-chain values
#' @param return_df Return a summary data.frame where each row is a group, if
#' `FALSE` each row will represent a cell
#' @param nest When `return_df` is `FALSE` should per-chain results be
#' concatenated using `sep`
#' @param sep Separator for storing per-chain data
#' @return data.frame containing clonotype abundances
#' @noRd
.calc_freq <- function(df_in, data_cols, cluster_col = NULL, chain = NULL,
                       chain_col = global$chain_col, prefix = "",
                       per_chain = FALSE, return_df = FALSE, nest = TRUE,
                       na_remove = TRUE, sep = global$sep) {

  # Format input data
  cluster       <- !is.null(cluster_col)
  multi_cluster <- cluster && length(cluster_col) > 1
  include_zeros <- return_df && cluster
  nest_results  <- per_chain && !return_df && nest
  filt_chains   <- !is.null(chain)

  if (filt_chains) vdj_cols <- c(data_cols, chain_col)
  else             vdj_cols <- data_cols

  res <- df_in

  if (na_remove) {
    res <- df_in <- dplyr::filter(df_in, dplyr::if_all(
      dplyr::all_of(data_cols), ~ !is.na(.x)
    ))
  }

  # Count total cells per cluster
  if (cluster) {
    clst_nm <- paste0(cluster_col, collapse = "_")

    if (multi_cluster) {
      res <- tidyr::unite(
        res, !!sym(clst_nm), !!!syms(cluster_col), remove = FALSE
      )
    }

    clst_key    <- dplyr::select(res, -all_of(data_cols))
    clsts       <- pull(clst_key, clst_nm)
    clst_key    <- dplyr::distinct(clst_key, !!!syms(c(clst_nm, cluster_col)))
    clst_counts <- table(clsts)
    clsts       <- unique(clsts)
  }

  res <- fetch_vdj(
    res,
    data_cols    = vdj_cols,
    filter_cells = FALSE,         # cells are already filtered
    unnest       = !filt_chains,
    per_chain    = per_chain,
    sep          = sep
  )

  # Filter by chain
  if (filt_chains) {
    res <- .filter_chains(
      res,
      data_cols  = data_cols,
      chain      = chain,
      chain_col  = chain_col,
      col_names  = "{.col}",
      allow_dups = TRUE,
      empty_val  = "None"
    )

    res <- dplyr::select(res, -!!sym(chain_col))
    res <- tidyr::unnest(res, dplyr::all_of(data_cols))
  }

  # Filter for unique rows, any duplicated rows are cells with the same
  # gene for multiple chains
  # this is NOT essential for calculations, but good practice
  res <- dplyr::distinct(res)

  # Count cells
  grp_cols <- data_cols

  if (cluster) grp_cols <- unique(c(grp_cols, clst_nm, cluster_col))

  tot_cells <- dplyr::n_distinct(df_in[[global$cell_col]])

  res <- dplyr::group_by(res, !!!syms(grp_cols))

  res <- dplyr::mutate(
    res,
    n_cells = tot_cells,
    freq    = dplyr::n_distinct(!!sym(global$cell_col))
  )

  res <- dplyr::ungroup(res)

  # If returning summary data.frame remove cell ids
  if (return_df) {
    res <- dplyr::select(res, -global$cell_col)
    res <- dplyr::distinct(res)
  }

  # Report zeros for missing groups
  # Need to remove original cluster_col
  # Need to make sure clst_nm is present, when cluster_col is length 1 it
  # is the same as clst_nm
  if (include_zeros) {
    res <- dplyr::select(res, -all_of(cluster_col), dplyr::all_of(clst_nm))

    res <- tidyr::complete(
      res,
      !!!syms(c(data_cols, clst_nm)),
      fill     = list(freq = 0, n_cells = tot_cells),
      explicit = FALSE
    )

    res <- dplyr::left_join(res, clst_key, by = clst_nm)
  }

  # Identify shared groups
  if (cluster) {
    res <- dplyr::group_by(res, !!!syms(data_cols))
    res <- dplyr::mutate(res, shared = length(freq[freq > 0]) > 1)
    res <- dplyr::ungroup(res)
  }

  # Calculate percentage used
  if (cluster) {
    res <- dplyr::mutate(
      res, n_cells = as.numeric(clst_counts[as.character(!!sym(clst_nm))])
    )
  }

  res <- dplyr::mutate(res, pct = (.data$freq / .data$n_cells) * 100)
  res <- dplyr::arrange(res, desc(.data$pct), desc(.data$freq))

  # Output columns
  stat_cols <- c("freq", "pct")

  if (cluster) stat_cols <- c(stat_cols, "shared")

  # Nest per-chain data
  # tidyr::chop() followed by .nest_vdj() is faster than dplyr::summarize()
  if (nest_results) {
    res <- tidyr::chop(res, dplyr::all_of(c(data_cols, stat_cols)))
    res <- .nest_vdj(res, stat_cols)
  }

  # Columns to use for joining final results
  # used when final result to user is a summary table
  if (return_df) {
    final_cols <- c(data_cols, cluster_col, "n_cells")

  # used when plotting results
  } else if (!nest_results) {
    final_cols <- c(global$cell_col, "n_cells")

    if (per_chain) final_cols <- c(final_cols, data_cols)

  # used when adding results to object
  } else {
    final_cols <- global$cell_col
  }

  stat_cols <- purrr::set_names(stat_cols, paste0(prefix, stat_cols))

  res <- dplyr::select(res, dplyr::all_of(final_cols), dplyr::all_of(stat_cols))

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
  # SHOULD REVISE TO USE cut()
  labs <- tibble::tibble(
    x   = uniq_x,
    grp = dplyr::ntile(uniq_x, n_grps)
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


#' Fetch top clonotypes
#'
#' WIP
#'
#' @param input Single cell object or data.frame containing V(D)J data. If a
#' data.frame is provided, the cell barcodes should be stored as row names.
#' @param data_cols meta.data column(s) to use for identifying top clonotypes
#' @param cluster_col meta.data column containing cluster IDs to use for
#' grouping cells.
#' @param clones Number of top clonotypes to identify.
#' @param sep Separator used for storing per-chain V(D)J data for each cell
#' @return data.frame containing information for top clonotypes
#' @noRd
fetch_top_clones <- function(input, data_cols, cluster_col = NULL,
                             clones = 10, sep = global$sep) {

  # Check that columns are present in object
  .check_obj_cols(input, data_cols, cluster_col)

  # Check input classes
  .check_args()

}


#' Plot clonotype frequency
#'
#' @param input Single cell object or data.frame containing V(D)J data. If a
#' data.frame is provided, the cell barcodes should be stored as row names.
#' @param data_col meta.data column containing clonotype IDs to use for
#' calculating clonotype abundance
#' @param cluster_col meta.data column containing cluster IDs to use for
#' grouping cells when calculating clonotype frequencies. Clonotypes will be
#' plotted separately for each cluster.
#' @param group_col meta.data column to use for grouping cluster IDs present in
#' cluster_col. This is useful when there are multiple replicates or patients
#' for each treatment condition. This is only applicable when `method` is
#' 'line'.
#' @param clones An integer specifying the number of clonotypes to show, or
#' a vector giving the names of clonotypes to include. If method
#' is set to 'line', this will specify the clonotypes to label.
#' @param method Method to use for plotting, possible values include:
#'
#' - 'bar', create a bargraph
#' - 'line', create a rank-abundance plot
#'
#' @param units Units to plot on the y-axis, either 'frequency' or 'percent'
#'
#'
#' @param plot_colors Character vector containing colors for plotting
#' @param plot_lvls Levels to use for ordering clusters
#' @param trans Transformation to use for plotting data, e.g. 'log10'. By
#' default values are not transformed, refer to [ggplot2::continuous_scale()]
#' @param panel_nrow The number of rows to use for arranging plot panels, use
#' this when separate bar graphs are created for each cell cluster
#' @param panel_scales Should scales for plot panels be fixed or free. This
#' passes a scales specification to ggplot2::facet_wrap, can be 'fixed', 'free',
#' 'free_x', or 'free_y'. 'fixed' will cause panels to share the same scales.
#' Use this when separate bar graphs are created for each cell cluster.
#' @param n_label Location on plot where n label should be added, this can be
#' any combination of the following:
#'
#' - 'corner', display the total number of cells plotted in the top right
#'   corner, the position of the label can be modified by passing `x` and `y`
#'   specifications with the `label_params` argument
#' - 'legend', display the number of cells plotted for each group shown in the
#'   plot legend
#' - 'none', do not display the number of cells plotted
#'
#' @param label_params Named list providing additional parameters to modify
#' clonotype and n label aesthetics, e.g. list(size = 4, color = "red")
#' @param ... Additional arguments to pass to ggplot2, e.g. color, fill, size,
#' linetype, etc.
#' @return ggplot object
#' @importFrom ggrepel geom_text_repel
#' @importFrom methods formalArgs
#' @seealso [calc_frequency()], [plot_frequency()]
#'
#' @examples
#' # Plot clonotype abundance using all cells
#' plot_clone_frequency(
#'   vdj_so,
#'   data_col = "clonotype_id"
#' )
#'
#' # Plot clonotype abundance separately for each cell cluster
#' plot_clone_frequency(
#'   vdj_sce,
#'   cluster_col = "orig.ident"
#' )
#'
#' # Plot the frequency of each clonotype instead of percentage
#' plot_clone_frequency(
#'   vdj_sce,
#'   cluster_col = "orig.ident",
#'   units = "frequency"
#' )
#'
#' # Specify colors to use for each cell cluster
#' plot_clone_frequency(
#'   vdj_so,
#'   cluster_col = "orig.ident",
#'   plot_colors = c(avid_1 = "blue", avid_2 = "red")
#' )
#'
#' # Specify order to use for plotting cell clusters
#' plot_clone_frequency(
#'   vdj_sce,
#'   cluster_col = "orig.ident",
#'   plot_lvls = c("avid_2", "avid_1")
#' )
#'
#' # Specify the number of top clonotypes to plot
#' plot_clone_frequency(
#'   vdj_so,
#'   clones = 5
#' )
#'
#' #' # Create line graph
#' # use clones to set the number of clonotypes to label
#' plot_clone_frequency(
#'   vdj_so,
#'   cluster_col = "orig.ident",
#'   method = "line",
#'   clones = 3
#' )
#'
#' @export
plot_clone_frequency <- function(input, data_col = global$clonotype_col,
                                 cluster_col = NULL, group_col = NULL,
                                 clones = NULL, method = "bar",
                                 units = "percent",
                                 plot_colors = NULL,
                                 plot_lvls = names(plot_colors),
                                 trans = "identity", panel_nrow = NULL,
                                 panel_scales = "free_x", n_label = "corner",
                                 label_params = list(), ...) {

  # Check that columns are present in object
  .check_obj_cols(input, data_col, cluster_col, group_col)

  # Check input classes
  .check_args()

  # Check input values
  .check_possible_values(
    method = c("bar", "line"),
    units  = c("frequency", "percent")
  )

  clones <- clones %||% switch(method, bar = 10, line = 3)

  if (identical(method, "bar") && is.numeric(clones) && clones <= 0) {
    cli::cli_abort("`clones` must be >0")
  }

  if (identical(method, "line") && is.numeric(clones) && clones < 0) {
    cli::cli_abort("`clones` must be >=0")
  }

  # For bargraph allow user to facet plot using cluster_col or group_col
  if (identical(method, "bar")) {
    cluster_col <- cluster_col %||% group_col
    group_col   <- NULL
  }

  .check_group_cols(cluster_col, group_col, uniq = FALSE)

  abun_col <- switch(units, frequency = ".freq", percent = ".pct")
  y_lab    <- .get_axis_label(units)

  # Calculate clonotype abundance
  plt_dat <- calc_frequency(
    input       = input,
    cluster_col = c(cluster_col, group_col),
    data_col    = data_col,
    prefix      = ".",
    return_df   = TRUE
  )

  plt_dat <- tibble::as_tibble(plt_dat, rownames = global$cell_col)
  plt_dat <- dplyr::filter(plt_dat, !is.na(!!sym(data_col)))

  # Save data for n label
  n_lab_dat  <- plt_dat

  # Identify data columns that the user should have access to
  keep_cols <- .get_matching_clmns(plt_dat, c(data_col, cluster_col, group_col))
  keep_cols <- c(cluster_col, group_col, data_col, keep_cols)
  plt_dat   <- dplyr::distinct(plt_dat, !!!syms(keep_cols))

  # Add and format label column for plotting
  plt_dat <- dplyr::mutate(plt_dat, .lab = trim_lab(!!sym(data_col)))

  # Rank by abundance
  if (!is.null(cluster_col)) {
    rnk_cols <- cluster_col

    if (identical(method, "line")) rnk_cols <- c(rnk_cols, group_col)

    plt_dat <- .set_lvls(plt_dat, cluster_col, plot_lvls)
    plt_dat <- dplyr::group_by(plt_dat, !!!syms(rnk_cols))

    plt_dat <- dplyr::mutate(
      plt_dat,
      orig_data        = !!sym(data_col),
      !!sym(data_col) := paste0(!!sym(cluster_col), "_", !!sym(data_col))
    )
  }

  plt_dat <- dplyr::mutate(
    plt_dat, rank = dplyr::row_number(dplyr::desc(!!sym(abun_col)))
  )

  # Identify top clonotypes
  if (is.numeric(clones)) {
    top_clones <- dplyr::slice_min(
      plt_dat,
      order_by  = rank,
      n         = clones,
      with_ties = FALSE
    )

  } else {
    top_clones <- dplyr::filter(plt_dat, !!sym(data_col) %in% clones)
    clones     <- length(clones)
  }

  plt_dat    <- dplyr::ungroup(plt_dat)
  top_clones <- dplyr::ungroup(top_clones)

  # Set plot arguments
  gg_args <- list(
    grp          = cluster_col,
    clrs         = plot_colors,
    nrow         = panel_nrow,
    scales       = panel_scales,
    n_label      = n_label,
    label_params = label_params,
    label_data   = n_lab_dat,
    ...
  )

  # Create bar graph
  if (identical(method, "bar")) {
    plt_labs   <- purrr::set_names(top_clones$.lab, top_clones[[data_col]])
    top_clones <- dplyr::arrange(top_clones, desc(!!sym(abun_col)))
    lvls       <- unique(top_clones[[data_col]])

    top_clones <- .set_lvls(top_clones, data_col, lvls)

    gg_args$df_in  <- top_clones
    gg_args$x      <- data_col
    gg_args$y      <- abun_col
    gg_args$y_ttl  <- y_lab
    gg_args$.fill  <- gg_args$.color <- cluster_col
    gg_args$x_ang  <- 45
    gg_args$x_hjst <- 1
    gg_args$trans  <- trans

    res <- lift(.create_bars)(gg_args)

    # Format clonotype labels
    res <- res +
      ggplot2::scale_x_discrete(labels = plt_labs)

  } else {
    gg_args$df_in     <- plt_dat
    gg_args$fn        <- ggplot2::geom_line
    gg_args$x         <- "rank"
    gg_args$y         <- abun_col
    gg_args$.color    <- cluster_col
    gg_args$.fill     <- FALSE        # fill is not accepted aes for geom_line
    gg_args$grp       <- group_col
    gg_args$trans_y   <- trans
    gg_args$linewidth <- gg_args$linewidth %||% 1

    if (clones > 0) {
      n_params             <- .parse_label_params(label_params)$n
      gg_args$label_params <- .get_uniq_text_args(n_params, "geom_text")
    }

    res <- lift(.create_plot)(gg_args) +
      labs(y = .get_axis_label(units))

    # Add clonotype labels with ggrepel
    lab_args <- label_params

    if (clones > 0) {
      if (!all(n_label == "none")) {
        label_params <- .get_uniq_text_args(label_params, "geom_text_repel")
      }

      label_params$mapping          <- ggplot2::aes(label = .data$.lab)
      label_params$data             <- top_clones
      label_params$nudge_x          <- label_params$nudge_x %||% 10
      label_params$direction        <- label_params$direction %||% "y"
      label_params$segment.colour   <- label_params$segment.colour %||% "black"
      label_params$segment.size     <- label_params$segment.size %||% 0.2
      label_params$segment.alpha    <- label_params$segment.alpha %||% 0.2
      label_params$segment.linetype <- label_params$segment.linetype %||% 2
      label_params$show.legend      <- label_params$show.legend %||% FALSE

      if (!is.null(label_params$size)) {
        label_params$size <- label_params$size / ggplot2::.pt
      }

      res <- res +
        lift(ggrepel::geom_text_repel)(label_params)
    }
  }

  res
}


#' Plot frequency of cell groups present in object
#'
#' Plot the frequency of each cell label present in the provided meta.data
#' column. This is useful for comparing the proportion of cells belonging to
#' different samples, cell types, clonotypes, isotypes, etc.
#'
#' @param input Single cell object or data.frame containing V(D)J data. If a
#' data.frame is provided, the cell barcodes should be stored as row names.
#' @param data_col meta.data column containing cell labels to use for
#' calculating frequency, e.g. cell types, isotypes etc.
#' @param cluster_col meta.data column containing cluster IDs (or patients,
#' treatment conditions, etc.) to use when calculating frequency. Calculations
#' will be performed separately for each cluster.
#' @param group_col meta.data column to use for grouping cluster IDs present in
#' cluster_col. This is useful when there are multiple replicates or patients
#' for each treatment condition.
#' @param method Method to use for plotting when `group_col` is provided,
#' possible values are 'bar' or 'boxplot'
#' @param stack If `TRUE`, stacked bargraphs will be generated, otherwise grouped
#' bargraphs will be generated
#' @param units Units to plot on the y-axis, either 'frequency' or 'percent'
#' @param top To only show the top cell groups present in `data_col`, provide
#' one of the following, all other cells will be labeled using the value
#' provided to the `other_label` argument. If `NULL` this will be automatically
#' set.
#'
#' - Integer specifying the number of top groups to show
#' - Vector specifying the names of cell groups to show
#'
#' @param other_label Label to use for 'other' cells when `top` is specified, if
#' `NULL` all cell groups present in data_col will be displayed on the plot.
#'
#'
#' @param plot_colors Character vector containing colors for plotting
#' @param plot_lvls Levels to use for ordering clusters or groups
#' @param na_color Color to use for missing values
#' @param trans Transformation to use for plotting data, e.g. 'log10'. By
#' default values are not transformed, refer to [ggplot2::continuous_scale()]
#' for more options. Values can only be transformed when stack is `FALSE`
#' @param show_points If `TRUE` data points will be shown on boxplots, the point
#' size can be adjusted using the `point.size` parameter
#' @param show_zeros If `TRUE` cell labels that are missing from a cluster will
#' still be shown on the plot
#' @param n_label Location on plot where n label should be added, this can be
#' any combination of the following:
#'
#' - 'corner', display the total number of cells plotted in the top right
#'   corner, the position of the label can be modified by passing `x` and `y`
#'   specifications with the `label_params` argument
#' - 'axis', display the number of cells plotted for each group shown on the
#'   x-axis
#' - 'legend', display the number of cells plotted for each group shown in the
#'   plot legend
#' - 'none', do not display the number of cells plotted
#'
#' @param p_label If `TRUE` p-values <0.05 will be shown on plot when
#' `group_col` is specified. A named vector can also be passed to include custom
#' labels for different p-value cutoffs,
#' e.g. `c('*' = 0.05, '**' = 0.01, '***' = 0.001)`.
#' When comparing two groups a t-test will be performed, when
#' comparing more than two groups the Kruskal-Wallis test will be performed.
#' p-values are adjusted for multiple testing using the Bonferroni correction.
#' @param p_file File path to save table containing p-values for each
#' comparison.
#' @param label_params Named list providing additional parameters to modify
#' n label aesthetics, e.g. list(size = 4, color = "red")
#' @param ... Additional arguments to pass to ggplot2, e.g. color, fill, size,
#' linetype, etc.
#'
#'
#' @param per_chain If `TRUE` the frequency of each per-chain value will be
#' calculated. If `FALSE` per-chain data will not be parsed and the values
#' present in `data_col` will be used as is.
#' @param chain Chain(s) to use for calculating frequency. Set to `NULL` to
#' include all chains.
#' @param chain_col meta.data column(s) containing chains for each cell
#' @param sep Separator for storing per-chain data
#' @return ggplot object
#' @seealso [calc_frequency()], [plot_clone_frequency()]
#'
#' @examples
#' # Plot frequency of different isotypes
#' plot_frequency(
#'   vdj_so,
#'   data_col = "isotype"
#' )
#'
#' # Plot frequency separately for cell clusters
#' plot_frequency(
#'   vdj_sce,
#'   data_col    = "isotype",
#'   cluster_col = "orig.ident"
#' )
#'
#' # Create grouped bar graphs
#' plot_frequency(
#'   vdj_sce,
#'   data_col    = "isotype",
#'   cluster_col = "orig.ident",
#'   stack       = FALSE
#' )
#'
#' # Plot number of cells on the y-axis
#' plot_frequency(
#'   vdj_so,
#'   data_col    = "seurat_clusters",
#'   cluster_col = "orig.ident",
#'   units       = "frequency"
#' )
#'
#' @export
plot_frequency <- function(input, data_col, cluster_col = NULL,
                           group_col = NULL, method = "bar",
                           stack = NULL, units = "percent",
                           top = NULL, other_label = "other",
                           plot_colors = NULL,
                           plot_lvls = NULL, na_color = "grey80",
                           trans = "identity", show_points = TRUE,
                           show_zeros = TRUE,
                           n_label = NULL, p_label = TRUE, p_file = NULL,
                           label_params = list(),
                           ...,
                           per_chain = FALSE,
                           chain = NULL,
                           chain_col = global$chain_col, sep = global$sep) {

  # Check that columns are present in object
  .check_obj_cols(input, data_col, cluster_col)

  # Check input classes
  .check_args()

  # Check input values
  .check_group_cols(cluster_col, group_col, input)

  .check_possible_values(units = c("frequency", "percent"))

  if (!identical(trans, "identity") || per_chain) stack <- stack %||% FALSE
  else stack <- stack %||% !is.null(cluster_col)

  if (stack && !identical(trans, "identity")) {
    cli::cli_abort("Values can only be transformed when `stack` is `FALSE`")
  }

  if (stack && per_chain) {
    cli::cli_warn(
      "Stacked bargraphs can only be created when `per_chain` is `FALSE`"
    )

    stack <- FALSE
  }

  abun_col <- switch(units, frequency = ".freq", percent = ".pct")

  y_lab <- .get_axis_label(units)

  # Set plot variables
  if      (!is.null(group_col))   x_col <- data_col
  else if (!is.null(cluster_col)) x_col <- cluster_col
  else if (stack)                 x_col <- NULL
  else                            x_col <- data_col

  clr_col <- group_col %||% data_col

  # Calculate frequency
  vdj_cols <- c(global$cell_col, data_col, cluster_col)

  if (!is.null(chain)) vdj_cols <- c(vdj_cols, chain_col)

  meta    <- .get_meta(input)
  plt_dat <- dplyr::select(meta, all_of(vdj_cols))

  plt_dat <- .calc_freq(
    df_in       = plt_dat,
    data_cols   = data_col,
    cluster_col = cluster_col,
    chain       = chain,
    chain_col   = chain_col,
    per_chain   = per_chain,
    nest        = !per_chain,
    na_remove   = FALSE,
    return_df   = FALSE,
    prefix      = ".",
    sep         = sep
  )

  freq_cols <- c("n_cells", ".freq", ".pct", ".shared")

  # Format plot data
  # When per_chain, data_col is included in .calc_freq output
  # This is because we need the unnested data_col values to determine keep_cols
  if (per_chain) meta <- dplyr::select(meta, -all_of(data_col))

  meta    <- dplyr::select(meta, -any_of(freq_cols))
  plt_dat <- dplyr::left_join(meta, plt_dat, by = global$cell_col)

  # Identify data columns that the user should have access to when plotting
  # This should only include columns that do not split data_col or cluster_col
  # i.e. variables that are not present in multiple groups found in data_col or
  # cluster_col
  keep_cols <- .get_matching_clmns(plt_dat, c(data_col, cluster_col))
  keep_cols <- c(cluster_col, data_col, keep_cols)
  plt_dat   <- dplyr::distinct(plt_dat, !!!syms(keep_cols))

  # Rank values in data_col
  # Need to sum .freq so values for 'other' group are combined
  plt_dat <- .set_other_grps(
    plt_dat, data_col, abun_col, plot_lvls = plot_lvls, top = top,
    other_label = other_label, method = max
  )

  keep_cols <- .get_matching_clmns(plt_dat, c(data_col, cluster_col))
  keep_cols <- c(cluster_col, data_col, keep_cols)
  keep_cols <- keep_cols[!keep_cols %in% freq_cols]

  plt_dat <- dplyr::group_by(plt_dat, .data$n_cells, !!!syms(keep_cols))

  plt_dat <- dplyr::summarize(
    plt_dat, dplyr::across(all_of(c(".freq", ".pct")), sum),
    .groups = "drop"
  )

  # Set plot colors
  plot_colors <- .set_colors(
    plt_dat,
    data_col    = clr_col,
    plot_colors = plot_colors,
    plot_lvls   = levels(plt_dat[[clr_col]]),
    other_label = other_label
  )

  # Set n label data
  # if cluster_col is NULL use data_col to calculate n,
  # use .freq for cell counts
  n_lab_dat <- list()

  grp_cols <- cluster_col %||% data_col
  n_col    <- ifelse(is.null(cluster_col), ".freq", "n_cells")

  grp_cols <- c(.n = n_col, grp_cols, group_col)

  n_lab_dat$corner <- dplyr::distinct(plt_dat, !!!syms(grp_cols))
  n_lab_dat$legend <- n_lab_dat$axis <- n_lab_dat$corner

  if (!is.null(group_col)) {
    n_lab_dat$axis <- dplyr::rename(plt_dat, .n = ".freq")

  } else if (!is.null(cluster_col)) {
    n_lab_dat$legend <- dplyr::rename(plt_dat, .n = ".freq")
  }

  # Plot arguments
  gg_args <- list(
    df_in        = plt_dat,
    x            = x_col,
    y            = abun_col,
    .color       = clr_col,
    .fill        = clr_col,
    clrs         = plot_colors,
    y_ttl        = y_lab,
    na_clr       = na_color,
    trans        = trans,
    n_label      = n_label,
    label_params = label_params,
    label_data   = n_lab_dat,
    n_fn         = sum,
    ...
  )

  # Create grouped boxplot
  if (!is.null(group_col)) {
    gg_args$grp         <- group_col
    gg_args$method      <- method
    gg_args$show_points <- show_points
    gg_args$show_zeros  <- show_zeros
    gg_args$p_label     <- p_label
    gg_args$p_file      <- p_file

    res <- lift(.create_grouped_plot)(gg_args)

  # Create bar graph
  } else {
    gg_args$x_ang  <- 45
    gg_args$x_hjst <- 1

    # Set bar position
    # fill in missing values with 0s
    if (stack) {
      gg_pos <- ggplot2::position_stack()

    } else {
      if (show_zeros) {
        plt_dat <- .add_missing_zeros(plt_dat, abun_col, c(x_col, clr_col))
        gg_args$df_in <- plt_dat
      }

      gg_pos <- ggplot2::position_dodge2(preserve = "single")
    }

    gg_args$position <- gg_args$position %||% gg_pos

    res <- lift(.create_bars)(gg_args)
  }

  res
}
