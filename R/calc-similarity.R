#' Calculate repertoire similarity
#'
#' @param input Object containing V(D)J data. If a data.frame is provided, the
#' cell barcodes should be stored as row names.
#' @param data_col meta.data column containing values to use for
#' calculating pairwise similarity between clusters, e.g. 'clonotype_id'
#' @param cluster_col meta.data column containing cluster IDs to use for
#' calculating repertoire overlap
#' @param method Method to use for comparing clusters, possible values include:
#'
#' - 'count', count the number of clonotypes overlapping between each cluster
#' - A function that takes two numeric vectors containing counts for each
#' unique value in data_col, e.g. [abdiv::jaccard()]
#'
#' @param chain Chain to use for comparing clusters. To perform calculations
#' using a single chain, the column passed to the data_col argument must contain
#' per-chain data such as CDR3 sequences. Set to NULL to include all chains.
#' @param chain_col meta.data column containing chains for each cell
#' @param prefix Prefix to add to new columns
#' @param return_mat Return a matrix with similarity values. If set to
#' FALSE, results will be added to the input object.
#' @param sep Separator used for storing per-chain V(D)J data for each cell
#' @return Single cell object or data.frame with similarity values
#' @seealso [plot_similarity()], [calc_mds()], [plot_mds()]
#'
#' @examples
#' # Calculate repertoire overlap
#' res <- calc_similarity(
#'   vdj_sce,
#'   data_col    = "clonotype_id",
#'   cluster_col = "orig.ident",
#'   method      = abdiv::jaccard
#' )
#'
#' head(slot(res, "colData"), 1)
#'
#' # Add a prefix to the new columns
#' # this is useful if multiple calculations are stored in the meta.data
#' res <- calc_similarity(
#'   vdj_sce,
#'   data_col    = "clonotype_id",
#'   cluster_col = "orig.ident",
#'   prefix      = "bcr_"
#' )
#'
#' head(slot(res, "colData"), 1)
#'
#' # Return a matrix instead of adding the results to the input object
#' calc_similarity(
#'   vdj_sce,
#'   data_col    = "clonotype_id",
#'   cluster_col = "orig.ident",
#'   return_mat  = TRUE
#' )
#'
#' @export
calc_similarity <- function(input, data_col, cluster_col,
                            method = abdiv::jaccard, chain = NULL,
                            chain_col = global$chain_col, prefix = NULL,
                            return_mat = FALSE, sep = global$sep) {

  # Check that columns are present in object
  .check_obj_cols(
    input, data_col, cluster_col, chain = chain, chain_col = chain_col
  )

  # Check input classes
  .check_args(method = list(Class = list(c("character", "function"))))

  # Check input values
  if (identical(method, "count")) {
    method <- function(x, y) length(x[x > 0 & y > 0])
    prefix <- prefix %||% "count_"
  }

  # If no prefix provided, use method name
  if (is.null(prefix)) {
    prefix <- as.character(substitute(method))
    prefix <- dplyr::last(prefix)
    prefix <- paste0(prefix, "_")
  }

  # filter chains if provided
  # if all chains removed, NA will be returned
  meta <- vdj <- .get_meta(input)

  if (!is.null(chain)) {
    vdj <- fetch_vdj(
      input,
      data_cols     = c(data_col, chain_col),
      clonotype_col = data_col,
      filter_cells  = TRUE,
      unnest        = FALSE,
      sep           = sep
    )

    vdj <- .filter_chains(
      vdj,
      data_cols  = data_col,
      chain      = chain,
      chain_col  = chain_col,
      col_names  = "{.col}",
      allow_dups = FALSE
    )
  }

  # Check number of clusters after filtering
  vdj <- dplyr::filter(vdj, !is.na(!!sym(data_col)))
  vdj <- dplyr::select(vdj, all_of(c(global$cell_col, data_col, cluster_col)))

  n_clsts <- dplyr::n_distinct(vdj[[cluster_col]])

  if (n_clsts < 2) {
    cli::cli_abort("`cluster_col` must contain at least two unique clusters")
  }

  # Count number of occurrences of each value in data_col
  vdj <- dplyr::group_by(vdj, !!!syms(c(cluster_col, data_col)))

  vdj <- dplyr::summarize(
    vdj,
    n = dplyr::n_distinct(!!sym(global$cell_col)), .groups = "drop"
  )

  clsts <- unique(vdj[[cluster_col]])

  vdj <- tidyr::pivot_wider(
    vdj,
    names_from  = all_of(cluster_col),
    values_from = "n",
    values_fill = 0
  )

  # Calculate similarity index
  combs <- utils::combn(clsts, 2, simplify = FALSE)

  res <- purrr::map_dfr(combs, ~ {
    x <- pull(vdj, .x[1])
    y <- pull(vdj, .x[2])

    tibble::tibble(
      Var1 = .x[1],
      Var2 = .x[2],
      sim  = method(x, y)
    )
  })

  # Calculate self similarity
  res_s <- purrr::map_dfr(clsts, ~ {
    d <- pull(vdj, .x)

    tibble::tibble(
      Var1 = .x,
      Var2 = .x,
      sim  = method(d, d)
    )
  })

  # Combine with inverse combinations
  res_i <- dplyr::rename(res, Var1 = "Var2", Var2 = "Var1")
  res   <- dplyr::bind_rows(res, res_i, res_s)

  # Format data.frame
  clmns <- sort(unique(res$Var2))
  res   <- dplyr::arrange(res, .data$Var2)

  res <- tidyr::pivot_wider(
    res,
    names_from  = "Var1",
    values_from = "sim"
  )

  res <- dplyr::select(res, !!sym(cluster_col) := "Var2", all_of(clmns))

  # Return matrix
  if (return_mat) {
    res <- tibble::column_to_rownames(res, cluster_col)

    return(as.matrix(res))
  }

  # Add column prefixes
  res <- dplyr::rename_with(res, ~ paste0(prefix, .x), all_of(clmns))

  # Add results to input
  res <- dplyr::left_join(meta, res, by = cluster_col)

  # Format results
  # join will cause factor levels to be lost, add these back
  clst_col <- pull(meta, cluster_col)

  if (is.factor(clst_col)) {
    res <- dplyr::mutate(
      res,
      !!sym(cluster_col) := factor(
        !!sym(cluster_col),
        levels = levels(clst_col)
      )
    )
  }

  res <- .add_meta(input, meta = res)

  res
}

#' Plot repertoire similarity
#'
#' @param input Single cell object or data.frame containing V(D)J data. If a
#' data.frame is provided, the cell barcodes should be stored as row names.
#' @param data_col meta.data column containing values to use for
#' calculating pairwise similarity between clusters, e.g. 'clonotype_id'
#' @param cluster_col meta.data column containing cluster IDs to use for
#' calculating overlap
#' @param group_col meta.data column to use for grouping cluster IDs present in
#' cluster_col. This is useful when there are multiple replicates or patients
#' for each treatment condition.
#' @param method Method to use for comparing clusters, possible values are:
#'
#' - A function that takes two numeric vectors containing counts for each
#' clonotype in the object, such as most beta diversity functions provided by
#' the abdiv package. This will generate a heatmap.
#' - 'count', count the number of clonotypes overlapping between each cluster,
#' this will generate a heatmap.
#' - 'circos', create a circos plot summarizing the overlap between clusters
#'
#' @param chain Chain to use for comparing clusters. To perform calculations
#' using a single chain, the column passed to the data_col argument must contain
#' per-chain data such as CDR3 sequences. Set to NULL to include all chains.
#' @param chain_col meta.data column containing chains for each cell
#' @param cluster_heatmap If FALSE, rows and columns of heatmap will not be
#' clustered.
#' @param sep Separator used for storing per-chain V(D)J data for each cell
#'
#'
#' @param plot_colors Character vector containing colors for plotting
#' @param plot_lvls Levels to use for ordering clusters
#' @param rotate_labels Should labels on circos plot be rotated to reduce
#' overlapping text
#' @param remove_upper_triangle If TRUE, upper triangle for heatmap will not
#' be shown.
#' @param remove_diagonal If TRUE, diagonal for heatmap will not be shown.
#' @param ... Additional arguments to pass to plotting function,
#' [ComplexHeatmap::Heatmap()] for heatmap, [circlize::chordDiagram()] for
#' circos plot
#' @return heatmap or circos plot
#' @seealso [calc_similarity()], [calc_mds()], [plot_mds()]
#'
#' @examples
#' # Plot repertoire overlap
#' # use clonotype IDs present in 'clonotype_id' column for calculations
#' plot_similarity(
#'   vdj_sce,
#'   data_col    = "clonotype_id",
#'   cluster_col = "orig.ident"
#' )
#'
#' # Specify method to use for calculating repertoire overlap
#' plot_similarity(
#'   vdj_sce,
#'   data_col    = "clonotype_id",
#'   cluster_col = "orig.ident",
#'   method      = abdiv::morisita
#' )
#'
#' # Specify colors to use for heatmap
#' plot_similarity(
#'   vdj_sce,
#'   data_col    = "clonotype_id",
#'   cluster_col = "orig.ident",
#'   plot_color  = c("white", "red")
#' )
#'
#' # Create circos plot
#' plot_similarity(
#'   vdj_sce,
#'   data_col    = "clonotype_id",
#'   cluster_col = "orig.ident",
#'   method      = "circos"
#' )
#'
#' @export
plot_similarity <- function(input, data_col, cluster_col, group_col = NULL,
                            method = abdiv::jaccard,
                            chain = NULL,
                            chain_col = global$chain_col,
                            cluster_heatmap = TRUE,
                            sep = global$sep,
                            plot_colors = NULL,
                            plot_lvls = names(plot_colors),
                            rotate_labels = FALSE,
                            remove_upper_triangle = FALSE,
                            remove_diagonal = remove_upper_triangle,
                            ...) {

  # Check that columns are present in object
  .check_obj_cols(
    input,
    data_col, cluster_col, group_col, chain = chain, chain_col = chain_col
  )

  # Check input classes
  .check_args(method = list(Class = list(c("character", "function"))))

  # Check input values
  is_circ <- identical(method, "circos")

  if (is_circ) method <- "count"

  .check_group_cols(cluster_col, group_col, input)

  # Calculate similarity
  plt_dat <- calc_similarity(
    input       = input,
    data_col    = data_col,
    cluster_col = cluster_col,
    method      = method,
    chain       = chain,
    chain_col   = chain_col,
    prefix      = "",
    return_mat  = TRUE,
    sep         = sep
  )

  # Create circos plot
  # if level order not provided, sort columns/rows
  if (is_circ) {
    grps <- NULL

    if (!is.null(group_col)) {
      meta <- .get_meta(input)
      grps <- dplyr::distinct(meta, !!!syms(c(cluster_col, group_col)))
      grps <- dplyr::arrange(grps, !!sym(cluster_col))
      grps <- purrr::set_names(grps[[group_col]], grps[[cluster_col]])
    }

    all_nms   <- union(colnames(plt_dat), rownames(plt_dat))
    plot_lvls <- plot_lvls %||% sort(all_nms)

    .create_circos(
      plt_dat,
      clrs          = plot_colors,
      lvls          = plot_lvls,
      grps          = grps,
      rotate_labels = rotate_labels,
      ...
    )

    return(invisible())
  }

  # Similarity column
  sim_col <- as.character(substitute(method))
  sim_col <- dplyr::last(sim_col)

  # Create heatmap
  res <- .create_heatmap(
    plt_dat,
    clrs     = plot_colors,
    lvls     = plot_lvls,
    lgd_ttl  = sim_col,
    cluster  = cluster_heatmap,
    rm_upper = remove_upper_triangle,
    rm_diag  = remove_diagonal,
    ...
  )

  res
}

#' Perform multidimensional scaling
#'
#' Calculate MDS coordinates based on a beta diversity metric.
#'
#' @param input Object containing V(D)J data. If a data.frame is provided, the
#' cell barcodes should be stored as row names.
#' @param data_col meta.data column containing values to use for
#' calculating pairwise similarity between clusters, e.g. 'clonotype_id'
#' @param cluster_col meta.data column containing cluster IDs to use for
#' calculating repertoire overlap
#' @param method Method to use for comparing clusters and calculating MDS
#' coordinates, available methods include:
#'
#' - 'jaccard', Jaccard dissimilarity index implemented with [abdiv::jaccard()]
#' - 'horn_morisita', Horn-Morisita index implemented with
#' [abdiv::horn_morisita()]
#'
#' @param chain Chain to use for comparing clusters. To perform calculations
#' using a single chain, the column passed to the data_col argument must contain
#' per-chain data such as CDR3 sequences. Set to NULL to include all chains.
#' @param chain_col meta.data column containing chains for each cell
#' @param prefix Prefix to add to new columns
#' @param return_df Return results as a data.frame. If set to FALSE, results
#' will be added to the input object.
#' @param sep Separator used for storing per-chain V(D)J data for each cell
#' @return Single cell object or data.frame with MDS coordinates
#' @seealso [plot_mds()], [calc_similarity()], [plot_similarity()]
#'
#' @examples
#' # Calculate MDS coordinates
#' res <- calc_mds(
#'   vdj_sce,
#'   data_col = "clonotype_id",
#'   cluster_col = "isotype"
#' )
#'
#' # Calculate MDS coordinates based on IGK CDR3 sequences
#' res <- calc_mds(
#'   vdj_sce,
#'   data_col    = "cdr3",
#'   cluster_col = "isotype",
#'   chain       = "IGK"
#' )
#'
#' # Change the method used for calculating repertoire similarity
#' res <- calc_mds(
#'   vdj_sce,
#'   data_col    = "clonotype_id",
#'   cluster_col = "isotype",
#'   method      = "horn_morisita"
#' )
#'
#' @export
calc_mds <- function(input, data_col, cluster_col, method = "jaccard",
                     chain = NULL, chain_col = global$chain_col, prefix = "",
                     return_df = FALSE, sep = global$sep) {

  # Check for installed packages
  .check_packages("MASS")

  # Check that columns are present in object
  .check_obj_cols(
    input, data_col, cluster_col, chain = chain, chain_col = chain_col
  )

  # Check input classes
  .check_args()

  # Check input values
  mets <- c(
    "jaccard" = abdiv::jaccard, "horn_morisita" = abdiv::horn_morisita
  )

  if (!method %in% names(mets)) {
    cli::cli_abort("`method` must be {.or {names(mets)}}")
  }

  method <- mets[[method]]

  # Calculate similarity
  res <- calc_similarity(
    input,
    data_col    = data_col,
    cluster_col = cluster_col,
    method      = method,
    chain       = chain,
    chain_col   = chain_col,
    return_mat  = TRUE,
    sep         = sep
  )

  # Must have at least 3 clusters
  if (nrow(res) < 3) {
    cli::cli_abort("`cluster_col` must contain at least three unique clusters")
  }

  # Calculate MDS
  clmns <- c("MDS_1", "MDS_2")

  mds_fn <- purrr::quietly(MASS::isoMDS)
  res    <- mds_fn(stats::as.dist(res))
  res    <- res$result$points

  colnames(res) <- clmns

  res <- tibble::as_tibble(res, rownames = cluster_col)

  # Return data.frame
  if (return_df) return(res)

  # Add column prefixes
  res <- dplyr::rename_with(
    res,
    ~ paste0(prefix, .x),
    all_of(clmns)
  )

  # Add results to input
  meta <- .get_meta(input)

  res <- dplyr::left_join(meta, res, by = cluster_col)

  # Format results
  # join will cause factor levels to be lost, add these back
  clst_col <- pull(meta, cluster_col)

  if (is.factor(clst_col)) {
    res <- dplyr::mutate(
      res,
      !!sym(cluster_col) := factor(
        !!sym(cluster_col),
        levels = levels(clst_col)
      )
    )
  }

  res <- .add_meta(input, meta = res)

  res
}

#' Create MDS plot
#'
#' Calculate MDS coordinates based on a beta diversity metric and plot results.
#'
#' @param input Single cell object or data.frame containing V(D)J data. If a
#' data.frame is provided, the cell barcodes should be stored as row names.
#' @param data_col meta.data column containing values to use for
#' calculating pairwise similarity between clusters, e.g. 'clonotype_id'
#' @param cluster_col meta.data column containing cluster IDs to use for
#' calculating overlap
#' @param method Method to use for comparing clusters and calculating MDS
#' coordinates, available methods include:
#'
#' - 'jaccard', Jaccard dissimilarity index implemented with [abdiv::jaccard()]
#' - 'horn_morisita', Horn-Morisita index implemented with
#' [abdiv::horn_morisita()]
#'
#' @param chain Chain to use for comparing clusters. To perform calculations
#' using a single chain, the column passed to `data_col` must contain
#' per-chain data such as CDR3 sequences. Set to `NULL` to include all chains.
#' @param chain_col meta.data column containing chains for each cell
#' @param plot_colors Character vector containing colors for plotting
#' @param plot_lvls Levels to use for ordering clusters
#' @param label_points Label points on plot
#' @param n_label Location on plot where n label should be added, this can be
#' any combination of the following:
#'
#' - 'corner', display the total number of points plotted in the top right
#'   corner, the position of the label can be modified by passing `x` and `y`
#'   specifications with the `label_params` argument
#' - 'none', do not display the number of points plotted
#'
#' @param label_params Named list providing additional parameters to modify
#' n label aesthetics, e.g. list(size = 4, color = "red")
#' @param sep Separator used for storing per-chain V(D)J data for each cell
#' @param ... Additional arguments to pass to [ggplot2::geom_point()]
#' @return ggplot object
#' @seealso [calc_mds()], [calc_similarity()], [plot_similarity()],
#' [MASS::isoMDS()]
#'
#' @examples
#' # Calculate MDS coordinates
#' plot_mds(
#'   vdj_sce,
#'   data_col = "clonotype_id",
#'   cluster_col = "isotype"
#' )
#'
#' # Calculate MDS coordinates based on IGK CDR3 sequences
#' plot_mds(
#'   vdj_sce,
#'   data_col    = "cdr3",
#'   cluster_col = "isotype",
#'   chain       = "IGK"
#' )
#'
#' # Calculate repertoire similarity using the Horn-Morisita index
#' plot_mds(
#'   vdj_sce,
#'   data_col    = "clonotype_id",
#'   cluster_col = "isotype",
#'   method      = "horn_morisita"
#' )
#'
#' @export
plot_mds <- function(input, data_col, cluster_col,
                     method = "jaccard", chain = NULL,
                     chain_col = global$chain_col, plot_colors = NULL,
                     plot_lvls = names(plot_colors), label_points = TRUE,
                     n_label = "none", label_params = list(),
                     sep = global$sep, ...) {

  # Check that columns are present in object
  .check_obj_cols(
    input, data_col, cluster_col, chain = chain, chain_col = chain_col
  )

  # Check input classes
  .check_args()

  # Calculate MDS
  plt_dat <- calc_mds(
    input       = input,
    data_col    = data_col,
    cluster_col = cluster_col,
    method      = method,
    chain       = chain,
    chain_col   = chain_col,
    prefix      = "",
    return_df   = TRUE,
    sep         = sep
  )

  # Create MDS plot
  lab_args <- .get_uniq_text_args(label_params, "geom_text")

  res <- plot_scatter(
    plt_dat,
    x = "MDS_1",
    y = "MDS_2",
    data_col     = cluster_col,
    plot_colors  = plot_colors,
    plot_lvls    = plot_lvls,
    n_label      = n_label,
    label_params = lab_args,
    ...
  )

  if (label_points) {
    lab_args <- .get_uniq_text_args(label_params, "geom_text_repel")

    lab_args$mapping <- ggplot2::aes(label = !!sym(cluster_col))

    if (!is.null(lab_args$size)) lab_args$size <- lab_args$size / ggplot2::.pt

    res <- res +
      .lift(ggrepel::geom_text_repel)(lab_args)
  }

  res
}
