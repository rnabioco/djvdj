#' Calculate cluster similarity
#'
#' @param input Object containing V(D)J data. If a data.frame is provided, the
#' cell barcodes should be stored as row names.
#' @param cluster_col meta.data column containing cluster IDs to use for
#' calculating repertoire overlap
#' @param method Method to use for comparing clusters, possible values are:
#'
#' - 'count', count the number of clonotypes overlapping between each cluster
#' - A function that takes two numeric vectors containing counts for each
#' unique value in the column provided to the data_col column, e.g.
#' [abdiv::jaccard()]
#'
#' @param data_col meta.data column containing clonotype IDs to use for
#' calculating overlap
#' @param chain Chain to use for calculating gene usage. Set to NULL to include
#' all chains.
#' @param chain_col meta.data column containing chains for each cell
#' @param prefix Prefix to add to new columns
#' @param return_mat Return a matrix with similarity values. If set to
#' FALSE, results will be added to the input object.
#' @param sep Separator used for storing per-chain V(D)J data for each cell
#' @return Single cell object or data.frame with similarity values
#' @importFrom abdiv jaccard
#' @seealso [plot_similarity()]
#'
#' @examples
#' # Calculate repertoire overlap
#' res <- calc_similarity(
#'   vdj_so,
#'   data_col = "clonotype_id",
#'   cluster_col = "orig.ident",
#'   method = abdiv::jaccard
#' )
#'
#' head(res@meta.data, 1)
#'
#' # Add a prefix to the new columns
#' # this is useful if multiple calculations are stored in the meta.data
#' res <- calc_similarity(
#'   vdj_sce,
#'   data_col = "clonotype_id",
#'   cluster_col = "orig.ident",
#'   prefix = "bcr_"
#' )
#'
#' head(res@colData, 1)
#'
#' # Return a matrix instead of adding the results to the input object
#' calc_similarity(
#'   vdj_so,
#'   data_col = "clonotype_id",
#'   cluster_col = "orig.ident",
#'   return_mat = TRUE
#' )
#'
#' @export
calc_similarity <- function(input, data_col, cluster_col, method = abdiv::jaccard,
                            chain = NULL, chain_col = "chains", prefix = NULL,
                            return_mat = FALSE, sep = ";") {

  # Check inputs
  is_counts <- identical(method, "count")

  if (!is.function(method) && !is_counts) {
    stop(
      "method must be 'count' or a function to use for comparing ",
      "clusters, e.g. method = abdiv::jaccard."
    )
  }

  if (is_counts) {
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
      vdj_cols      = c(data_col, chain_col),
      clonotype_col = data_col,
      filter_cells  = TRUE,
      unnest        = FALSE,
      sep           = sep
    )

    vdj <- .filter_chains(
      vdj,
      vdj_cols   = data_col,
      chain      = chain,
      chain_col  = chain_col,
      col_names  = "{.col}",
      allow_dups = FALSE
    )
  }

  # Check number of clusters after filtering
  vdj <- dplyr::filter(vdj, !is.na(!!sym(data_col)))
  vdj <- dplyr::select(vdj, all_of(c(CELL_COL, data_col, cluster_col)))

  n_clsts <- dplyr::n_distinct(vdj[[cluster_col]])

  if (n_clsts < 2) {
    stop("cluster_col must contain at least two unique groups.")
  }

  # Count number of occurrences of each value in data_col
  vdj <- dplyr::group_by(vdj, !!!syms(c(cluster_col, data_col)))

  vdj <- dplyr::summarize(
    vdj,
    n = dplyr::n_distinct(!!sym(CELL_COL)), .groups = "drop"
  )

  clsts <- unique(vdj[[cluster_col]])

  vdj <- tidyr::pivot_wider(
    vdj,
    names_from  = all_of(cluster_col),
    values_from = .data$n,
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
  res_i <- dplyr::rename(res, Var1 = .data$Var2, Var2 = .data$Var1)
  res   <- dplyr::bind_rows(res, res_i, res_s)

  # Format data.frame
  clmns <- sort(unique(res$Var2))
  res   <- dplyr::arrange(res, .data$Var2)

  res <- tidyr::pivot_wider(
    res,
    names_from  = .data$Var1,
    values_from = .data$sim
  )

  res <- dplyr::select(res, !!sym(cluster_col) := .data$Var2, all_of(clmns))

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

#' Plot cluster similarity
#'
#' @param input Single cell object or data.frame containing V(D)J data. If a
#' data.frame is provided, the cell barcodes should be stored as row names.
#' @param data_col meta.data column containing clonotype IDs to use for
#' calculating overlap
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
#' - 'circos', create a circos plot
#'
#' @param chain Chain to use for calculating gene usage. Set to NULL to include
#' all chains.
#' @param chain_col meta.data column containing chains for each cell
#' @param plot_colors Character vector containing colors for plotting
#' @param plot_lvls Levels to use for ordering clusters
#' @param cluster_heatmap If FALSE, rows and columns of heatmap will not be
#' clustered.
#' @param remove_upper_triangle If TRUE, upper triangle for heatmap will not
#' be shown.
#' @param remove_diagonal If TRUE, diagonal for heatmap will not be shown.
#' @param sep Separator used for storing per-chain V(D)J data for each cell
#' @param ... Additional arguments to pass to plotting function,
#' [ComplexHeatmap::Heatmap()] for heatmap, [circlize::chordDiagram()] for
#' circos plot
#' @importFrom abdiv jaccard
#' @seealso [calc_similarity()]
#'
#' @examples
#' # Plot repertoire overlap
#' plot_similarity(
#'   vdj_so,
#'   data_col = "clonotype_id",
#'   cluster_col = "orig.ident"
#' )
#'
#' # Specify method to use for calculating repertoire overlap
#' plot_similarity(
#'   vdj_sce,
#'   data_col = "clonotype_id",
#'   cluster_col = "orig.ident",
#'   method = abdiv::morisita
#' )
#'
#' # Specify colors to use for heatmap
#' plot_similarity(
#'   vdj_so,
#'   data_col = "clonotype_id",
#'   cluster_col = "orig.ident",
#'   plot_color = c("white", "red")
#' )
#'
#' # Create circos plot
#' plot_similarity(
#'   vdj_so,
#'   data_col = "clonotype_id",
#'   cluster_col = "orig.ident",
#'   method = "circos"
#' )
#'
#' @export
plot_similarity <- function(input, data_col, cluster_col, group_col = NULL,
                            method = abdiv::jaccard, chain = NULL,
                            chain_col = "chains", plot_colors = NULL,
                            plot_lvls = names(plot_colors),
                            cluster_heatmap = TRUE,
                            remove_upper_triangle = FALSE,
                            remove_diagonal = remove_upper_triangle, sep = ";",
                            ...) {

  # Check inputs
  is_circ <- identical(method, "circos")

  if (is_circ) method <- "count"

  .chk_group_cols(cluster_col, group_col)

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
      clrs = plot_colors,
      lvls = plot_lvls,
      grps = grps
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
#' @param input Object containing V(D)J data. If a data.frame is provided, the
#' cell barcodes should be stored as row names.
#' @param cluster_col meta.data column containing cluster IDs to use for
#' calculating repertoire overlap
#' @param method Method to use for comparing clusters, possible values are:
#'
#' - A function that takes two numeric vectors containing counts for each
#' unique value in the column provided to the data_col column, e.g.
#' abdiv::jaccard()
#'
#' @param data_col meta.data column containing clonotype IDs to use for
#' calculating overlap
#' @param chain Chain to use for calculating gene usage. Set to NULL to include
#' all chains.
#' @param chain_col meta.data column containing chains for each cell
#' @param prefix Prefix to add to new columns
#' @param return_df Return results as a data.frame. If set to FALSE, results
#' will be added to the input object.
#' @param sep Separator used for storing per-chain V(D)J data for each cell
#' @return Single cell object or data.frame with MDS coordinates
#' @importFrom MASS isoMDS
#' @seealso [plot_mds()], [MASS::isoMDS()]
#' @export
calc_mds <- function(input, data_col, cluster_col, method = abdiv::jaccard,
                     chain = NULL, chain_col = "chains", prefix = "",
                     return_df = FALSE, sep = ";") {

  # Check inputs
  if (!is.function(method)) {
    stop(
      "method must be a function to use for comparing ",
      "clusters, e.g. method = abdiv::jaccard."
    )
  }

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
    stop("cluster_col must contain at least three unique groups.")
  }

  # Calculate MDS
  clmns <- c("MDS_1", "MDS_2")

  mds_fn <- purrr::quietly(MASS::isoMDS)
  res    <- mds_fn(as.dist(res))
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
#' Perform multidimensional scaling and plot results
#'
#' @param input Single cell object or data.frame containing V(D)J data. If a
#' data.frame is provided, the cell barcodes should be stored as row names.
#' @param data_col meta.data column containing clonotype IDs to use for
#' calculating overlap
#' @param cluster_col meta.data column containing cluster IDs to use for
#' calculating overlap
#' @param method Method to use for comparing clusters, possible values are:
#'
#' - A function that takes two numeric vectors containing counts for each
#' clonotype in the object, such as most beta diversity functions provided by
#' the abdiv package. This will generate a heatmap.
#'
#' @param chain Chain to use for calculating gene usage. Set to NULL to include
#' all chains.
#' @param chain_col meta.data column containing chains for each cell
#' @param plot_colors Character vector containing colors for plotting
#' @param plot_lvls Levels to use for ordering clusters
#' @param label_points Label points on plot
#' @param sep Separator used for storing per-chain V(D)J data for each cell
#' @param ... Additional arguments to pass to ([ggplot2::geom_point()]
#' @seealso [calc_mds()]
#' @export
plot_mds <- function(input, data_col, cluster_col,
                     method = abdiv::jaccard, chain = NULL,
                     chain_col = "chains", plot_colors = NULL,
                     plot_lvls = names(plot_colors), label_points = TRUE,
                     sep = ";", ...) {

    # Calculate MDS
    plt_dat <- calc_mds(
      input       = input,
      data_col    = data_col,
      cluster_col = cluster_col,
      method      = method,
      chain       = chain,
      chain_col   = chain_col,
      prefix      = "",
      return_df  = TRUE,
      sep         = sep
    )

    # Create MDS plot
    res <- plot_features(
      plt_dat,
      x = "MDS_1",
      y = "MDS_2",
      feature     = cluster_col,
      plot_colors = plot_colors,
      plot_lvls   = plot_lvls,
      ...
    )

    if (label_points) {
      res <- res +
        ggrepel::geom_text_repel(ggplot2::aes(label = !!sym(cluster_col)))
    }

    res
  }
