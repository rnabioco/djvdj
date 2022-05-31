#' Cluster cells based on CDR3 sequences
#'
#' @param input Single cell object or data.frame containing V(D)J data. If a
#' data.frame is provided, the cell barcodes should be stored as row names.
#' @param data_col meta.data column containing CDR3 sequences to use for
#' calculating Levenshtein distance.
#' @param chain Chain to use for clustering CDR3 sequences. Cells with more
#' than one of the provided chain will be excluded from the analysis.
#' @param method Method to use for clustering, possible values are:
#'
#' - 'louvain', multi-level optimization of modality implemented with
#' igraph::cluster_louvain()
#' - 'leiden', the Leiden clustering algorithm implemented with
#' igraph::cluster_leiden()
#'
#' @param resolution Resolution (coarseness) of clusters
#' @param k Number of neighbors for k-nearest neighbors algorithm
#' @param chain_col meta.data column containing chains for each cell.
#' @param prefix Prefix to add to graph name
#' @param return_df Return results as a data.frame. If FALSE, results will be
#' added to the input object.
#' @param sep Separator used for storing per cell V(D)J data
#' @param ... Additional parameters to pass to clustering method
#' @return Single cell object or data.frame with clustering results
#' @importFrom uwot umap
#' @importFrom dbscan kNN
#' @importFrom igraph graph_from_data_frame cluster_louvain cluster_leiden membership
#' @export
cluster_seqs <- function(input, data_col = "cdr3", chain, method = "louvain",
                         k = 10, resolution = 0.5, chain_col = "chains",
                         prefix = paste0(data_col, "_"), return_df = FALSE,
                         sep = ";", ...) {

  # Check inputs
  if (!method %in% c("louvain", "leiden")) {
    stop("method must be either 'louvain' or 'leiden'.")
  }

  # Fetch data
  vdj <- .fetch_seqs(
    input,
    seq_col   = data_col,
    chain     = chain,
    chain_col = chain_col,
    sep       = sep
  )

  vdj <- dplyr::select(vdj, all_of(c(CELL_COL, data_col)))

  # Create Levenshtein distance matrix
  seqs     <- purrr::set_names(vdj[[data_col]], vdj[[CELL_COL]])
  lev_dist <- adist(seqs, ignore.case = TRUE)
  lev_dist <- stats::as.dist(lev_dist)

  # Calculate UMAP
  res <- uwot::umap(lev_dist)

  colnames(res) <- c("UMAP_1", "UMAP_2")
  rownames(res) <- vdj[[CELL_COL]]

  res <- tibble::as_tibble(res, rownames = CELL_COL)

  # Run KNN
  knn_res <- dbscan::kNN(lev_dist, k = k)

  make_adj_df <- function(mat) {
    res <- tibble::as_tibble(mat, rownames = "Var1")
    res <- tidyr::pivot_longer(res, -.data$Var1, values_to = "Var2")
    res
  }

  adj_df <- make_adj_df(knn_res$id)
  adj_df <- dplyr::mutate(adj_df, Var2 = vdj[[CELL_COL]][Var2])

  wt_df <- make_adj_df(knn_res$dist)

  adj_df <- dplyr::left_join(
    adj_df, wt_df,
    by     = c("Var1", "name"),
    suffix = c("", "_wt")
  )

  # Create adjacency graph
  adj_grph <- dplyr::select(adj_df, all_of(c("Var1", "Var2")))
  adj_grph <- igraph::graph_from_data_frame(adj_grph, directed = FALSE)

  # Clustering parameters
  clst_args <- list(
    graph   = adj_grph,
    weights = adj_df$Var2_wt,
    ...
  )

  browser()

  # Louvain clustering
  if (identical(method, "louvain")) {
    clst_method <- igraph::cluster_louvain

    clst_args$resolution <- resolution

  # Leiden clustering
  } else if (identical(method, "leiden")) {
    clst_method <- igraph::cluster_leiden

    clst_args$resolution_parameter <- resolution
    clst_args$objective_function   <- clst_args$objective_function %||% "modularity"
  }

  # Run clustering
  clsts <- purrr::lift_dl(clst_method)(clst_args)
  clsts <- igraph::membership(clsts)
  res   <- dplyr::mutate(res, cluster = as.character(clsts))

  # Format results
  res <- dplyr::rename_with(res, ~ paste0(prefix, .x), -all_of(CELL_COL))

  meta <- .get_meta(input)
  res  <- dplyr::left_join(meta, res, by = CELL_COL)

  if (return_df) input <- meta

  res <- .add_meta(input, meta = res)

  res
}


#' Create sequence logos for clusters
#'
#'@importFrom stringr str_trunc
#'@export
plot_seq_motifs <- function(input, data_col = "cdr3", cluster_col = NULL,
                            chain, plot_colors = NULL, plot_lvls = NULL,
                            chain_col = "chains", width = 0.75,
                            align_end = "5", facet_rows = NULL, sep = ";",
                            ...) {

  # Fetch sequences
  seqs <- .fetch_seqs(
    input,
    chain     = chain,
    seq_col   = data_col,
    chain_col = chain_col,
    sep       = sep
  )

  vdj_cols <- c(CELL_COL, data_col, cluster_col)
  seqs     <- dplyr::select(seqs, all_of(vdj_cols))

  # Trim/filter sequences
  trim_fn <- function(x) .trim_seq(x[[data_col]], width, align_end)

  if (!is.null(cluster_col)) {
    seqs <- .set_lvls(seqs, cluster_col, plot_lvls)
    seqs <- split(seqs, seqs[[cluster_col]])
    seqs <- purrr::map(seqs, trim_fn)

  } else {
    seqs <- trim_fn(seqs)
  }

  # Create logos
  res <- ggplot() +
    ggseqlogo::geom_logo(seqs, ...) +
    facet_wrap(~ seq_group, scales = "free", nrow = facet_rows) +
    djvdj_theme()

  if (!is.null(plot_colors)) {
    suppressMessages({
      res <- res +
        scale_fill_manual(values = plot_colors)
    })
  }

  res
}

#' Fetch chain sequences
#'
#' @noRd
.fetch_seqs <- function(input, seq_col, chain, chain_col, sep = ";") {
  res <- fetch_vdj(
    input,
    vdj_cols      = c(seq_col, chain_col),
    unnest        = FALSE,
    filter_cells  = TRUE,
    clonotype_col = seq_col,
    sep           = sep
  )

  res <- .filter_chains(
    res,
    vdj_cols  = seq_col,
    chain     = chain,
    chain_col = chain_col,
    col_names = "{.col}",
    allow_dups = FALSE
  )

  res <- dplyr::select(res, -all_of(chain_col))
  res <- tidyr::unnest(res, cols = all_of(seq_col))
  res <- dplyr::filter(res, !is.na(!!sym(seq_col)))

  res
}

#' Trim and filter sequences
#'
#' @noRd
.trim_seq <- function(seqs, width = 0.75, end = "5") {

  lens <- nchar(seqs)

  if (width < 1) {
    pct   <- dplyr::percent_rank(lens)
    width <- min(lens[pct >= width])
  }

  res <- seqs[lens >= width]

  s   <- ifelse(end == "5", "right", "left")
  res <- stringr::str_trunc(res, width, side = s, ellipsis = "")

  res
}






# DBScan clustering
# dbscan_clsts <- as.dist(vdj_dist) %>%
#   dbscan::dbscan(eps = 5)
#
# dbscan_clsts <- as.dist(vdj_dist) %>%
#   dbscan::sNNclust(k = 10, eps = 5)
#
# u <- u %>%
#   mutate(dbscan_clsts = as.character(dbscan_clsts$cluster))

# KNN by hand
# knn_res <- tibble()
#
# k <- 10
#
# for (i in seq_len(nrow(vdj_dist))) {
#   nns <- vdj_dist[i, ]
#   nns <- head(sort(nns), k)
#   nns <- names(nns)
#
#   df <- tibble(
#     Var1 = rownames(vdj_dist)[i],
#     Var2 = list(nns)
#   )
#
#   knn_res <- bind_rows(knn_res, df)
# }
#
# knn_res <- unnest(knn_res, Var2)
#
# adj_grph <- knn_res %>%
#   igraph::graph_from_data_frame(directed = FALSE)

# RANN KNN - NOT CORRECT
# knn_res <- RANN::nn2(vdj_dist, k = 3)
#
# adj_df <- knn_res$nn.idx
#
# colnames(adj_df) <- as.character(seq_len(ncol(adj_df)))
#
# adj_df <- as_tibble(adj_df)
#
# adj_grph <- adj_df %>%
#   pivot_longer(all_of(colnames(.)[-1]), values_to = "Var2") %>%
#   select(Var1 = `1`, Var2) %>%
#   igraph::graph_from_data_frame(directed = FALSE)

# ORIGINAL FUNCTION
# https://github.com/rnabioco/djvdj/blame/e9e201d0f3f54702fe74f51858f82c21124dd3b5/R/utils.R
# # Extract sequences
# # Only include cells with VDJ data
# seqs <- Seurat::FetchData(sobj_in, cdr3_col)
# seqs <- tibble::rownames_to_column(seqs, "cell_id")
# seqs <- na.omit(seqs)
#
# # Select chains to used for calculating distance
# if (!is.null(use_chains)) {
#   re       <- stringr::str_c("(?<=", use_chains, ":)[A-Z]+")
#   seq_cols <- stringr::str_c("V", seq_along(use_chains))
#
#   seqs <- purrr::map2(re, seq_cols, ~ {
#     n_col <- dplyr::sym(.y)
#     c_col <- dplyr::sym(cdr3_col)
#
#     dplyr::mutate(seqs, !!n_col := stringr::str_extract(!!c_col, .x))
#   })
#
#   seqs <- purrr::reduce(seqs, dplyr::left_join, by = c("cell_id", cdr3_col))
#   seqs <- na.omit(seqs)
#   seqs <- dplyr::mutate(seqs, seq_c = stringr::str_c(!!!dplyr::syms(seq_cols)))
#
# } else {
#   seqs <- mutate(
#     seqs,
#     seq_c = stringr::str_extract_all(!!dplyr::sym(cdr3_col), "(?<=:)[A-Z]+"),
#     seq_c = purrr::map(seq_c, purrr::reduce, stringr::str_c)
#   )
# }
#
# # Create Levenshtein distance matrix
# seqs     <- purrr::set_names(seqs$seq_c, seqs$cell_id)
# vdj_dist <- adist(seqs)
#
# # Create nearest neighbors graph
# # Add graph this way or error thrown due to differing number of cells
# vdj_snn  <- Seurat::FindNeighbors(vdj_dist, distance.matrix = T)
# snn_name <- stringr::str_c(prefix, "snn")
#
# sobj_in@graphs[[snn_name]] <- vdj_snn$snn
#
# # Find clusters
# res <- Seurat::FindClusters(
#   object     = sobj_in,
#   resolution = resolution,
#   graph.name = snn_name,
#   ...
# )
