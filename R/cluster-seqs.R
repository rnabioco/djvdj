#' Cluster cells based on CDR3 sequences
#'
#' @param input Single cell object or data.frame containing V(D)J data. If a
#' data.frame is provided, the cell barcodes should be stored as row names.
#' @param data_col meta.data column containing sequences to use for calculating
#' Levenshtein distance.
#' @param chain Chain to use for clustering sequences. Cells with more than one
#' of the provided chain will be excluded from the analysis.
#' @param method Method to use for clustering, possible values are:
#'
#' - 'louvain', multi-level optimization of modality implemented with
#' igraph::cluster_louvain()
#' - 'leiden', the Leiden clustering algorithm implemented with
#' igraph::cluster_leiden()
#'
#' @param resolution Resolution (coarseness) of clusters
#' @param k Number of neighbors for k-nearest neighbors algorithm
#' @param dist_method Method to use for computing distance between sequences.
#' If NULL, distance is calculated for amino acid sequences using the BLOSUM62
#' substitution matrix and levenshtein distance is calculated for nucleotide
#' sequences. Other possible values include:
#'
#' - 'levenshtein'
#' - 'hamming'
#' - The name of a substitution matrix available from the Biostrings package,
#' e.g. 'BLOSUM62'
#'
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
#' @importFrom Biostrings stringDist
#' @export
cluster_seqs <- function(input, data_col = "cdr3", chain, method = "louvain",
                         resolution = 0.5, k = 10, dist_method = NULL,
                         chain_col = "chains", prefix = paste0(data_col, "_"),
                         return_df = FALSE, sep = ";", ...) {

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

  vdj  <- dplyr::select(vdj, all_of(c(CELL_COL, data_col)))
  seqs <- vdj[[data_col]]

  # Calculate distance
  if (is.null(dist_method)) {
    typ <- .detect_seq_type(seqs, n = 100)
    dist_method <- ifelse(typ == "aa", "BLOSUM62", "levenshtein")
  }

  dist_args <- list(x = seqs)

  if (!dist_method %in% c("hamming", "levenshtein")) {
    dist_args$substitutionMatrix <- dist_method
    dist_method <- "substitutionMatrix"
  }

  dist_args$method <- dist_method

  dist_mat <- purrr::lift_dl(Biostrings::stringDist)(dist_args)

  # Calculate UMAP
  res <- uwot::umap(dist_mat)

  colnames(res) <- paste0(prefix, c("UMAP_1", "UMAP_2"))
  rownames(res) <- vdj[[CELL_COL]]

  res <- tibble::as_tibble(res, rownames = CELL_COL)

  # Run KNN
  knn_res <- dbscan::kNN(dist_mat, k = k)

  make_adj_df <- function(mat) {
    res <- tibble::as_tibble(mat, rownames = "Var1")
    res <- tidyr::pivot_longer(res, -.data$Var1, values_to = "Var2")
    res
  }

  adj_df <- make_adj_df(knn_res$id)

  adj_df <- dplyr::mutate(
    adj_df,
    Var1 = vdj[[CELL_COL]][as.integer(Var1)],
    Var2 = vdj[[CELL_COL]][Var2]
  )

  # Create adjacency graph
  adj_grph <- dplyr::select(adj_df, all_of(c("Var1", "Var2")))
  adj_grph <- igraph::graph_from_data_frame(adj_grph, directed = FALSE)

  # Clustering parameters
  clst_args <- list(graph = adj_grph, ...)

  # Louvain clustering
  if (identical(method, "louvain")) {
    clst_method <- igraph::cluster_louvain
    rsln_arg    <- "resolution"

  # Leiden clustering
  } else if (identical(method, "leiden")) {
    clst_method <- igraph::cluster_leiden
    rsln_arg    <- "resolution_parameter"

    clst_args$objective_function <- clst_args$objective_function %||% "modularity"
  }

  # Cluster sequences
  resolution <- purrr::set_names(
    resolution,
    paste0(prefix, "cluster_", resolution)
  )

  clsts <- purrr::imap_dfc(resolution, ~ {
    clst_args[rsln_arg] <- .x

    clsts <- purrr::lift_dl(clst_method)(clst_args)
    clsts <- igraph::membership(clsts)

    tibble(!!sym(.y) := as.character(clsts))
  })

  res <- dplyr::bind_cols(res, clsts)

  # Format results
  meta <- .get_meta(input)
  res  <- dplyr::left_join(meta, res, by = CELL_COL)

  if (return_df) input <- meta

  res <- .add_meta(input, meta = res)

  res
}

#' Create sequence logos for clusters
#'
#' @param input Single cell object or data.frame containing V(D)J data. If a
#' data.frame is provided, the cell barcodes should be stored as row names.
#' @param data_col meta.data column containing sequences to use for plotting.
#' @param cluster_col meta.data column containing cluster IDs to use for
#' grouping cells.
#' @param chain Chain to use for clustering CDR3 sequences. Cells with more
#' than one of the provided chain will be excluded from the analysis.
#' @param plot_colors Character vector containing colors for plotting
#' @param plot_lvls Character vector containing levels for ordering
#' @param chain_col meta.data column containing chains for each cell.
#' @param width Integer specifying how many residues to plot, sequences longer
#' than width will be get trimmed based on the align_end argument, sequences
#' shorter than width will get removed. If a fraction is provided, the width
#' cutoff will be set to include the specified fraction of sequences, e.g. a
#' value of 0.75 would set the cutoff so 75% of sequences are included in the
#' plot.
#' @param align_end End to use for aligning sequences, specify '5' or '3' to
#' align sequences at the 5' or 3' end when plotting.
#' @param facet_rows The number of facet rows for final plot
#' @param sep Separator used for storing per cell V(D)J data
#' @param ... Additional parameters to pass to ggseqlogo::geom_logo()
#' @importFrom stringr str_trunc
#' @export
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
    ggplot2::facet_wrap(~ seq_group, scales = "free", nrow = facet_rows) +
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
#' @param input Single cell object or data.frame containing V(D)J data. If a
#' data.frame is provided, the cell barcodes should be stored as row names.
#' @param seq_col meta.data column containing sequences
#' @param chain Chain to use for clustering sequences. Cells with more
#' than one of the provided chain will be excluded from the analysis.
#' @param chain_col meta.data column containing chains for each cell.
#' @param sep Separator used for storing per cell V(D)J data
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
#' @param seqs Vector of sequences
#' @param width Integer specifying how many residues to include, sequences
#' longer than width will be get trimmed based on the end argument, sequences
#' shorter than width will get removed. If a fraction is provided, the width
#' cutoff will be set to include the specified fraction of sequences, e.g. a
#' value of 0.75 would set the cutoff so 75% of sequences are included in the
#' output.
#' @param end End to use for aligning sequences, specify '5' or '3' to
#' align sequences at the 5' or 3' end.
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

#' Detect sequence type
#'
#' @param seqs Vector of sequences
#' @param n Number of sequences to check
#' @noRd
.detect_seq_type <- function(seqs, n = 100) {

  dna <- c("A", "T", "G", "C")
  rna <- c("A", "U", "G", "C")

  aa  <- c(
    "A", "R", "N", "D", "C", "Q", "E", "G", "H", "I",
    "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"
  )

  # Take first n sequences and split into residues
  seqs <- na.omit(seqs)
  seqs <- head(seqs, n)
  nts  <- purrr::map(seqs, strsplit, "")
  nts  <- unlist(nts)

  # Check intersection with known characters
  int <- intersect(nts, c(dna, rna, aa))

  if (length(int) == 0) stop("Could not determine sequence type")

  # Check if any characters overlap with aa and not dna, rna
  int <- setdiff(intersect(nts, aa), c(dna, rna))

  res <- dplyr::case_when(
    length(int) > 0 ~ "aa",
    "U" %in% nts    ~ "rna",
    TRUE            ~ "dna"
  )

  res
}

# seqs <- vdj_so %>%
#   filter_vdj(chains == "IGH") %>%
#   fetch_vdj(c("cdr3", "chains")) %>%
#   group_by(.cell_id) %>%
#   filter(n() == 1) %>%
#   pull(cdr3) %>%
#   na.omit()











