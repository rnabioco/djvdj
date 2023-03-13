#' Cluster cells based on CDR3 sequences
#'
#' @param input Single cell object or data.frame containing V(D)J data. If a
#' data.frame is provided, the cell barcodes should be stored as row names.
#' @param data_col meta.data column containing sequences to use for calculating
#' Levenshtein distance.
#' @param chain Chain to use for clustering sequences. Cells with more than one
#' of the provided chain will be excluded from the analysis. If NULL, sequences
#' for cells with multiple chains will be concatenated.
#' @param method Method to use for clustering, possible values include:
#'
#' - 'louvain', multi-level optimization of modality implemented with
#' [igraph::cluster_louvain()]
#' - 'leiden', the Leiden clustering algorithm implemented with
#' [igraph::cluster_leiden()]
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
#' @param run_umap Should the Uniform Manifold Approximation and Projection
#' (UMAP) dimensional reduction method be performed. This will add UMAP
#' coordinates to the meta.data.
#' @param chain_col meta.data column containing chains for each cell.
#' @param prefix Prefix to add to graph name
#' @param return_df Return results as a data.frame. If FALSE, results will be
#' added to the input object.
#' @param sep Separator used for storing per cell V(D)J data
#' @param ... Additional parameters to pass to clustering method
#' @return Single cell object or data.frame with clustering results
#' @importFrom uwot umap
#' @importFrom dbscan kNN
#' @importFrom igraph graph_from_data_frame cluster_louvain cluster_leiden
#' membership
#' @importFrom Biostrings stringDist
#' @seealso [plot_motifs()]
#'
#' @examples
#' # Cluster cells based on CDR3 amino acid sequences and plot results
#' res <- cluster_sequences(
#'   vdj_so,
#'   data_col = "cdr3"
#' )
#'
#' plot_scatter(
#'   res,
#'   x = "cdr3_UMAP_1",
#'   y = "cdr3_UMAP_2",
#'   data_col = "cdr3_cluster_0.5"
#' )
#'
#' # Cluster cells based on sequences for a specific chain
#' res <- cluster_sequences(
#'   vdj_sce,
#'   data_col = "cdr3",
#'   chain    = "IGK"
#' )
#'
#' # Use Levenschtein distance for clustering
#' res <- cluster_sequences(
#'   vdj_so,
#'   data_col    = "cdr3",
#'   dist_method = "levenshtein"
#' )
#'
#' # Cluster cells using the Leiden algorithm
#' res <- cluster_sequences(
#'   vdj_so,
#'   data_col = "cdr3",
#'   method   = "leiden"
#' )
#'
#' # Adjust clustering resolution
#' res <- cluster_sequences(
#'   vdj_so,
#'   data_col   = "cdr3",
#'   resolution = c(0.5, 1)
#' )
#'
#' @export
cluster_sequences <- function(input, data_col = "cdr3", chain = NULL,
                              method = "louvain", resolution = 0.5, k = 10,
                              dist_method = NULL, run_umap = TRUE,
                              chain_col = global$chain_col,
                              prefix = paste0(data_col, "_"),
                              return_df = FALSE, sep = global$sep, ...) {

  # Check that columns are present in object
  .check_obj_cols(
    input, data_col, chain = chain, chain_col = chain_col
  )

  # Check input classes
  .check_args()

  # Check input values
  mets <- c("louvain", "leiden")

  if (!method %in% mets) cli::cli_abort("`method` must be {.or {mets}}")

  # Fetch data
  vdj <- .fetch_seqs(
    input,
    seq_col   = data_col,
    chain     = chain,
    chain_col = chain_col,
    sep       = sep
  )

  vdj  <- dplyr::distinct(vdj, !!!syms(c(global$cell_col, data_col)))
  seqs <- vdj[[data_col]]
  seqs <- unique(sort(seqs))

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

  dist_mat <- lift(Biostrings::stringDist)(dist_args)

  # Calculate UMAP
  if (run_umap) {
    umap_res <- uwot::umap(dist_mat)

    colnames(umap_res) <- paste0(prefix, c("UMAP_1", "UMAP_2"))

    umap_res <- tibble::as_tibble(umap_res)
  }

  # Run KNN
  knn_res <- dbscan::kNN(dist_mat, k = k)

  make_adj_df <- function(mat) {
    res <- tibble::as_tibble(mat, rownames = "Var1")
    res <- tidyr::pivot_longer(res, -"Var1", values_to = "Var2")
    res
  }

  adj_df <- make_adj_df(knn_res$id)

  adj_df <- dplyr::mutate(
    adj_df,
    Var1 = seqs[as.integer(.data$Var1)],
    Var2 = seqs[.data$Var2]
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

    clst_args$objective_function <- clst_args$objective_function %||%
      "modularity"
  }

  # Cluster sequences
  resolution <- purrr::set_names(
    resolution,
    paste0(prefix, "cluster_", resolution)
  )

  res <- purrr::imap_dfc(resolution, ~ {
    clst_args[rsln_arg] <- .x

    clsts <- lift(clst_method)(clst_args)
    clsts <- igraph::membership(clsts)

    tibble(!!sym(.y) := as.character(clsts))
  })

  # Combine cluster and UMAP results, add sequences
  if (run_umap) res <- dplyr::bind_cols(umap_res, res)

  res <- dplyr::mutate(res, !!sym(data_col) := seqs)

  # Join by CDR3 to match clusters with cell IDs
  res <- dplyr::left_join(vdj, res, by = data_col)
  res <- dplyr::select(res, -all_of(data_col))

  # Format results
  meta <- .get_meta(input)
  res  <- dplyr::left_join(meta, res, by = global$cell_col)

  if (return_df) input <- meta

  res <- .add_meta(input, meta = res)

  res
}

#' Plot CDR3 sequence motifs
#'
#' @param input Single cell object or data.frame containing V(D)J data. If a
#' data.frame is provided, the cell barcodes should be stored as row names.
#' @param data_col meta.data column containing sequences to use for plotting.
#' @param cluster_col meta.data column containing cluster IDs to use for
#' grouping cells.
#' @param chain Chain to use for plotting sequences. Cells with more than one
#' of the provided chain will be excluded from the analysis.
#' @param plot_colors Character vector containing colors for plotting
#' @param plot_lvls Character vector containing levels for ordering
#' @param chain_col meta.data column containing chains for each cell.
#' @param width Integer specifying how many residues to include, sequences
#' longer than width will be get trimmed based on the align_end argument,
#' sequences shorter than width will get removed. If a fraction is provided,
#' the width cutoff is set based on percent rank, i.e. a value of 0.75 would
#' select a width where at least 75% of sequences are longer than the cutoff.
#' @param align_end End to use for aligning sequences, specify '5' or '3' to
#' align sequences at the 5' or 3' end when plotting.
#' @param panel_nrow The number of rows to use for arranging plot panels
#' @param panel_scales Should scales for plot panels be fixed or free. This
#' passes a scales specification to `ggplot2::facet_wrap()`, can be 'fixed', 'free',
#' 'free_x', or 'free_y'. 'fixed' will cause panels to share the same scales.
#' Use this when separate bar graphs are created for each cell cluster.
#' @param n_label Location on plot where n label should be added, this can be
#' one of the following:
#'
#' - 'corner', display the total number of cells plotted in the top right
#'   corner, the position of the label can be modified by passing `x` and `y`
#'   specifications with the `label_params` argument
#' - 'none', do not display the number of cells plotted
#'
#' @param label_params Named list providing additional parameters to modify
#' n label aesthetics, e.g. list(size = 4, color = "red")
#' @param quiet If `TRUE` messages will not be displayed
#' @param sep Separator used for storing per cell V(D)J data
#' @param ... Additional parameters to pass to [ggseqlogo::geom_logo()]
#' @return ggplot object
#' @importFrom stringr str_trunc
#' @seealso [cluster_sequences()]
#'
#' @examples
#' # Cluster cells based on CDR3 amino acid sequences and plot sequence motifs
#' res <- cluster_sequences(
#'   vdj_so,
#'   data_col = "cdr3"
#' )
#'
#' plot_motifs(
#'   res,
#'   data_col    = "cdr3",
#'   cluster_col = "cdr3_cluster_0.5",
#'   chain       = "IGK"
#' )
#'
#' @export
plot_motifs <- function(input, data_col = global$cdr3_col, cluster_col = NULL,
                        chain, plot_colors = NULL,
                        plot_lvls = names(plot_colors),
                        chain_col = global$chain_col, width = 0.75,
                        align_end = "5", panel_nrow = NULL,
                        panel_scales = "free", n_label = "corner",
                        label_params = list(), quiet = FALSE, sep = global$sep,
                        ...) {

  # Check that columns are present in object
  .check_obj_cols(
    input, data_col, cluster_col, chain = chain, chain_col = chain_col
  )

  # Check input classes
  .check_args()

  # Check input values
  if (width <= 0) cli::cli_abort("`width` must be >0")

  # Fetch sequences
  seqs <- .fetch_seqs(
    input,
    chain     = chain,
    seq_col   = data_col,
    chain_col = chain_col,
    sep       = sep
  )

  vdj_cols <- c(global$cell_col, data_col, cluster_col)
  seqs     <- dplyr::select(seqs, all_of(vdj_cols))

  n_seqs <- nrow(seqs)

  # Trim/filter sequences
  trim_fn <- function(x) .trim_seq(x[[data_col]], width, align_end)

  if (!is.null(cluster_col)) {
    seqs <- .set_lvls(seqs, cluster_col, plot_lvls)
    seqs <- split(seqs, seqs[[cluster_col]])
    seqs <- purrr::map(seqs, trim_fn)
    seqs <- purrr::discard(seqs, is.null)

    plt_n_seqs <- sum(purrr::map_dbl(seqs, length))

  } else {
    seqs <- trim_fn(seqs)
    plt_n_seqs <- length(seqs)
  }

  # Check number of removed sequences
  if (plt_n_seqs == 0) {
    cli::cli_abort(
      "There are no sequences longer than {width} residues,
       try selecting a shorter cutoff for `width`"
    )
  }

  if (plt_n_seqs < n_seqs && !quiet) {
    pct <- 1 - (plt_n_seqs / n_seqs)
    pct <- round(pct * 100, 1)

    cli::cli_alert_info(
      "{n_seqs - plt_n_seqs} sequences ({pct}%) are shorter
       than `width` and were removed",
      wrap = TRUE
    )
  }

  # Set n label data
  # If cluster_col is provided, seqs will be named list with a vector for each
  # cluster
  if (!is.list(seqs)) seqs <- list(seq = seqs)

  n_lab_dat <- tibble::tibble(
    seq_group = names(seqs),
    .n        = purrr::map_dbl(seqs, length)
  )

  # Create logos
  res <- ggplot() +
    ggseqlogo::geom_logo(seqs, ...) +
    djvdj_theme()

  if (!all(n_label == "none")) {
    res <- .add_n_label(
      res, n_lab_dat,
      n_label  = n_label,
      crnr_col = "seq_group",
      n_fn     = sum,
      lab_args = label_params
    )
  }

  if (!is.null(cluster_col)) {
    res <- res +
      ggplot2::facet_wrap(~ seq_group, scales = panel_scales, nrow = panel_nrow)
  }

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
#' @importFrom stringr str_remove_all
#' @noRd
.fetch_seqs <- function(input, seq_col, chain, chain_col, sep = global$sep) {

  per_chain <- !is.null(chain)

  res <- fetch_vdj(
    input,
    data_cols     = c(seq_col, chain_col),
    unnest        = FALSE,
    filter_cells  = TRUE,
    per_chain     = per_chain,
    clonotype_col = seq_col,
    sep           = sep
  )

  if (!per_chain) {
    if (is.null(sep)) {
      cli::cli_abort("`sep` must be provided when `per_chain` is`FALSE`")
    }

    res <- dplyr::mutate(
      res,
      !!sym(seq_col) := stringr::str_remove_all(!!sym(seq_col), sep)
    )

    return(res)
  }

  res <- .filter_chains(
    res,
    data_cols  = seq_col,
    chain      = chain,
    chain_col  = chain_col,
    col_names  = "{.col}",
    allow_dups = FALSE
  )

  # Remove chain_col since it does not get filtered by .filter_chains() and
  # will be included in the results as a list-col
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
#' cutoff is set based on percent rank, i.e. a value of 0.75 would select a
#' width where at least 75% of sequences are longer than the cutoff.
#' @param end End to use for aligning sequences, specify '5' or '3' to
#' align sequences at the 5' or 3' end.
#' @noRd
.trim_seq <- function(seqs, width = 0.75, end = "5") {

  lens <- nchar(seqs)

  if (width < 1) {
    width <- 1 - width
    pct   <- dplyr::cume_dist(lens)
    width <- min(lens[pct > width])
  }

  res <- seqs[lens >= width]

  if (purrr::is_empty(res)) {
    cli::cli_warn(
      "There are no sequences present that are >{width} residues long"
    )

    return(NULL)
  }

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
  dna <- Biostrings::DNA_BASES
  rna <- Biostrings::RNA_BASES
  aa  <- Biostrings::AA_STANDARD

  # Take first n sequences and split into residues
  seqs <- na.omit(seqs)
  seqs <- utils::head(seqs, n)
  nts  <- purrr::map(seqs, strsplit, "")
  nts  <- unlist(nts)

  # Check intersection with known characters
  int <- intersect(nts, c(dna, rna, aa))

  if (length(int) == 0) cli::cli_abort("Could not determine sequence type")

  # Check if any characters overlap with aa and not dna, rna
  int <- setdiff(intersect(nts, aa), c(dna, rna))

  res <- dplyr::case_when(
    length(int) > 0 ~ "aa",
    "U" %in% nts    ~ "rna",
    TRUE            ~ "dna"
  )

  res
}

