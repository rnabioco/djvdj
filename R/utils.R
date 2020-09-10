#' Add VDJ data to Seurat object
#'
#' @param sobj_in Seurat object
#' @param vdj_dir cellranger vdj output directory
#' @param include_chains Only inlude clonotypes that have the indicated chains
#' @param prefix Prefix to add to new meta.data columns
#' @param cell_prefix Prefix to add to cell barcodes
#' @return Seurat object with VDJ data added to meta.data
#' @export
import_vdj <- function(sobj_in, vdj_dir, include_chains = NULL, prefix = "",
                       cell_prefix = "") {

  # Load contigs
  vdj_cols <- c(
    "barcode",    "raw_clonotype_id",
    "v_gene",     "d_gene",
    "j_gene",     "c_gene",
    "productive", "full_length"
  )

  contigs <- readr::read_csv(file.path(vdj_dir, "filtered_contig_annotations.csv"))
  contigs <- dplyr::select(contigs, all_of(vdj_cols))
  contigs <- filter(contigs, productive, full_length)
  contigs <- unique(contigs)
  contigs <- mutate(contigs, barcode = str_c(cell_prefix, barcode))

  contigs <- group_by(contigs, barcode, raw_clonotype_id)

  contigs <- summarize(
    contigs,
    dplyr::across(
      dplyr::ends_with("_gene"),
      ~ purrr::reduce(.x, stringr::str_c, sep = ";")
    ),
    .groups = "drop"
  )

  # Load CDR3 sequences for clonotypes
  ctypes <- readr::read_csv(file.path(vdj_dir, "clonotypes.csv"))

  # Merge clonotype and contig info
  meta_df <- dplyr::left_join(
    contigs, ctypes,
    by = c("raw_clonotype_id" = "clonotype_id")
  )

  # Filter for cells present in sobj_in
  cells   <- Seurat::Cells(sobj_in)
  meta_df <- dplyr::filter(meta_df, barcode %in% cells)
  meta_df <- dplyr::rename(meta_df, clonotype_id = raw_clonotype_id)
  meta_df <- dplyr::select(meta_df, -frequency, -proportion)

  # Filter for clonotypes that include the given chains
  if (!is.null(include_chains)) {
    re <- purrr::map_chr(include_chains, ~ stringr::str_c("(?=.*", .x, ":)"))
    re <- purrr::reduce(re, stringr::str_c)

    meta_df <- dplyr::filter(meta_df, stringr::str_detect(cdr3s_aa, re))
  }

  # Calculate stats
  meta_df <- dplyr::mutate(
    meta_df,
    n_chains = stringr::str_split(cdr3s_aa, ";"),
    n_chains = purrr::map_dbl(n_chains, length)
  )

  meta_df <- dplyr::group_by(meta_df, clonotype_id)
  meta_df <- dplyr::mutate(meta_df, clone_freq = dplyr::n_distinct(barcode))
  meta_df <- dplyr::ungroup(meta_df)
  meta_df <- dplyr::mutate(meta_df, clone_prop = clone_freq / nrow(meta_df))

  # Add meta.data to Seurat object
  meta_df <- tibble::column_to_rownames(meta_df, "barcode")
  meta_df <- dplyr::rename_with(meta_df, ~ stringr::str_c(prefix, .x))

  res <- Seurat::AddMetaData(sobj_in, metadata = meta_df)

  res
}


#' Calculate receptor diversity (inverse Simpson Index)
#'
#' @param sobj_in Seurat object
#' @param clonotype_col meta.data column containing clonotype ids
#' @param cluster_col meta.data column containing cluster ids used for
#' calculating receptor diversity. If cluster_col is omitted, diversity index
#' will be calculated for all clonotypes.
#' @param prefix Prefix to add to new meta.data columns
#' @return Seurat object with inverse Simpson index added to meta.data
#' @export
calc_diversity <- function(sobj_in, clonotype_col = "clonotype_id",
                           cluster_col = NULL, prefix = "") {

  # meta.data
  vdj_cols      <- clonotype_col
  clonotype_col <- dplyr::sym(clonotype_col)
  meta_df       <- tibble::as_tibble(sobj_in@meta.data, rownames = "cell_id")
  vdj_df        <- dplyr::filter(meta_df, !is.na(!!clonotype_col))

  # Count clonotypes
  if (!is.null(cluster_col)) {
    vdj_cols    <- c(cluster_col, vdj_cols)
    cluster_col <- dplyr::sym(cluster_col)
    vdj_df      <- dplyr::group_by(vdj_df, !!cluster_col)
  }

  vdj_df <- dplyr::group_by(vdj_df, !!clonotype_col, .add = T)
  vdj_df <- dplyr::summarize(
    .data   = vdj_df,
    num     = dplyr::n_distinct(cell_id),
    .groups = "drop"
  )

  # Calculate diversity
  if (!is.null(cluster_col)) {
    vdj_df <- dplyr::group_by(vdj_df, !!cluster_col)
  }

  div_col <- str_c(prefix, "diversity")

  vdj_df <- dplyr::mutate(
    vdj_df,
    frac     = num / sum(num),
    sum_frac = sum(frac ^ 2),

    !!dplyr::sym(div_col) := 1 - sum_frac
    # !!dplyr::sym(div_col) := 1 / sum_frac
  )

  vdj_df <- dplyr::select(vdj_df, all_of(c(vdj_cols, div_col)))

  # Add resuts to meta.data
  meta_df <- dplyr::left_join(meta_df, vdj_df, by = vdj_cols)
  meta_df <- tibble::column_to_rownames(meta_df, "cell_id")
  meta_df <- as.data.frame(meta_df)

  res <- Seurat::AddMetaData(sobj_in, metadata = meta_df)

  res
}


#' Calculate Jaccard index for cell identities
#'
#' @param sobj_in Seurat object
#' @param clonotype_col meta.data column containing clonotype ids
#' @param cluster_col meta.data column containing cell clusters
#' @param ref_cluster Cluster id to use as a reference for calculating Jaccard
#' index. If ref_cluster is omitted, Jaccard index will be calculated for all
#' combinations of clusters.
#' @param return_matrix Return matrix instead of Seurat object
#' @param prefix Prefix to add to new meta.data columns
#' @return Seurat object with Jaccard index added to meta.data
#' @export
calc_jaccard <- function(sobj_in, clonotype_col = "clonotype_id", cluster_col,
                         ref_cluster = NULL, prefix = "") {

  # Helper to calculate jaccard index
  calc_jidx <- function(df_in, comparison, clonotype_col) {

    if (length(comparison) != 2) {
      stop("comparison must be a character vector containing two elements")
    }

    uniq_vars     <- unique(comparison)
    clonotype_col <- dplyr::sym(clonotype_col)

    if (length(uniq_vars) == 1) {
      dup_var       <- stringr::str_c(uniq_vars, "_dup")
      comparison[3] <- dup_var
      dup_var       <- dplyr::sym(dup_var)
      uniq_vars     <- dplyr::sym(uniq_vars)

      df_in <- dplyr::mutate(df_in, !!dup_var := !!uniq_vars)
    }

    row_sums <- dplyr::select(df_in, !!clonotype_col, dplyr::all_of(comparison))
    row_sums <- tidyr::pivot_longer(row_sums, cols = -!!clonotype_col)
    row_sums <- dplyr::group_by(row_sums, !!clonotype_col)
    row_sums <- dplyr::summarize(row_sums, row_sum = sum(value), .groups = "drop")
    row_sums <- row_sums$row_sum

    a     <- length(row_sums[row_sums == 2])      # Intersection
    a_b_c <- length(row_sums[row_sums == 1]) + a  # Union

    jaccard <- a / a_b_c

    res <- tibble::tibble(
      Var1 = comparison[1],
      Var2 = comparison[2],
      jaccard
    )

    res
  }

  # Fetch clonotypes and cell identities
  # so_idents   <- Seurat::Idents(sobj_in)
  so_idents   <- Seurat::FetchData(sobj_in, cluster_col)
  so_idents   <- so_idents[, cluster_col]
  so_idents   <- as.character(so_idents)
  uniq_idents <- unique(so_idents)
  uniq_idents <- na.omit(uniq_idents)

  ctypes   <- Seurat::FetchData(sobj_in, clonotype_col)
  vdj_meta <- dplyr::bind_cols(ctypes, idents = so_idents)
  vdj_meta <- filter(vdj_meta, !is.na(!!dplyr::sym(clonotype_col)))

  # Create data.frame for calculating Jaccard index
  j_df <- dplyr::mutate(vdj_meta, num = 1)
  j_df <- dplyr::mutate(j_df, num = 1)

  j_df <- dplyr::group_by(j_df, idents)
  j_df <- dplyr::group_split(j_df)
  j_df <- purrr::map(
    .x          = j_df,
    .f          = pivot_wider,
    names_from  = idents,
    values_from = num,
    values_fn   = list
  )

  j_df <- purrr::map(j_df, tidyr::unnest, cols = -!!dplyr::sym(clonotype_col))
  j_df <- purrr::map(j_df, unique)
  j_df <- purrr::reduce(j_df, dplyr::full_join, by = clonotype_col)
  j_df <- dplyr::mutate_all(j_df, tidyr::replace_na, replace = 0)

  # Create data.frame of comparisons
  # Currently combinations include duplicates
  # comps <- combinations(
  #   n = length(uniq_idents),
  #   r = 2,
  #   v = uniq_idents
  # ) %>%
  #   as.data.frame(stringsAsFactors = F) %>%
  #   dplyr::bind_rows(data.frame(V1 = uniq_idents, V2 = uniq_idents))

  comps <- expand.grid(
    uniq_idents, uniq_idents,
    stringsAsFactors = F
  )

  if (!is.null(ref_cluster)) {
    comps <- data.frame(
      Var1 = ref_cluster,
      Var2 = uniq_idents,
      stringsAsFactors = F
    )
  }

  # Calculate Jaccard index for comparisons
  res <- purrr::map2_dfr(
    .x = comps$Var1,
    .y = comps$Var2,
    .f = ~ calc_jidx(j_df, c(.x, .y), clonotype_col = clonotype_col)
  )

  res <- dplyr::mutate(res, Var1 = str_c(prefix, Var1, "_jaccard"))
  res <- tidyr::pivot_wider(res, names_from = Var1, values_from = jaccard)

  # Add jaccard index to meta.data
  vdj_meta <- tibble::as_tibble(vdj_meta, rownames = "cell_id")

  vdj_meta <- dplyr::left_join(vdj_meta, res, by = c("idents" = "Var2"))
  vdj_meta <- dplyr::select(vdj_meta, -idents)
  vdj_meta <- tibble::column_to_rownames(vdj_meta, "cell_id")
  vdj_meta <- as.data.frame(vdj_meta)

  res <- Seurat::AddMetaData(sobj_in, vdj_meta)

  res
}


#' Cluster cells based on receptor sequence
#'
#' @param sobj_in Seurat object
#' @param cdr3_col meta.data column containing CDR3 sequences to use for
#' calculating Levenshtein distance
#' @param resolution Clustering resolution to pass to FindClusters
#' @param use_chains Chains to use for calculating Levenshtein distance. If
#' multiple sequences are present for a chain, the first sequence is used.
#' @param prefix Prefix to add to graph name
#' @param ... Additional parameters to pass to FindClusters
#' @return Seurat object with an added shared nearest neighbors graph (vdj_snn)
#' and a meta.data column containing cluster ids
#' @export
cluster_vdj <- function(sobj_in, cdr3_col = "cdr3s_aa", resolution = 0.1,
                        use_chains = NULL, prefix = "vdj_") {

  # Extract sequences
  # Only include cells with VDJ data
  orig_idents <- Idents(sobj_in)

  seqs <- Seurat::FetchData(sobj_in, cdr3_col)
  seqs <- tibble::rownames_to_column(seqs, "cell_id")
  seqs <- na.omit(seqs)

  # Select chains to used for calculating distance
  if (!is.null(use_chains)) {
    re       <- stringr::str_c("(?<=", use_chains, ":)[A-Z]+")
    seq_cols <- stringr::str_c("V", seq_along(use_chains))

    seqs <- purrr::map2(re, seq_cols, ~ {
      n_col <- dplyr::sym(.y)
      c_col <- dplyr::sym(cdr3_col)

      dplyr::mutate(seqs, !!n_col := stringr::str_extract(!!c_col, .x))
    })

    seqs <- purrr::reduce(seqs, dplyr::left_join, by = c("cell_id", cdr3_col))
    seqs <- na.omit(seqs)
    seqs <- dplyr::mutate(seqs, seq_c = stringr::str_c(!!!dplyr::syms(seq_cols)))

  } else {
    seqs <- mutate(
      seqs,
      seq_c = stringr::str_extract_all(!!dplyr::sym(cdr3_col), "(?<=:)[A-Z]+"),
      seq_c = purrr::map(seq_c, purrr::reduce, stringr::str_c)
    )
  }

  # Create Levenshtein distance matrix
  seqs     <- purrr::set_names(seqs$seq_c, seqs$cell_id)
  vdj_dist <- adist(seqs)

  # Create nearest neighbors graph
  # Add graph this way or error thrown due to differing number of cells
  vdj_snn  <- Seurat::FindNeighbors(vdj_dist, distance.matrix = T)
  snn_name <- stringr::str_c(prefix, "snn")

  sobj_in@graphs[[snn_name]] <- vdj_snn$snn

  # Find clusters
  res <- Seurat::FindClusters(
    object     = sobj_in,
    resolution = resolution,
    graph.name = snn_name,
    ...
  )

  Idents(res) <- orig_idents

  res
}


#' Run UMAP using Seurat object containing VDJ nearest neighbors graph
#'
#' @param sobj_in Seurat object containing shared nearest neighbors graph for
#' VDJ data
#' @param umap_key Key to use for UMAP columns in meta.data
#' @param vdj_graph Name of shared nearest neighbors graph stored in Seurat
#' object
#' @return Seurat object containing UMAP coordinates in meta.data
#' @export
run_umap_vdj <- function(sobj_in, umap_key = "vdjUMAP_", vdj_graph = "vdj_snn") {

  # Subset sobj_in to only include VDJ cells and add vdj_snn graph
  # RunUMAP does not like running with a graph that does not include results
  # for all cells in the object
  vdj_cells <- rownames(sobj_in[[vdj_graph]])

  vdj_so <- subset(sobj_in, cells = vdj_cells)
  vdj_so[[vdj_graph]] <- sobj_in[[vdj_graph]]

  # Run UMAP and add reduction object back to original object
  vdj_so <- Seurat::RunUMAP(
    object         = vdj_so,
    reduction.name = "vdj_umap",
    reduction.key  = umap_key,
    graph          = vdj_graph
  )

  umap_coords <- Seurat::Embeddings(vdj_so, reduction = "vdj_umap")
  umap_cols   <- str_c(umap_key, c("1", "2"))

  res <- Seurat::AddMetaData(
    object   = sobj_in,
    metadata = umap_coords,
    col.name = umap_cols
  )

  res
}


#' Subset Seurat object based on VDJ meta.data
#'
#' @param sobj_in Seurat object containing CDR3 sequences
#' @param ... Expression to use for filtering object. To filter based on
#' receptor chains and CDR3 sequences use the terms `.chains` and `.seqs`.
#' @param cdr3_col meta.data column containing CDR3 sequences to use for
#' filtering
#' @return Subsetted Seurat object
#' @export
filter_vdj <- function(sobj_in, ..., cdr3_col = "cdr3s_aa") {

  cdr3_col <- dplyr::sym(cdr3_col)

  # Format meta.data for filtering
  meta_df <- tibble::as_tibble(sobj_in@meta.data, rownames = ".cell_id")

  vdj_df <- dplyr::mutate(
    meta_df,
    .chains = stringr::str_extract_all(!!cdr3_col, "[A-Z]+(?=:)"),
    .seqs   = stringr::str_extract_all(!!cdr3_col, "(?<=:)[A-Z]+")
  )
  vdj_df <- tidyr::unnest(vdj_df, cols = c(.chains, .seqs))
  vdj_df <- dplyr::group_by(vdj_df, .cell_id)
  vdj_df <- dplyr::filter(vdj_df, ...)

  # Subset Seurat object
  vdj_cells <- unique(vdj_df$.cell_id)
  meta_df   <- dplyr::filter(meta_df, .cell_id %in% vdj_cells | is.na(!!cdr3_col))

  res <- subset(sobj_in, cells = meta_df$.cell_id)

  res
}




