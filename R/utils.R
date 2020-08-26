#' Add VDJ data to Seurat object
#'
#' @param sobj_in Seurat object
#' @param vdj_dir cellranger vdj output directory
#' @param productive_pair Only add clonotypes with at least one productive V-J
#' spanning contig for each chain of receptor pair
#' @param prefix Prefix to add to new meta.data columns
#' @return Seurat object with VDJ data added to meta.data
import_vdj <- function(sobj_in, vdj_dir, productive_pair = F, prefix = "") {

  # Filtered contigs
  ctypes <- readr::read_csv(file.path(vdj_dir, "filtered_contig_annotations.csv"))

  # Only include clonotypes with productive V-J spanning contigs
  if (productive_pair) {
    TCR_chains <- c("TRA", "TRB")
    BCR_chains <- c("IGH", "IGL", "IGK")

    if (any(TCR_chains %in% ctypes$chain)) {
      all_chains <- TCR_chains

    } else {
      all_chains <- BCR_chains
    }

    ctypes <- dplyr::filter(ctypes, productive)
    ctypes <- dplyr::group_by(ctypes, raw_clonotype_id)
    ctypes <- dplyr::filter(ctypes, all(all_chains %in% chain))  # FIX THIS FOR BCRs
    ctypes <- dplyr::ungroup(ctypes)
  }

  ctypes <- dplyr::select(ctypes, barcode, raw_clonotype_id)
  ctypes <- unique(ctypes)

  # Retrieve CDR3 sequences for clonotypes passing filters
  consensus_ctypes <- readr::read_csv(file.path(vdj_dir, "clonotypes.csv"))

  ctypes <- dplyr::left_join(
    ctypes, consensus_ctypes,
    by = c("raw_clonotype_id" = "clonotype_id")
  )

  ctypes <- dplyr::rename(ctypes, clonotype_id = raw_clonotype_id)
  ctypes <- dplyr::select(ctypes, -frequency, -proportion)

  # Create data.frame with clonotype info
  cells   <- data.frame(cell_id = Cells(sobj_in))
  meta_df <- dplyr::left_join(cells, ctypes, by = c("cell_id" = "barcode"))
  meta_df <- tibble::column_to_rownames(meta_df, "cell_id")
  meta_df <- dplyr::rename_all(meta_df, ~ str_c(prefix, .))

  # Add to meta.data
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
    !!dplyr::sym(div_col) := 1 / sum_frac
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
#' @param prefix Prefix to add to new meta.data columns
#' @return Seurat object with Jaccard index added to meta.data
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

    a     <- length(row_sums[row_sums == 2])      # Similarity
    a_b_c <- length(row_sums[row_sums == 1]) + a  # Total

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
#' @param ... Additional parameters to pass to FindClusters
#' @return Seurat object with an added shared nearest neighbor graph (vdj_snn)
#' and a meta.data column containing cluster ids
cluster_vdj <- function(sobj_in, cdr3_col = "cdr3s_aa", resolution = 0.1, ...) {

  seqs <- Seurat::FetchData(sobj_in, cdr3_col)
  seqs <- purrr::set_names(pull(seqs, cdr3_col), rownames(seqs))
  seqs <- purrr::map_chr(seqs, stringr::str_remove_all, ";*(IGH|IGK|IGL):")

  vdj_dist <- adist(seqs)
  vdj_snn  <- Seurat::FindNeighbors(vdj_dist, distance.matrix = T)

  sobj_in[["vdj_snn"]] <- vdj_snn$snn

  res <- Seurat::FindClusters(
    object     = sobj_in,
    resolution = resolution,
    graph.name = "vdj_snn",
    ...
  )

  res
}












# STANDARD WORKFLOW ----

# # Load data
# library(tidyverse)
# library(Seurat)
#
# # data_dir <- "~/Projects/Rincon_scVDJseq/results/KI_DN4_GE/outs"
# data_dir    <- "~/Projects/Smith_AVIDseq/2020-07-17"
# so_list     <- Read10X(file.path(data_dir, "JH191_GEX/outs/filtered_feature_bc_matrix"))
# so          <- CreateSeuratObject(so_list$`Gene Expression`)
# so[["ADT"]] <- CreateAssayObject(so_list$`Antibody Capture`)
#
# # QC filtering
# so <- so %>%
#   PercentageFeatureSet(
#     pattern  = "^mt-",
#     col.name = "Percent_mito"
#   ) %>%
#   subset(
#     nFeature_RNA > 200 &
#     nFeature_RNA < 5000 &
#     Percent_mito < 15
#   )
#
# # Normalize
# so <- so %>%
#   NormalizeData() %>%
#   NormalizeData(
#     assay = "ADT",
#     normalization.method = "CLR"
#   ) %>%
#   FindVariableFeatures() %>%
#   ScaleData()
#
# # Cluster
# so <- so %>%
#   RunPCA() %>%
#   RunUMAP(dims = 1:40) %>%
#   FindNeighbors(dims = 1:40) %>%
#   FindClusters(resolution = 0.2)


# VDJ WORKFLOW ----

# # Add VDJ data to meta.data
# so_vdj <- import_vdj(
#   sobj_in         = so,
#   vdj_dir         = file.path(data_dir, "BCR/outs"),
#   productive_pair = F,
#   prefix          = ""
# )
#
# # Calculate repertoire diversity
# so_vdj <- calc_diversity(
#   sobj_in       = so_vdj,
#   clonotype_col = "clonotype_id",
#   cluster_col   = "seurat_clusters",
#   prefix        = ""
# )
#
# # Calculate repertoire overlap
# so_vdj <- calc_jaccard(
#   sobj_in       = so_vdj,
#   clonotype_col = "clonotype_id",
#   cluster_col   = "seurat_clusters",
#   ref_cluster   = "5",
#   prefix        = "x"
# )
#
# # Cluster based on receptor sequence
# so_vdj <- cluster_vdj(
#   sobj_in    = so_vdj,
#   cdr3_col   = "cdr3s_aa",
#   resolution = 0.1
# )
#
# so_vdj <- so_vdj %>%
#   RunUMAP(
#     reduction.name = "vdj_umap",
#     reduction.key  = "vdjUMAP_",
#     graph          = "vdj_snn",
#     verbose        = F
#   )

