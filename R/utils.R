#' Add VDJ data to Seurat object
#'
#' @param sobj_in Seurat object
#' @param vdj_dir cellranger vdj output directory
#' @param productive_pair Only add clonotypes with at least one
#' productive V-J spanning contig for each chain of receptor pair
#' @return Seurat object with VDJ data added to meta.data
import_vdj <- function(sobj_in, vdj_dir, productive_pair = F) {

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

  # Add to meta.data
  res <- Seurat::AddMetaData(sobj_in, metadata = meta_df)

  res
}


#' Calculate receptor diversity (inverse Simpson Index)
#'
#' @param sobj_in Seurat object
#' @param clonotype_col meta.data column containing clonotype ids
#' @param cluster_col meta.data column containing cluster ids used for calculating
#' diversity. If cluster_col is omitted, diversity index will be calculated for
#' all clonotypes
#' @return Seurat object with inverse Simpson index added to meta.data
calc_diversity <- function(sobj_in, clonotype_col = "clonotype_id", cluster_col = NULL) {

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
  vdj_df <- dplyr::summarize(vdj_df, num = dplyr::n_distinct(cell_id))

  # Calculate diversity
  if (!is.null(cluster_col)) {
    vdj_df <- dplyr::group_by(vdj_df, !!cluster_col)
  }

  vdj_df <- dplyr::mutate(
    vdj_df,
    frac      = num / sum(num),
    sum_frac  = sum(frac ^ 2),
    diversity = 1 / sum_frac
    # diversity = sum_frac
  )

  vdj_df <- dplyr::select(vdj_df, all_of(c(vdj_cols, "diversity")))

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
#' @param cell_ident Cell identity to use as a reference for calculating Jaccard index
#' @return Seurat object with Jaccard index added to meta.data
calc_jaccard <- function(sobj_in, clonotype_col = "clonotype_id", cell_ident = NULL) {

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
  so_idents   <- Seurat::Idents(sobj_in)
  uniq_idents <- levels(so_idents)

  ctypes   <- Seurat::FetchData(sobj_in, clonotype_col)
  vdj_meta <- dplyr::bind_cols(ctypes, idents = so_idents)
  vdj_meta <- filter(vdj_meta, !is.na(!!dplyr::sym(clonotype_col)))

  # Create data.frame for calculating Jaccard index
  j_df <- dplyr::mutate(vdj_meta, num = 1)
  j_df <- tidyr::pivot_wider(
    data        = j_df,
    names_from  = idents,
    values_from = num,
    values_fn   = list
  )

  j_df <- tidyr::unnest(j_df, cols = -!!dplyr::sym(clonotype_col))
  j_df <- dplyr::mutate_all(j_df, replace_na, replace = 0)

  # Create data.frame of comparisons
  comps <- expand.grid(
    uniq_idents, uniq_idents,
    stringsAsFactors = F
  )

  if (!is.null(cell_ident)) {
    comps <- data.frame(
      Var1 = cell_ident,
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

  res <- dplyr::mutate(res, Var1 = str_c(Var1, "_jaccard"))
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













# TEST ----

# so <- Read10X("~/Projects/Rincon_scVDJseq/results/KI_DN4_GE/outs/filtered_feature_bc_matrix/") %>%
#   CreateSeuratObject()
#
# so_vdj <- so %>%
#   import_vdj("~/Projects/Rincon_scVDJseq/results/KI_DN4_TCR/outs/", productive_pair = T)
#
# so_vdj@meta.data <- so_vdj@meta.data %>%
#   rownames_to_column("cell_id") %>%
#   mutate(orig.ident = if_else(row_number() < 1000, "grp1", "grp2")) %>%
#   mutate(orig.ident = if_else(row_number() > 2000, "grp3", orig.ident)) %>%
#   column_to_rownames("cell_id")
#
# Idents(so_vdj) <- FetchData(so_vdj, "orig.ident")

# trb_aa <- so_vdj@meta.data %>%
#   as_tibble() %>%
#   mutate(
#     TRA = str_extract(cdr3s_aa, "TRA:[A-Z]+"),
#     TRA = str_remove(TRA, "^TRA:"),
#     TRB = str_extract(cdr3s_aa, "TRB:[A-Z]+"),
#     TRB = str_remove(TRB, "^TRB:")
#   ) %>%
#   pull(TRB) %>%
#   na.omit()

# ctypes <- so_vdj@meta.data %>%
#   select(raw_clonotype_id) %>%
#   na.omit() %>%
#   as_tibble(rownames = "cell_id") %>%
#
#   group_by(raw_clonotype_id) %>%
#
#   summarize(n = n_distinct(cell_id)) %>%
#
#   arrange(desc(n))
#
# ctypes %>%
#   mutate(
#     freq  = n / sum(n),
#     sq    = freq ^ 2,
#     total = sum(sq),
#     D     = 1 / total
#   )
