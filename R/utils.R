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
  )

  vdj_df <- dplyr::select(vdj_df, all_of(c(vdj_cols, "diversity")))

  # Add resuts to meta.data
  meta_df <- dplyr::left_join(meta_df, vdj_df, by = vdj_cols)
  meta_df <- tibble::column_to_rownames(meta_df, "cell_id")
  meta_df <- as.data.frame(meta_df)

  res <- Seurat::AddMetaData(sobj_in, metadata = meta_df)

  res
}





# TEST ----

# so <- Read10X("~/Projects/Rincon_scVDJseq/results/KI_DN4_GE/outs/filtered_feature_bc_matrix/") %>%
#   CreateSeuratObject()
#
# so_vdj <- so %>%
#   import_vdj("~/Projects/Rincon_scVDJseq/results/KI_DN4_TCR/outs/", productive_vj_pair = T)

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
