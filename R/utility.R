#' Add VDJ data to Seurat object
#'
#' @param sobj_in Seurat object
#' @param vdj_dir cellranger vdj output directory
#' @param productive_vj_pair Only include clonotypes with at least one
#' productive V-J spanning contig for each chain of receptor pair
#' @return Seurat object with VDJ data added to meta.data
import_vdj <- function(sobj_in, vdj_dir, productive_vj_pair = F) {

  # Filtered contigs
  ctypes <- readr::read_csv(file.path(vdj_dir, "filtered_contig_annotations.csv"))

  # Only include clonotypes with productive V-J spanning contigs
  if (productive_vj_pair) {
    TCR_chains <- c("TRA", "TRB")
    BCR_chains <- c("IGH", "IGL", "IGK")

    if (any(TCR_chains %in% ctypes$chain)) {
      all_chains <- TCR_chains

    } else {
      all_chains <- BCR_chains
    }

    ctypes <- dplyr::filter(ctypes, productive)
    ctypes <- dplyr::group_by(ctypes, raw_clonotype_id)
    ctypes <- dplyr::filter(ctypes, all(all_chains %in% chain))
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

  # Create data.frame with clonotype info
  cells   <- data.frame(cell_id = Cells(sobj_in))
  meta_df <- left_join(cells, ctypes, by = c("cell_id" = "barcode"))
  meta_df <- column_to_rownames(meta_df, "cell_id")

  # Add to meta.data
  res <- AddMetaData(sobj_in, metadata = meta_df)

  res
}


