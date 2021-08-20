#' Calculate clonotype abundance
#'
#' @export
calc_abundance <- function(input, ...) {
  UseMethod("calc_abundance", input)
}

#' @rdname calc_abundance
#' @param input Single cell object or data.frame containing V(D)J data. If a
#' data.frame is provided, the cell barcodes should be stored as row names.
#' @param clonotype_col meta.data column containing clonotype IDs
#' @param cluster_col meta.data column containing cluster IDs to use for
#' grouping cells when calculating clonotype abundance
#' @param prefix Prefix to add to new meta.data columns
#' @param return_seurat Return a Seurat object. If set to FALSE, a tibble
#' summarizing the results is returned.
#' @return Seurat object with clonotype abundance added to meta.data
#' @export
calc_abundance.default <- function(input, clonotype_col = "cdr3_nt", cluster_col = NULL,
                                   prefix = "", ...) {

  # Format meta.data
  input   <- tibble::as_tibble(input, rownames = ".cell_id")
  meta_df <- dplyr::filter(input, !is.na(!!sym(clonotype_col)))

  meta_df <- dplyr::select(
    meta_df,
    .data$.cell_id, all_of(c(cluster_col, clonotype_col))
  )

  # Calculate clonotype abundance
  meta_df <- .calc_abund(
    df_in     = meta_df,
    cell_col  = ".cell_id",
    clone_col = clonotype_col,
    clust_col = cluster_col
  )

  new_cols <- c("freq", "pct")

  if (!is.null(cluster_col)) {
    new_cols <- c(new_cols, "shared")
  }

  new_cols <- purrr::set_names(
    paste0(".", new_cols),
    paste0(prefix, "clone_", new_cols)
  )

  meta_df <- select(meta_df, .cell_id, !!!syms(new_cols))

  res <- dplyr::left_join(input, meta_df, by = ".cell_id")
  res <- tibble::column_to_rownames(res, ".cell_id")

  res
}

#' @rdname calc_abundance
#' @export
calc_abundance.Seurat <- function(input, clonotype_col = "cdr3_nt", cluster_col = NULL,
                                  prefix = "", return_seurat = TRUE) {

  res <- calc_abundance(
    input         = input@meta.data,
    clonotype_col = clonotype_col,
    cluster_col   = cluster_col,
    prefix        = prefix
  )

  if (return_seurat) {
    res <- Seurat::AddMetaData(input, metadata = res)
  }

  res
}

#' Calculate clonotype abundance
#'
#' @param df_in Input data.frame
#' @param cell_col Column containing cell IDs
#' @param clone_col Column containing clonotype IDs
#' @param clust_col Column containing cluster IDs
#' @return data.frame
.calc_abund <- function(df_in, cell_col, clone_col, clust_col = NULL) {

  # Count number of cells in each group
  if (!is.null(clust_col)) {
    df_in <- dplyr::group_by(df_in, !!sym(clust_col))
  }

  df_in <- dplyr::mutate(
    df_in,
    .n_cells = dplyr::n_distinct(!!sym(cell_col))
  )

  # Calculate frequency
  res <- dplyr::group_by(df_in, !!sym(clone_col), .add = TRUE)

  res <- dplyr::mutate(
    res,
    .freq = dplyr::n_distinct(!!sym(cell_col)),
    .pct  = (.data$.freq / .data$.n_cells) * 100
  )

  # Identify shared clonotypes
  if (!is.null(clust_col)) {
    res <- dplyr::group_by(res, !!sym(clone_col))

    res <- dplyr::mutate(
      res,
      .shared = dplyr::n_distinct(!!sym(clust_col)) > 1
    )
  }

  res <- dplyr::ungroup(res)

  res
}
