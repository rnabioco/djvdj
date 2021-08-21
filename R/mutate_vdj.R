#' Mutate V(D)J data
#'
#'@export
mutate_vdj <- function(input, ...) {
  UseMethod("mutate_vdj", input)
}

#' @rdname mutate_vdj
#' @param input Single cell object or data.frame containing V(D)J data. If a
#' data.frame is provided, the cell barcodes should be stored as row names.
#' @param ... Name-value pairs to use for creating or modifying meta.data
#' columns
#' @param clonotype_col meta.data column containing clonotype IDs. This column
#' is used to determine which cells have V(D)J data
#' @param sep Separator to use for expanding meta.data columns
#' @param vdj_cols meta.data columns containing VDJ data to use for filtering.
#' If set to NULL (recommended) columns are automatically selected.
#' @param return_df Return results as a data.frame. If set to FALSE, results
#' will be added to the input object.
#' @return Modified single cell object or data.frame
#' @export
mutate_vdj.default <- function(input, ..., clonotype_col = "cdr3_nt", sep = ";", vdj_cols = NULL) {

  # Format input data
  meta_df <- tibble::as_tibble(input, rownames = ".cell_id")

  # Identify columns with VDJ data
  col_list <- .get_vdj_cols(
    df_in     = meta_df,
    clone_col = clonotype_col,
    cols_in   = vdj_cols,
    sep       = sep
  )

  vdj_cols <- col_list$vdj
  sep_cols <- col_list$sep

  # Create list-cols for VDJ columns that contain sep
  if (!purrr::is_empty(sep_cols)) {
    sep_cols <- set_names(
      x  = sep_cols,
      nm = paste0(".", sep_cols)
    )

    meta_df <- .split_vdj(
      df_in    = meta_df,
      sep_cols = sep_cols,
      sep      = sep
    )

    meta_df <- dplyr::rowwise(meta_df)
  }

  # Mutate meta.data
  meta_df <- dplyr::mutate(meta_df, ...)
  meta_df <- dplyr::ungroup(meta_df)

  # Remove columns created for mutate
  if (!purrr::is_empty(sep_cols)) {
    sep_cols <- purrr::set_names(
      x  = names(sep_cols),
      nm = unname(sep_cols)
    )

    meta_df <- dplyr::select(
      meta_df,
      !all_of(names(sep_cols))
    )

    meta_df <- dplyr::rename(meta_df, !!!syms(sep_cols))
  }

  # Format results
  in_cols <- colnames(input)

  res <- dplyr::relocate(meta_df, all_of(in_cols), .after = .cell_id)
  res <- tibble::column_to_rownames(res, ".cell_id")

  res
}

#' @rdname mutate_vdj
#' @export
mutate_vdj.Seurat <- function(input, ..., clonotype_col = "cdr3_nt", sep = ";", vdj_cols = NULL,
                              return_df = FALSE) {

  res <- mutate_vdj(
    input         = input@meta.data,
    ...,
    clonotype_col = clonotype_col,
    sep           = sep,
    vdj_cols      = vdj_cols
  )

  if (!return_df) {
    res <- Seurat::AddMetaData(input, metadata = res)
  }

  res
}
