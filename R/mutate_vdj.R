#' Mutate V(D)J data
#'
#' @param input Object containing V(D)J data. If a data.frame is provided, the
#' cell barcodes should be stored as row names.
#' @param ... Name-value pairs to use for creating or modifying meta.data
#' columns
#' @param clonotype_col meta.data column containing clonotype IDs. This column
#' is used to determine which cells have V(D)J data
#' @param vdj_cols meta.data columns containing VDJ data to use for filtering.
#' If set to NULL (recommended) columns are automatically selected.
#' @param return_df Return results as a data.frame. If set to FALSE, results
#' will be added to the input object.
#' @param sep Separator used for storing per cell V(D)J data
#' @return Object with mutated meta.data
#' @export
mutate_vdj <- function(input, ..., clonotype_col = "cdr3_nt", vdj_cols = NULL, return_df = FALSE, sep = ";") {

  # Format input data
  meta <- .get_meta(input)
  vdj  <- meta

  # Identify columns with VDJ data
  col_list <- .get_vdj_cols(
    df_in     = vdj,
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

    vdj <- .split_vdj(
      df_in    = vdj,
      sep_cols = sep_cols,
      sep      = sep
    )

    vdj <- dplyr::rowwise(vdj)
  }

  # Mutate meta.data
  vdj <- dplyr::mutate(vdj, ...)
  vdj <- dplyr::ungroup(vdj)

  # Remove columns created for mutate
  if (!purrr::is_empty(sep_cols)) {
    sep_cols <- purrr::set_names(
      x  = names(sep_cols),
      nm = unname(sep_cols)
    )

    vdj <- dplyr::select(vdj, !all_of(names(sep_cols)))
    vdj <- dplyr::rename(vdj, !!!syms(sep_cols))
  }

  # Format results
  in_cols <- colnames(meta)

  res <- dplyr::relocate(vdj, all_of(in_cols), .after = .data$.cell_id)

  if (return_df) {
    input <- meta
  }

  res <- .add_meta(input, meta = res)

  res
}


#' Mutate object meta.data
#'
#' @param input Object containing single cell data
#' @param .fun Function or formula to use for modifying the meta.data. If a
#' formula is provided, use .x to refer to the meta.data table.
#' @param ... Arguments to pass to the provided function
#' @return Object with mutated meta.data
#' @export
mutate_meta <- function(input, .fun, ...) {

  if (!is_function(.fun) && !is_formula(.fun)) {
    stop(".fun must be either a function or a formula")
  }

  meta <- .get_meta(input)

  if (is_formula(.fun)) {
    .fun <- as_mapper(.fun, ...)
  }

  res <- .fun(meta, ...)

  res <- .add_meta(input, meta = res)

  res
}


