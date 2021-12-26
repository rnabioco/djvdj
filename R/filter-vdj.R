#' Filter V(D)J data
#'
#' @param input Single cell object or data.frame containing V(D)J data. If a
#' data.frame is provided, the cell barcodes should be stored as row names.
#' @param filt Condition to use for filtering meta.data
#' @param clonotype_col meta.data column containing clonotype IDs. This column
#' is used to determine which cells have V(D)J data. If clonotype_col is set to
#' NULL, filtering is performed regardless of whether V(D)J data is present for
#' the cell
#' @param filter_cells Should cells be removed from the object? If FALSE
#' (default) V(D)J data will be removed from the meta.data but no cells will be
#' removed from the object.
#' @param vdj_cols meta.data columns containing V(D)J data to use for filtering.
#' If set to NULL (recommended) columns are automatically selected by
#' identifying columns that have NAs in the same rows as clonotype_col.
#' @param sep Separator used for storing per cell V(D)J data
#' @return Object with filtered meta.data
#' @export
filter_vdj <- function(input, filt, clonotype_col = "clonotype_id", filter_cells = FALSE,
                       vdj_cols = NULL, sep = ";") {

  if (!filter_cells && is.null(clonotype_col)) {
    stop("clonotype_col must be provided when filter_cells is set to FALSE")
  }

  # Format input data
  meta <- .get_meta(input)
  vdj  <- meta

  # Identify columns with V(D)J data
  col_list <- .get_vdj_cols(
    df_in     = vdj,
    clone_col = clonotype_col,
    cols_in   = vdj_cols,
    sep       = sep
  )

  vdj_cols <- col_list$vdj
  sep_cols <- col_list$sep

  if (purrr::is_empty(sep_cols)) {
    warning("The separator '", sep, "' is not present in the data")
  }

  # Create list-cols for V(D)J columns that contain sep
  if (!purrr::is_empty(sep_cols)) {
    sep_cols <- set_names(
      x  = sep_cols,
      nm = paste0(".", sep_cols)
    )

    vdj <- .unnest_vdj(
      df_in    = vdj,
      sep      = sep,
      sep_cols = sep_cols,
      unnest   = FALSE
    )

    vdj <- dplyr::rowwise(vdj)
  }

  # Store results for filtering
  vdj <- dplyr::mutate(vdj, .KEEP = {{filt}})
  vdj <- dplyr::ungroup(vdj)

  # Remove V(D)J data from meta.data (without filtering cells)
  if (!filter_cells) {
    vdj <- dplyr::mutate(
      vdj,
      across(
        all_of(c(vdj_cols, names(sep_cols))),
        ~ ifelse(.KEEP, .x, NA)
      )
    )

  # Filter cells from object
  # when clonotype_col != NULL only filter cells with V(D)J data
  } else {
    if (!is.null(clonotype_col)) {
      vdj <- dplyr::mutate(
        vdj,
        .KEEP = dplyr::if_else(
          is.na(!!sym(clonotype_col)),
          TRUE,
          .data$.KEEP
        )
      )
    }

    vdj <- dplyr::filter(vdj, .data$.KEEP)
  }

  # Remove columns created for filtering
  if (!purrr::is_empty(sep_cols)) {
    sep_cols <- purrr::set_names(
      x  = names(sep_cols),
      nm = unname(sep_cols)
    )

    vdj <- dplyr::select(vdj, -all_of(names(sep_cols)))
    vdj <- dplyr::rename(vdj, !!!syms(sep_cols))
  }

  vdj <- dplyr::select(vdj, -.data$.KEEP)

  # Format results
  in_cols <- colnames(meta)

  vdj <- dplyr::relocate(vdj, all_of(in_cols), .after = .data$.cell_id)

  cells <- vdj$.cell_id
  res   <- .subset_cells(input, cells)
  res   <- .add_meta(res, meta = vdj)

  res
}

#' Subset object so it only contains given cells
#'
#' @param input Single cell object or data.frame. If a data.frame is provided, the
#' cell barcodes should be stored as row names.
#' @param cells Vector with cell barcodes to use for subsetting
#' @return Subsetted object
#' @noRd
.subset_cells <- function(input, cells) {
  UseMethod(".subset_cells", input)
}

.subset_cells.default <- function(input, cells) {
  res <- input[, cells]

  res
}

.subset_cells.data.frame <- function(input, cells) {
  res <- input[rownames(input) %in% cells, ]

  res
}


