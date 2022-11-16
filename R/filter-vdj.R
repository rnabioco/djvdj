#' Filter V(D)J data in object
#'
#' Remove V(D)J data for chains/cells that do not satisfy the provided
#' condition
#'
#' @param input Single cell object or data.frame containing V(D)J data. If a
#' data.frame is provided, cell barcodes should be stored as row names.
#' @param filt Condition to use for filtering V(D)J data. To allow for
#' filtering of per-chain V(D)J data, the data for each cell is converted into
#' a vector, e.g. 'IGH;IGK' is equivalent to c('IGH', 'IGK'). This allows R
#' vector operations to be performed on the per-chain values. The filtering
#' condition must return TRUE/FALSE for each chain or a single TRUE/FALSE for
#' each cell. Data can be filtered based on cell barcodes by referring to the
#' '.cell_id' column.
#' @param data_cols meta.data columns containing V(D)J data to use for
#' filtering. If NULL, V(D)J data are automatically selected by identifying
#' columns that have NAs in the same rows as clonotype_col.
#' @param clonotype_col meta.data column containing clonotype IDs. This column
#' is used to determine which columns contain V(D)J data.
#' @param per_cell If TRUE, per-chain data will not be parsed and values
#' in each meta.data column will be filtered as is.
#' @param sep Separator used for storing per cell V(D)J data. If NULL, per_cell
#' will be set TRUE.
#' @return Object with filtered meta.data
#'
#' @examples
#' # Only include V(D)J data for productive chains
#' res <- filter_vdj(vdj_so, productive)
#'
#' head(res@meta.data, 1)
#'
#' # Only include V(D)J data for cells with paired chains
#' res <- filter_vdj(vdj_sce, paired)
#'
#' head(res@colData, 3)
#'
#' # Only include V(D)J data for cells with at least one heavy and one light chain
#' res <- filter_vdj(
#'   vdj_so,
#'   "IGH" %in% chains && any(c("IGK", "IGL") %in% chains)
#' )
#'
#' head(res@meta.data, 3)
#'
#' # Only include V(D)J data for cells that have an IGH, IGK, and IGL chain
#' res <- filter_vdj(
#'   vdj_sce,
#'   all(c("IGH", "IGK", "IGL") %in% chains)
#' )
#'
#' head(res@colData, 1)
#'
#' # Only include V(D)J data for heavy chains
#' res <- filter_vdj(vdj_so, chains == "IGH")
#'
#' head(res@meta.data, 3)
#'
#' # Remove chains that only have 1 UMI for support
#' res <- filter_vdj(vdj_sce, umis > 1)
#'
#' head(res@colData, 1)
#'
#' # Filter based on cell barcode
#' res <- filter_vdj(vdj_so, .cell_id == "1_AAACCTGAGTAGGTGC-1")
#'
#' head(res@meta.data, 3)
#'
#' @export
filter_vdj <- function(input, filt, data_cols = NULL,
                       clonotype_col = "clonotype_id", sep = ";",
                       per_cell = FALSE) {

  if (per_cell) sep <- NULL

  # Format input data
  meta <- vdj <- .get_meta(input)

  # Identify columns with V(D)J data
  col_list <- .get_vdj_cols(
    df_in     = vdj,
    clone_col = clonotype_col,
    cols_in   = data_cols,
    sep       = sep
  )

  vdj_cols <- col_list$vdj
  sep_cols <- col_list$sep

  if (!is.null(sep) && purrr::is_empty(sep_cols)) {
    warning("The separator '", sep, "' is not present in the data")
  }

  # Create list-cols for V(D)J columns that contain sep
  if (purrr::is_empty(sep_cols)) {
    vdj <- dplyr::mutate(vdj, .KEEP = {{filt}})
    length_one <- TRUE

  } else {
    vdj <- .unnest_vdj(
      df_in    = vdj,
      sep      = sep,
      sep_cols = sep_cols,
      unnest   = FALSE
    )

    vdj <- dplyr::rowwise(vdj)
    vdj <- dplyr::mutate(vdj, .KEEP = list({{filt}}))
    vdj <- dplyr::ungroup(vdj)

    length_one <- map_lgl(vdj$.KEEP, ~ length(.x) == 1)
    length_one <- all(length_one)
  }

  keep_rows <- vdj$.KEEP
  vdj       <- dplyr::select(vdj, -.data$.KEEP)

  # If vectors in keep_rows are all length 1, filter cells
  if (length_one) {
    keep_rows <- unlist(keep_rows)
    vdj_cols  <- vdj_cols[vdj_cols != CELL_COL]

    vdj <- mutate(meta, across(all_of(vdj_cols), .add_na, !keep_rows))

    res <- .add_meta(input, vdj)

    return(res)
  }

  # If vectors in keep_rows > length 1, filter chains
  for (clmn in sep_cols) {
    x <- vdj[[clmn]]

    x <- purrr::map2_chr(x, keep_rows, ~ {
      if (!all(is.na(.x)) && length(.x) != length(.y)) {
        stop(
          "Filtering condition must return TRUE/FALSE for each chain, ",
          "or a single TRUE/FALSE for each cell."
        )
      }

      ifelse(
        all(is.na(.x)) || purrr::is_empty(.x[.y]),
        as.character(NA),
        paste0(.x[.y], collapse = sep)
      )
    })

    vdj[clmn] <- x
  }

  # Set NAs for vdj_cols not in sep_cols
  other_cols <- vdj_cols[!vdj_cols %in% sep_cols]
  na_rows    <- purrr::map_lgl(keep_rows, ~ !any(.x))

  vdj <- mutate(vdj, across(all_of(other_cols), .add_na, na_rows))

  # Format results
  res <- .add_meta(input, vdj)

  res
}

#' Insert NAs based on logical index
#' @noRd
.add_na <- function(x, lgl_idx) {
  typ <- typeof(x)

  x[lgl_idx] <- as(NA, typ)

  x
}




# vdj <- dplyr::mutate(
#   vdj,
#   across(all_of(other_cols), ~ {
#     ifelse(
#       if_all(all_of(sep_cols), is.na),
#       NA,
#       .x
#     )
#   })
# )

# for (clmn in vdj_cols) {
#   x <- vdj[[clmn]]
#
#   # Collapse sep_cols
#   if (!purrr::is_empty(sep_cols) && clmn %in% sep_cols) {
#     x <- purrr::map2_chr(x, keep_rows, ~ {
#       ifelse(
#         all(is.na(.x)) | !.y,
#         as.character(NA),
#         paste0(.x, collapse = sep)
#       )
#     })
#
#   } else x <- ifelse(keep_rows, x, NA)
#
#   vdj <- dplyr::mutate(vdj, !!sym(clmn) := x)
# }
