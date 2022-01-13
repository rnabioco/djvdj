#' Filter V(D)J data
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
#' each cell.
#' @param vdj_cols meta.data columns containing V(D)J data to use for
#' filtering. If NULL, V(D)J data are automatically selected by identifying
#' columns that have NAs in the same rows as clonotype_col.
#' @param clonotype_col meta.data column containing clonotype IDs. This column
#' is used to determine which columns contain V(D)J data.
#' @param sep Separator used for storing per cell V(D)J data
#' @return Object with filtered meta.data
#'
#' @examples
#' # Only include V(D)J data for productive chains
#' filter_vdj(vdj_so, productive)
#'
#' # Only include V(D)J data for cells with paired chains
#' filter_vdj(vdj_sce, paired)
#'
#' # Only include V(D)J data for cells with at least one heavy and one light chain
#' filter_vdj(
#'   vdj_so,
#'   "IGH" %in% chains && any(c("IGK", "IGL") %in% chains)
#' )
#'
#' # Only include V(D)J data for cells that have an IGH, IGK, and IGL chain
#' filter_vdj(
#'   vdj_sce,
#'   all(c("IGH", "IGK", "IGL") %in% chains)
#' )
#'
#' # Remove chains that only have 1 UMI for support
#' filter_vdj(vdj_sce, umis > 1)
#'
#' # Only include V(D)J data for heavy chains
#' filter_vdj(vdj_so, chains == "IGH")
#'
#' @export
filter_vdj <- function(input, filt, vdj_cols = NULL, clonotype_col = "clonotype_id",
                       sep = ";") {

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

  } else {
    length_one <- TRUE

    vdj <- dplyr::mutate(vdj, .KEEP = {{filt}})
  }

  keep_rows <- vdj$.KEEP

  vdj <- dplyr::select(vdj, -.data$.KEEP)

  # If vectors in keep_rows > length 1, filter chains
  if (!length_one) {
    for (clmn in sep_cols) {
      x <- vdj[[clmn]]

      x <- purrr::map2_chr(x, keep_rows, ~ {
        if (!all(is.na(.x))) {
          if (length(.x) != length(.y)) {
            stop(
              "Filtering condition must return TRUE/FALSE for each chain, ",
              "or a single TRUE/FALSE for each cell."
            )
          }
        }

        ifelse(
          all(is.na(.x)) || purrr::is_empty(.x[.y]),
          as.character(NA),
          paste0(.x[.y], collapse = sep)
        )
      })

      vdj <- dplyr::mutate(vdj, !!sym(clmn) := x)
    }

    # Set NAs for vdj_cols not in sep_cols
    other_cols <- vdj_cols[!vdj_cols %in% sep_cols]

    vdj <- dplyr::mutate(
      vdj,
      across(all_of(other_cols), ~ {
        ifelse(
          if_all(all_of(sep_cols), is.na),
          NA,
          .x
        )
      })
    )

  # If vectors in keep_rows == length 1, filter V(D)J data per cell
  } else {
    keep_rows <- unlist(keep_rows)

    for (clmn in vdj_cols) {
      x <- vdj[[clmn]]

      # Collapse sep_cols
      if (!purrr::is_empty(sep_cols) && clmn %in% sep_cols) {
        x <- purrr::map2_chr(x, keep_rows, ~ {
          ifelse(
            all(is.na(.x)) | !.y,
            as.character(NA),
            paste0(.x, collapse = sep)
          )
        })

      } else {
        x <- ifelse(keep_rows, x, NA)
      }

      vdj <- dplyr::mutate(vdj, !!sym(clmn) := x)
    }
  }

  # Format results
  res <- .add_meta(input, vdj)

  res
}

