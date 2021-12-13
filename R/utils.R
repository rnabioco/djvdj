#' Helper to test all combinations of provided arguments
#'
#' @param arg_lst Named list of arguments to test
#' @param .fn Function to test
#' @param desc Description to pass to test_that
#' @param chk Function or expression to using for testing
#' @param dryrun Do not run tests, just return table of arguments that will be
#' tested
#' @return Output from test_that
#' @noRd
test_all_args <- function(arg_lst, .fn, desc, chk, dryrun = FALSE) {

  arg_lst <- expand.grid(arg_lst, stringsAsFactors = FALSE)

  if (dryrun) {
    return(arg_lst)
  }

  n <- 1

  pwalk(arg_lst, ~ {
    test_that(paste(desc, n), {

      if (is.call(chk)) {
        .res <- .fn(...)

        return(eval(chk))
      }

      chk(.fn(...))
    })

    n <<- n + 1
  })
}


#' Summarize values for chains
#'
#' Summarize values present for each column provided to the data_cols argument.
#' For each cell, the function(s) provided will be applied to each unique
#' chain.
#'
#' @param input Single cell object or data.frame containing V(D)J data. If a
#' data.frame is provided, the cell barcodes should be stored as row names.
#' @param data_cols meta.data columns to summarize
#' @param fn Function to use for summarizing data_cols
#' @param chain_col meta.data column(s) containing labels for each chain
#' expressed in the cell. These labels are used for grouping the summary
#' output. Set chain_col to NULL to group solely based on the cell barcodes.
#' @param include_cols Additional columns to include in the output data.frame
#' @param sep Separator used for storing per cell V(D)J data
#' @return data.frame containing summary results
#' @export
summarize_chains <- function(input, data_cols = c("umis", "reads"), fn, chain_col = "chains",
                             include_cols = NULL, sep = ";") {

  # Format input data
  fetch_cols <- c(".cell_id", data_cols, chain_col, include_cols)

  meta <- .get_meta(input)
  meta <- dplyr::select(meta, all_of(fetch_cols))

  meta <- dplyr::filter(
    meta,
    across(all_of(data_cols), ~ !is.na(.x))
  )

  # Expand meta.data
  coerce_cols <- purrr::set_names(
    rep("numeric", length(data_cols)),
    data_cols
  )

  res <- .split_vdj(
    meta,
    sep         = sep,
    sep_cols    = c(data_cols, chain_col),
    coerce_cols = coerce_cols,
    expand      = TRUE
  )

  # Summarize chains
  grp_cols <- c(".cell_id", chain_col, include_cols)
  res      <- dplyr::group_by(res, !!!syms(grp_cols))

  res <- dplyr::summarize(
    res,
    across(all_of(data_cols), fn),
    .groups = "drop"
  )

  res
}


#' Add meta.data to single cell object
#'
#' @param input Object containing single cell data
#' @param meta meta.data to add to object
#' @param row_col Column containing meta.data rownames
#' @return Object with added meta.data
#' @noRd
.add_meta <- function(input, meta, row_col) {

  UseMethod(".add_meta", input)

}

.add_meta.default <- function(input, meta, row_col = ".cell_id") {

  tibble::column_to_rownames(meta, row_col)

}

.add_meta.Seurat <- function(input, meta, row_col = ".cell_id") {

  meta <- tibble::column_to_rownames(meta, row_col)

  input@meta.data <- meta

  input
}

.add_meta.SingleCellExperiment <- function(input, meta, row_col = ".cell_id") {

  meta <- tibble::column_to_rownames(meta, row_col)

  input@colData <- S4Vectors::DataFrame(meta)

  input
}


#' Pull meta.data from single cell object
#'
#' @param input Object containing single cell data
#' @param row_col New column to store meta.data rownames
#' @return tibble containing meta.data pulled from object
#' @noRd
.get_meta <- function(input, row_col) {

  UseMethod(".get_meta", input)

}

.get_meta.default <- function(input, row_col = ".cell_id") {

  .to_tibble(input, row_col)

}

.get_meta.Seurat <- function(input, row_col = ".cell_id") {

  .to_tibble(input@meta.data, row_col)

}

.get_meta.SingleCellExperiment <- function(input, row_col = ".cell_id") {

  .to_tibble(input@colData, row_col)

}


#' Convert to tibble
#'
#' The point of this function is to handle situations where the name passed to
#' rownames already exists in the data.frame. This occurs when running
#' .get_meta() multiple times on the same object. The default behavior of
#' as_tibble is to throw an error. This function just returns the tibble
#' without adding a new rowname column.
#'
#' @param input Object to coerce to tibble
#' @param row_col Name of new column to store rownames
#' @return tibble with rownames added to new column
#' @noRd
.to_tibble <- function(input, row_col) {

  res <- tibble::as_tibble(input, rownames = NA)

  if (!is.null(row_col) && !row_col %in% colnames(res)) {
    res <- tibble::rownames_to_column(res, row_col)
  }

  res
}


#' Merge new meta.data with object
#'
#' @param input Object containing single cell data
#' @param meta meta.data to merge with object
#' @param by Columns to use for merging
#' @return Object with added meta.data
#' @noRd
.merge_meta <- function(input, meta, by = ".cell_id") {

  if (is.null(input)) {
    return(meta)
  }

  # Join meta.data
  # remove columns already present in input object to prevent duplicates
  # this mimics behavior of Seurat::AddMetaData
  meta     <- .get_meta(meta, row_col = ".cell_id")
  rm_cols  <- colnames(meta)
  rm_cols  <- rm_cols[!rm_cols %in% by]
  obj_meta <- .get_meta(input)
  obj_meta <- dplyr::select(obj_meta, !any_of(rm_cols))

  meta <- dplyr::left_join(obj_meta, meta, by = by)
  res  <- .add_meta(input, meta = meta)

  res
}


#' Split V(D)J meta.data columns into vectors
#'
#' @param df_in data.frame to modify
#' @param sep Separator used for storing per cell V(D)J data
#' @param sep_cols Columns to split based on sep
#' @param expand Should columns be unnested after splitting into vectors
#' @param coerce_cols Named vector specifying columns that should be coerced to
#' a new type. The vector names should be the name of each column,
#' e.g., c(umis = "numeric")
#' @return data.frame
#' @noRd
.split_vdj <- function(df_in, sep = ";", sep_cols, expand = FALSE,
                       coerce_cols = c(
                         umis        = "numeric", reads          = "numeric",
                         cdr3_length = "numeric", cdr3_nt_length = "numeric",
                         productive  = "logical", full_length    = "logical"
                       )) {

  # Add new set of columns based on sep_cols names
  # this is useful if want to leave original columns unmodified
  res <- dplyr::mutate(df_in, !!!syms(sep_cols))

  # Split columns into vectors
  res <- dplyr::mutate(res, across(
    all_of(unname(sep_cols)),
    ~ strsplit(as.character(.x), sep)
  ))

  # Coerce columns to correct types
  if (!is.null(coerce_cols)) {
    coerce_cols <- coerce_cols[names(coerce_cols) %in% sep_cols]

    purrr::iwalk(coerce_cols, ~ {
      Class <- .x
      clmn <- sym(.y)

      res <<- dplyr::mutate(res, !!clmn := map(!!clmn, .convert_char, Class))
    })
  }

  # Unnest columns
  if (expand) {
    res <- tidyr::unnest(res, all_of(unname(sep_cols)))
  }

  res
}


#' Identify columns with V(D)J data
#'
#' @param df_in data.frame
#' @param clone_col Column containing clonotype IDs to use for identifying
#' columns with V(D)J data. If both clone_col and cols_in are set to NULL all
#' columns are used.
#' @param cols_in meta.data columns containing V(D)J data to use for filtering.
#' If set to NULL (the default) columns are automatically selected by
#' identifying columns that have NAs in the same rows as clone_col.
#' @param sep Separator used for storing per cell V(D)J data
#' @return List with two vectors, one containing columns with V(D)J data and
#' the other containing columns where separator has been detected.
#' @noRd
.get_vdj_cols <- function(df_in, clone_col, cols_in, sep) {

  # If clone_col and cols_in are both NULL, use all columns
  if (is.null(clone_col) && is.null(cols_in)) {
    cols_in <- colnames(df_in)
    cols_in <- cols_in[cols_in != ".cell_id"]
  }

  # If cols_in is NULL, identify columns with VDJ data based on NAs in
  # clone_col
  if (is.null(cols_in)) {
    cols_in <- dplyr::mutate(
      df_in,
      across(dplyr::everything(), is.na)
    )

    cols_in <- purrr::keep(
      cols_in,
      ~ identical(.x, pull(cols_in, clone_col))
    )

    cols_in <- colnames(cols_in)
  }

  # Identify columns to split based on sep
  sep_cols <- NULL

  if (!is.null(sep)) {
    sep_cols <- dplyr::select(df_in, all_of(cols_in))

    sep_cols <- purrr::keep(
      sep_cols,
      ~ any(purrr::map_lgl(na.omit(.x), grepl, pattern = sep))
    )

    sep_cols <- colnames(sep_cols)
  }

  # Return list of vectors
  res <- list(
    "vdj" = cols_in,
    "sep" = sep_cols
  )

  res
}


#' Attempt to coerce a character vector to a given class
#'
#' @param x Character vector to coerce
#' @param Class Name of the class to which x should be coerced
#' @return Coerced vector
#' @importFrom methods as
#' @noRd
.convert_char <- function(x, Class) {

  if (!is.character(x)) {
    return(x)
  }

  fn <- function() as(x, Class)

  suppressWarnings(ifelse(!is.na(fn()) | is.na(x), fn(), x))
}

