#' Add meta.data to single cell object
#'
#' @export
.add_meta <- function(input, ...) {
  UseMethod(".add_meta", input)
}

#' @rdname dot-add_meta
#' @param input Object containing single cell data
#' @param meta meta.data to add to object
#' @param row_col Column containing meta.data rownames
#' @param ... Arguments passed to other methods
#' @return Object with added meta.data
#' @export
.add_meta.default <- function(meta, ..., row_col = ".cell_id") {

  res <- tibble::column_to_rownames(meta, row_col)

  res
}

#' @rdname dot-add_meta
#' @export
.add_meta.Seurat <- function(input, meta, row_col = ".cell_id") {

  meta <- tibble::column_to_rownames(meta, row_col)

  input@meta.data <- meta

  input
}

#' @rdname dot-add_meta
#' @importFrom S4Vectors DataFrame
#' @export
.add_meta.SingleCellExperiment <- function(input, meta, row_col = ".cell_id") {

  meta <- tibble::column_to_rownames(meta, row_col)

  input@colData <- S4Vectors::DataFrame(meta)

  input
}


#' Add meta.data to single cell object
#'
#'@export
.get_meta <- function(input, ...) {
  UseMethod(".get_meta", input)
}

#' @rdname dot-get_meta
#' @param input Object containing single cell data
#' @param row_col New column to store meta.data rownames
#' @param ... Arguments passed to other methods
#' @return tibble containing meta.data pulled from object
#' @export
.get_meta.default <- function(input, row_col = ".cell_id") {

  meta <- tibble::as_tibble(input, rownames = row_col)

  meta
}

#' @rdname dot-get_meta
#' @export
.get_meta.Seurat <- function(input, row_col = ".cell_id") {

  meta <- tibble::as_tibble(input@meta.data, rownames = row_col)

  meta
}

#' @rdname dot-get_meta
#' @export
.get_meta.SingleCellExperiment <- function(input, row_col = ".cell_id") {

  meta <- tibble::as_tibble(input@colData, rownames = row_col)

  meta
}


#' Summarize values for chains
#'
#' Summarize values present for each column provided to the data_cols argument.
#' For each cell, the function(s) provided will be applied to each unique label
#' in chain_col.
#'
#' @param input Object containing single cell data
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
  res <- dplyr::mutate(meta, across(
    all_of(c(data_cols, chain_col)),
    ~ strsplit(as.character(.x), sep)
  ))

  res <- tidyr::unnest(res, cols = all_of(c(data_cols, chain_col)))

  # Summarize data_cols for each chain present for the cell
  res <- dplyr::mutate(res, across(
    all_of(data_cols),
    ~ .convert_char(.x, as.numeric)
  ))

  grp_cols <- c(".cell_id", chain_col, include_cols)
  res      <- dplyr::group_by(res, !!!syms(grp_cols))

  res <- dplyr::summarize(
    res,
    across(all_of(data_cols), fn),
    .groups = "drop"
  )

  res
}


#' Split V(D)J meta.data columns into vectors
#'
#' @param df_in data.frame to modify
#' @param sep Separator used for storing per cell V(D)J data
#' @param sep_cols Columns to split based on sep
#' @param num_cols Columns to convert to numeric
#' @param lgl_cols Columns to convert to logical
#' @return data.frame
.split_vdj <- function(df_in, sep = ";", sep_cols, num_cols = c("umis", "reads"),
                       lgl_cols = c("productive", "full_length")) {

  # Split columns into vectors
  res <- dplyr::mutate(df_in, !!!syms(sep_cols))

  res <- dplyr::mutate(res, across(
    all_of(unname(sep_cols)),
    ~ strsplit(as.character(.x), sep)
  ))

  # Convert to numeric or logical
  res <- dplyr::mutate(
    res,
    across(
      all_of(num_cols),
      map, ~ .convert_char(.x, as.numeric)
    ),
    across(
      all_of(lgl_cols),
      map, ~ .convert_char(.x, as.logical)
    )
  )

  res
}


#' Identify columns with V(D)J data
#'
#' @param df_in data.frame
#' @param clone_col Column containing clonotype IDs to use for identifying
#' columns with V(D)J data. If NULL all columns are used.
#' @param cols_in meta.data columns containing VDJ data to use for filtering.
#' If set to NULL (recommended) columns are automatically selected.
#' @param sep Separator used for storing per cell V(D)J data
#' @return list of vectors containing columns with V(D)J data and sep
.get_vdj_cols <- function(df_in, clone_col, cols_in, sep) {

  # Identify columns with VDJ data based on NAs in clonotype_col
  # If no clonotype_col is given use all columns
  if (is.null(clone_col)) {
    cols_in <- colnames(df_in)
    cols_in <- cols_in[cols_in != ".cell_id"]
  }

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


#' Attempt to convert character vector using provided function
#'
#' @param x Character vector to convert
#' @param fn Function to try
#' @return Value converted using fn
.convert_char <- function(x, fn) {
  if (!is.character(x)) {
    return(x)
  }

  suppressWarnings(ifelse(!is.na(fn(x)) | is.na(x), fn(x), x))
}


#' Helper to test all combinations of provided arguments
#'
#' @param arg_lst Named list of arguments to test
#' @param .fn Function to test
#' @param desc Description to pass to test_that
#' @param chk Function or expression to using for testing
#' @param dryrun Do not run tests, just return table of arguments that will be
#' tested
#' @return Output from test_that
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


