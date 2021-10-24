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


#' Add meta.data to single cell object
#'
.add_meta <- function(input, meta, row_col) {

  UseMethod(".add_meta", input)

}

#' @rdname dot-add_meta
#' @param input Object containing single cell data
#' @param meta meta.data to add to object
#' @param row_col Column containing meta.data rownames
#' @return Object with added meta.data
.add_meta.default <- function(input, meta, row_col = ".cell_id") {

  tibble::column_to_rownames(meta, row_col)

}

#' @rdname dot-add_meta
.add_meta.Seurat <- function(input, meta, row_col = ".cell_id") {

  meta <- tibble::column_to_rownames(meta, row_col)

  input@meta.data <- meta

  input
}

#' @rdname dot-add_meta
#' @importFrom S4Vectors DataFrame
.add_meta.SingleCellExperiment <- function(input, meta, row_col = ".cell_id") {

  meta <- tibble::column_to_rownames(meta, row_col)

  input@colData <- S4Vectors::DataFrame(meta)

  input
}


#' Pull meta.data from single cell object
#'
.get_meta <- function(input, row_col) {

  UseMethod(".get_meta", input)

}

#' @rdname dot-get_meta
#' @param input Object containing single cell data
#' @param row_col New column to store meta.data rownames
#' @return tibble containing meta.data pulled from object
.get_meta.default <- function(input, row_col = ".cell_id") {

  .to_tibble(input, row_col)

}

#' @rdname dot-get_meta
.get_meta.Seurat <- function(input, row_col = ".cell_id") {

  .to_tibble(input@meta.data, row_col)

}

#' @rdname dot-get_meta
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

