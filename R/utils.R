#' tibble imports
#'
#' @importFrom tibble tibble as_tibble column_to_rownames rownames_to_column
#' @noRd
NULL

#' purrr imports
#'
#' @importFrom purrr map imap map_dfr imap_dfr map_lgl map_chr map2_int map_dbl map_int iwalk pwalk
#' @importFrom purrr reduce keep is_empty is_function is_formula as_mapper
#' @noRd
NULL

#' stats imports
#'
#' @importFrom stats median complete.cases as.formula as.dist hclust cutree na.omit
#' @noRd
NULL


#' Global variables
#'
#' - CELL_COL is the column name to use for storing cell barcodes
#'
#' @noRd
CELL_COL <- ".cell_id"


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
    arg_lst <- tibble::as_tibble(arg_lst)

    return(arg_lst)
  }

  n <- 1

  purrr::pwalk(arg_lst, ~ {
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


#' Base R version of stringr::str_extract_all
#'
#' This is a little slower than stingr::str_extract_all
#'
#' @param string Input vector
#' @param pattern Regular expression
#' @return A character vector
#' @noRd
.str_extract_all <- function(string, pattern) {

  match_pos <- gregexpr(pattern, string, perl = TRUE)

  regmatches(string, match_pos)
}


#' Filter V(D)J data based on chain
#'
#' @param df_in data.frame
#' @param vdj_cols meta.data column(s) containing V(D)J data to filter based on
#' chain
#' @param chain Chain to use for filtering V(D)J data
#' @param chain_col meta.data column(s) containing chains for each cell
#' @param col_names A glue specification that describes how to name the output
#' columns, this can use {.col} to stand for the selected column name
#' @param empty_val Value to use when no chains match chain
#' @param allow_dups Allow cells to have multiple copies of the same chain. If
#' FALSE, cells must have exactly one copy of each chain specified. If a cell
#' has multiple copies of the same chain, V(D)J data for the cell is removed.
#' This filtering is only performed for list-cols.
#' @return filtered data.frame
#' @noRd
.filter_chains <- function(df_in, vdj_cols, chain, chain_col = "chains",
                           col_names = "{.col}", allow_dups = TRUE,
                           empty_val = NA) {

  if (is.null(chain)) return(df_in)

  chain <- unique(chain)

  # If all vdj_cols are not list-cols, filter data.frame normally
  is_lst <- purrr::map_lgl(df_in[, c(chain_col, vdj_cols)], is.list)

  if (all(!is_lst)) {
    res <- dplyr::filter(
      df_in,
      !!sym(chain_col) %in% chain
    )

    return(res)

  } else if (!all(is_lst)) {
    stop(
      "chain_col and vdj_cols cannot be a mix ",
      "of normal columns and list-cols."
    )
  }

  # Function to check/filter chains
  .map_fn <- function(x, chns) {

    if (length(x) != length(chns)) {
      stop(
        "Values in ", chain_col,  " are not the same length as vdj_cols, ",
        "are you using the correct chain_col?"
      )
    }

    r <- x[chns %in% chain]

    dup_chns <- chns[chns %in% chain]
    dup_chns <- any(duplicated(dup_chns))

    if (purrr::is_empty(r) || (!allow_dups && dup_chns)) r <- empty_val

    list(r)
  }

  # Filter chains
  res <- dplyr::rowwise(df_in)

  res <- dplyr::mutate(
    res,
    across(
      all_of(vdj_cols),
      ~ .map_fn(.x, !!sym(chain_col)),
      .names = col_names
    )
  )

  res <- dplyr::ungroup(res)

  if (!allow_dups) res <- tidyr::unnest(res, all_of(vdj_cols))

  res
}


#' Add meta.data to single cell object
#'
#' @param input Object containing single cell data
#' @param meta meta.data to add to object
#' @param row_col Column containing rownames to use for meta.data
#' @name .add_meta
#' @noRd
NULL

#' @rdname .add_meta
#' @return Object with added meta.data
#' @importFrom S4Vectors DataFrame
#' @noRd
.add_meta <- function(input, meta, row_col) {

  UseMethod(".add_meta", input)
}

.add_meta.default <- function(input, meta, row_col = CELL_COL) {
  if (!is.data.frame(meta)) {
    stop("meta.data must be a data.frame.")
  }

  tibble::column_to_rownames(meta, row_col)
}

.add_meta.Seurat <- function(input, meta, row_col = CELL_COL) {
  meta <- .prepare_meta(input, meta, row_col)

  input@meta.data <- meta

  input
}

.add_meta.SingleCellExperiment <- function(input, meta, row_col = CELL_COL) {
  meta <- .prepare_meta(input, meta, row_col)

  input@colData <- S4Vectors::DataFrame(meta)

  input
}

#' Prepare meta.data to add to object
#'
#' Use this when adding meta.data to a single-cell object. Do not need to
#' perform all of these checks when input object is a data.frame or tibble.
#'
#' @rdname .add_meta
#' @noRd
.prepare_meta <- function(input, meta, row_col = CELL_COL) {
  if (!is.data.frame(meta)) {
    stop("meta.data must be a data.frame.")
  }

  is_lst <- purrr::map_lgl(meta, is.list)

  if (any(is_lst)) {
    stop("meta.data cannot include list-cols.")
  }

  if (!is.null(input) && !identical(colnames(input), meta[[row_col]])) {
    stop(
      "To add meta.data to an object, the meta.data must contain the same ",
      "cells as the target object."
    )
  }

  meta <- tibble::column_to_rownames(meta, row_col)

  meta
}


#' Pull meta.data from single cell object
#'
#' @param input Object containing single cell data
#' @param row_col New column to store meta.data rownames
#' @return tibble containing meta.data pulled from object
#' @importFrom SingleCellExperiment colData
#' @noRd
.get_meta <- function(input, row_col) {

  UseMethod(".get_meta", input)
}

.get_meta.default <- function(input, row_col = CELL_COL) {

  .to_tibble(input, row_col)
}

.get_meta.Seurat <- function(input, row_col = CELL_COL) {

  .get_meta(input@meta.data, row_col)
}

.get_meta.SingleCellExperiment <- function(input, row_col = CELL_COL) {
  dat <- SingleCellExperiment::colData(input)
  dat <- as.data.frame(dat)

  .get_meta(dat, row_col)
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
.merge_meta <- function(input, meta, by = CELL_COL) {

  if (is.null(input)) {
    return(meta)
  }

  # Join meta.data
  # remove columns already present in input object to prevent duplicates
  # this mimics behavior of Seurat::AddMetaData
  meta <- .get_meta(meta, row_col = by)

  rm_cols <- colnames(meta)
  rm_cols <- rm_cols[!rm_cols %in% by]

  obj_meta <- .get_meta(input)
  obj_meta <- dplyr::select(obj_meta, !any_of(rm_cols))

  meta <- dplyr::left_join(obj_meta, meta, by = by)
  res  <- .add_meta(input, meta = meta)

  res
}


#' Collapse meta.data list-cols into strings
#'
#' @param df_in data.frame to modify
#' @param sep_cols list-cols to collapse into strings, if NULL all list-cols
#' in the data.frame will be collapsed
#' @param sep Separator to use for collapsing list-cols
#' @return data.frame
#' @noRd
.nest_vdj <- function(df_in, sep_cols = NULL, sep = ";") {

  if (is.null(sep_cols)) {
    sep_cols <- colnames(df_in)
  }

  nest_cols <- purrr::map_lgl(df_in[, sep_cols], is.list)
  nest_cols <- sep_cols[nest_cols]

  if (purrr::is_empty(nest_cols)) {
    return(df_in)
  }

  # paste0 will convert NA to "NA", to avoid this first check for NA before
  # collapsing
  res <- dplyr::rowwise(df_in)

  res <- dplyr::mutate(
    res,
    across(
      all_of(nest_cols),
      ~ ifelse(all(is.na(.x)), NA, paste0(.x, collapse = sep))
    )
  )

  res <- dplyr::ungroup(res)

  res
}


#' Split meta.data columns into list-cols
#'
#' @param df_in data.frame to modify
#' @param sep_cols Columns to split based on sep
#' @param unnest Should columns be unnested after splitting into vectors
#' @param sep Separator used for storing per cell V(D)J data
#' @return data.frame with V(D)J data
#' @importFrom readr guess_parser
#' @noRd
.unnest_vdj <- function(df_in, sep_cols, unnest = TRUE, sep = ";") {

  df_in <- as_tibble(df_in)

  # Get types to use for coercing columns
  # use first 100 rows containing V(D)J data
  typs <- dplyr::select(df_in, all_of(sep_cols))
  typs <- dplyr::filter(typs, if_all(all_of(sep_cols), ~ !is.na(.x)))
  typs <- utils::head(typs, 100)

  typs <- tidyr::separate_rows(typs, all_of(sep_cols), sep = sep)

  typs <- purrr::map(typs, readr::guess_parser)
  typs <- purrr::map_chr(typs, ~ paste0("as.", .x))

  # Split sep_cols and convert types
  # tidyr::separate_rows with convert = TRUE is slower than determining types
  # and then converting
  # strsplit is faster than str_split_n used by separate_rows
  .str_split_convert <- function(x, pattern, fn) {
    res <- map(x, ~ {
      .x <- strsplit(.x, pattern, fixed = TRUE)
      .x <- unlist(.x)

      do.call(fn, list(x = .x))
    })

    res
  }

  if (!unnest) {
    res <- dplyr::mutate(
      df_in,
      across(
        all_of(sep_cols),
        ~ .str_split_convert(.x, sep, typs[[cur_column()]])
      )
    )

    return(res)
  }

  res <- purrr::modify_at(df_in, sep_cols, strsplit, split = sep, fixed = TRUE)
  res <- tidyr::unchop(res, all_of(sep_cols))

  res <- mutate(
    res,
    across(
      all_of(sep_cols),
      ~ do.call(typs[[cur_column()]], list(x = .x))
    )
  )

  res
}

#' Identify columns with V(D)J data
#'
#' @param df_in data.frame
#' @param clone_col Column containing clonotype IDs to use for identifying
#' columns with V(D)J data. If both clone_col and cols_in are NULL, all columns
#' are checked. This hurts performance when df_in has a lot of columns.
#' @param cols_in Columns containing V(D)J data. If NULL data are selected by
#' identifying columns that have NAs in the same rows as clone_col.
#' @param sep Separator used for storing per cell V(D)J data. This is used to
#' identify columns containing per-chain data that can be unnested.
#' @param cell_col Column containing cell IDs
#' @return List with two vectors, one containing columns with V(D)J data and
#' the other containing columns where separator has been detected.
#' @noRd
.get_vdj_cols <- function(df_in, clone_col, cols_in, sep,
                          cell_col = CELL_COL) {

  # Check clone_col
  if (length(clone_col) > 1) {
    stop("Provide a single value for clonotype column.")
  }

  if (is.null(cols_in) && !is.null(clone_col) && !clone_col %in% colnames(df_in)) {
    clone_col <- NULL

    warning(
      "The provided clonotype column is not present in input data, provide a ",
      "column containing clonotype IDs to increase performance."
    )
  }

  # If clone_col and cols_in are both NULL, use all columns
  if (is.null(clone_col) && is.null(cols_in)) {
    cols_in <- colnames(df_in)
    cols_in <- cols_in[cols_in != cell_col]
  }

  # If cols_in is NULL, identify columns with V(D)J data based on NAs in
  # clone_col
  if (is.null(cols_in)) {
    cols_in <- dplyr::mutate(
      df_in,
      across(dplyr::everything(), is.na)
    )

    cols_in <- purrr::keep(
      cols_in,
      ~ identical(.x, dplyr::pull(cols_in, clone_col))
    )

    cols_in <- colnames(cols_in)
  }

  # Identify columns to unnest based on sep
  # check first 1000 non-NA rows of every column
  # detect is faster for columns containing sep but slower when column does not
  # contain sep
  # ~ !is.null(detect(x, ~ grepl(sep, .x)))
  sep_cols <- NULL

  if (!is.null(sep)) {
    sep_cols <- dplyr::select(df_in, all_of(cols_in))

    sep_cols <- purrr::keep(sep_cols, ~ {
      x <- head(stats::na.omit(.x), 1000)

      any(purrr::map_lgl(x, grepl, pattern = sep, fixed = TRUE))
    })

    sep_cols <- colnames(sep_cols)

    if (purrr::is_empty(sep_cols)) {
      sep_cols <- NULL
    }
  }

  # Return list of vectors
  res <- list(
    "vdj" = cols_in,
    "sep" = sep_cols
  )

  res
}

