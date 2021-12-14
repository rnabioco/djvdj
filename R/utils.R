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

  res <- .unnest_vdj(
    meta,
    sep         = sep,
    sep_cols    = c(data_cols, chain_col),
    coerce_cols = coerce_cols,
    expand      = TRUE
  )

  # Filter for chain
  if (!is.null(chain)) {
    plt_dat <- dplyr::filter(plt_dat, !!sym(chain_col) %in% chain)
  }

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


#' Fetch V(D)J data
#'
#' Fetch per-chain V(D)J data from object. Within the object meta.data, each
#' row represents a single cell and can include information for multiple
#' chains. This function returns an unnested data.frame where each row
#' represents a single chain. This is useful for plotting per-chain metrics
#' such as CDR3 length or the number of insertions/deletions.
#'
#' @param input Single cell object or data.frame containing V(D)J data. If a
#' data.frame is provided, the cell barcodes should be stored as row names.
#' @param clonotype_col meta.data column containing clonotype IDs. This column
#' is used to determine which cells have V(D)J data.
#' @param vdj_cols meta.data columns containing V(D)J data to fetch. If set to
#' NULL columns are automatically selected by identifying columns that have NAs
#' in the same rows as clonotype_col.
#' @param sep Separator used for storing per cell V(D)J data
#' @param unnest If FALSE, a nested data.frame is returned where each row
#' represents a cell and V(D)J data is stored as list-cols
#' @return data.frame containing per-chain V(D)J data
#' @export
fetch_vdj <- function(input, clonotype_col = "clonotype_id", vdj_cols = NULL, sep = ";",
                      unnest = TRUE) {

  # Format input data
  meta <- .get_meta(input)

  # Identify columns with V(D)J data
  col_list <- .get_vdj_cols(
    df_in     = meta,
    clone_col = clonotype_col,
    cols_in   = vdj_cols,
    sep       = sep
  )

  vdj_cols <- col_list$vdj
  sep_cols <- col_list$sep

  if (purrr::is_empty(sep_cols)) {
    warning(
      "The separator '", sep, "' was not identified in any columns specified ",
      "by vdj_cols, the unmodified meta.data will be returned"
    )

    return(meta)
  }

  # Unnest V(D)J data
  res <- .unnest_vdj(
    meta,
    sep      = sep,
    sep_cols = sep_cols,
    unnest   = unnest
  )

  res
}


#' Summarize numeric V(D)J data
#'
#' @param input Single cell object or data.frame containing V(D)J data. If a
#' data.frame is provided, the cell barcodes should be stored as row names.
#' @param vdj_cols meta.data columns containing V(D)J data to summarize for
#' each cell
#' @param fn Function or list of funcitons to use for summarizing vdj_cols
#' @param chain Chain to use for summarizing V(D)J data
#' @param chain_col meta.data column(s) containing chains for each cell
#' @param sep Separator used for storing per cell V(D)J data
#' @return data.frame containing V(D)J data summarized for each cell
#' @export
summarize_vdj <- function(input, vdj_cols, fn = mean, chain = NULL, chain_col = "chains",
                          sep = ";") {

  # Fetch V(D)J data
  fetch_cols <- vdj_cols

  if (!is.null(chain)) {
    fetch_cols <- c(vdj_cols, chain_col)
  }

  res <- fetch_vdj(
    input,
    clonotype_col = NULL,
    vdj_cols      = fetch_cols,
    sep           = sep,
    unnest        = FALSE
  )

  # Column names
  fn_nm <- as.character(substitute(fn))
  nms   <- "{fn_nm}_{.col}"

  # Filter chains
  if (!is.null(chain)) {
    .map_fn <- function(x, chns) {
      purrr::map2(x, chns, ~ {
        r <- .x[.y %in% chain]

        ifelse(is_empty(r), as.numeric(NA), r)
      })
    }

    res <- dplyr::mutate(
      res,
      across(
        all_of(vdj_cols),
        ~ .map_fn(.x, !!sym(chain_col)),
        .names = nms
      )
    )

    vdj_cols <- paste0(fn_nm, "_", vdj_cols)
    nms      <- NULL
  }

  # Summarize columns
  .map_fn <- function(x) {
    purrr::map_dbl(x, fn)
  }

  res <- dplyr::mutate(
    res,
    across(
      all_of(vdj_cols),
      .map_fn,
      .names = nms
    )
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
.unnest_vdj <- function(df_in, sep = ";", sep_cols, unnest = FALSE,
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

      res <<- dplyr::mutate(
        res,
        !!clmn := map(!!clmn, .convert_char, Class)
      )
    })
  }

  # Unnest columns
  if (unnest) {
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

  # If cols_in is NULL, identify columns with V(D)J data based on NAs in
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

  # Identify columns to unnest based on sep
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
.convert_char <- function(x, Class = NULL) {

  if (!is.character(x)) {
    return(x)
  }

  fn <- function() as(x, Class)

  suppressWarnings(ifelse(!is.na(fn()) | is.na(x), fn(), x))
}

