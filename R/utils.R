#' tibble imports
#'
#' @importFrom tibble tibble as_tibble column_to_rownames rownames_to_column
#' @noRd
NULL

#' purrr imports
#'
#' @importFrom purrr map imap map_dfr imap_dfr map_lgl map_chr map2_int map_dbl
#' @importFrom purrr map_int iwalk pwalk reduce keep is_empty is_function
#' @importFrom purrr is_formula as_mapper
#' @noRd
NULL

#' stats imports
#'
#' @importFrom stats median complete.cases as.formula as.dist hclust cutree
#' @importFrom stats na.omit sd
#' @noRd
NULL

#' cli imports
#'
#' @importFrom cli cli_abort cli_warn cli_alert
#' @noRd
NULL


#' Helper to test all combinations of provided arguments
#'
#' @param arg_lst Named list of arguments to test
#' @param .fn Function to test
#' @param desc Description to pass to test_that
#' @param chk Function or expression to use for testing. If an expression is
#' passed, results from .fn can be referred to with .res.
#' @param dryrun Do not run tests, just return table of arguments that will be
#' tested
#' @return Output from test_that
#' @noRd
test_all_args <- function(arg_lst, .fn, desc, chk, dryrun = FALSE) {
  arg_lst    <- expand.grid(arg_lst, stringsAsFactors = FALSE)
  arg_lst$.n <- seq_len(nrow(arg_lst))

  if (dryrun) {
    arg_lst <- tibble::as_tibble(arg_lst)

    return(arg_lst)
  }

  purrr::pwalk(arg_lst, ~ {
    test_args    <- list(...)
    n            <- test_args$.n
    test_args$.n <- NULL

    test_that(paste(desc, n), {
      if (is.call(chk)) {
        .res <- lift(.fn)(test_args)

        return(eval(chk))
      }

      chk(lift(.fn)(test_args))
    })
  })
}


#' Lift the domain of a function
#'
#' To replace `purrr::lift_dl()`, which is deprecated in purrr v1.0.0
#'
#' @param ..f A function to lift.
#' @param ... Default arguments for `..f`. These will be
#' evaluated only once, when the lifting factory is called.
#' @param .unnamed This prevents matching of arguments by name, arguments will
#' be be matched by position instead.
#' @return A function.
#' @noRd
lift <- function(..f, ..., .unnamed = FALSE) {
  force(..f)

  defaults <- list(...)

  function(.x = list(), ...) {
    if (.unnamed) .x <- unname(.x)

    do.call("..f", c(.x, defaults, list(...)))
  }
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
#' This will remove V(D)J data for chains that do not match the provided chain.
#' NAs will be added for cells that do not have any matching chains. Filtered
#' data_cols will be returned as list-cols. The chain_col will not be filtered,
#' but will be included in the results as a list-col.
#'
#' @param df_in data.frame
#' @param data_cols meta.data column(s) containing V(D)J data to filter based on
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
.filter_chains <- function(df_in, data_cols, chain,
                           chain_col = global$chain_col, col_names = "{.col}",
                           allow_dups = TRUE, empty_val = NA) {

  if (is.null(chain)) return(df_in)

  chain <- unique(chain)

  # If all data_cols are not list-cols, filter data.frame normally
  is_lst <- purrr::map_lgl(df_in[, c(chain_col, data_cols)], is.list)

  if (all(!is_lst)) {
    res <- dplyr::filter(df_in, !!sym(chain_col) %in% chain)

    return(res)

  } else if (!all(is_lst)) {
    bad_cols <- names(is_lst[!is_lst])

    cli::cli_abort(
      "Some columns do not contain per-chain V(D)J data, can only filter
       based on `chain` if all `data_cols` contain per-chain data: {bad_cols}"
    )
  }

  # Function to check/filter chains
  .filt_fn <- function(x, chns) {
    if (length(x) != length(chns)) {
      cli::cli_abort(
        "The number of chains per cell present in the {chain_col} column differs
         from the number of per-chain values present in `data_cols`, are you
         using the correct `chain_col`?"
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
      all_of(data_cols),
      ~ .filt_fn(.x, !!sym(chain_col)),
      .names = col_names
    )
  )

  res <- dplyr::ungroup(res)

  if (!allow_dups) res <- tidyr::unnest(res, all_of(data_cols))

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
#' @importFrom methods slot
#' @noRd
.add_meta <- function(input, meta, row_col) {

  UseMethod(".add_meta", input)
}

.add_meta.default <- function(input, meta, row_col = global$cell_col) {
  if (!is.data.frame(meta)) {
    cli::cli_abort("meta.data must be a data.frame")
  }

  tibble::column_to_rownames(meta, row_col)
}

.add_meta.Seurat <- function(input, meta, row_col = global$cell_col) {
  meta <- .prepare_meta(input, meta, row_col)

  methods::slot(input, "meta.data") <- meta

  input
}

.add_meta.SingleCellExperiment <- function(input, meta,
                                           row_col = global$cell_col) {
  meta <- .prepare_meta(input, meta, row_col)

  methods::slot(input, "colData") <- S4Vectors::DataFrame(meta)

  input
}

#' Prepare meta.data to add to object
#'
#' Use this when adding meta.data to a single-cell object. Do not need to
#' perform all of these checks when input object is a data.frame or tibble.
#'
#' @rdname .add_meta
#' @noRd
.prepare_meta <- function(input, meta, row_col = global$cell_col) {

  if (!is.data.frame(meta)) {
    cli::cli_abort("meta.data must be a data.frame")
  }

  is_lst <- purrr::map_lgl(meta, is.list)

  if (any(is_lst)) {
    cli::cli_abort("meta.data cannot include list-cols")
  }

  if (!is.null(input) && !identical(colnames(input), meta[[row_col]])) {
    cli::cli_abort(
      "meta.data does not contain the same cells as the object, check your cell
       barcodes"
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

.get_meta.default <- function(input, row_col = global$cell_col) {

  .to_tibble(input, row_col)
}

.get_meta.Seurat <- function(input, row_col = global$cell_col) {

  .get_meta(input@meta.data, row_col)
}

.get_meta.SingleCellExperiment <- function(input,
                                           row_col = global$cell_col) {
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
.merge_meta <- function(input, meta, by = global$cell_col) {

  if (is.null(input)) return(meta)

  # Join meta.data
  # remove columns already present in input object to prevent duplicates
  # this mimics behavior of Seurat::AddMetaData()
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
.nest_vdj <- function(df_in, sep_cols = NULL, sep = global$sep) {

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
#' @importFrom utils head
#' @noRd
.unnest_vdj <- function(df_in, sep_cols, unnest = TRUE, sep = global$sep) {

  df_in <- tibble::as_tibble(df_in)

  # Check first 100 rows containing V(D)J data
  typs <- dplyr::select(df_in, all_of(sep_cols))
  typs <- dplyr::filter(typs, if_all(all_of(sep_cols), ~ !is.na(.x)))
  typs <- utils::head(typs, 100)

  # Check number of chains
  if (length(sep_cols) > 1) {
    chain_count <- purrr::map(typs, stringr::str_count, sep)

    i1 <- chain_count[[1]]
    rest <- chain_count[-1]

    for (x in rest) {
      if (!identical(i1, x)) {
        cli::cli_abort(c(
          "!" = "The provided columns ({sep_cols}) do not contain matching
                 per-chain data, check that all columns contain data for the
                 same set of chains. You may receive this error for several
                 reasons:",
          "*" = "The wrong `chain_col` was provided",
          "*" = "The selected columns correspond to different V(D)J data types
                 (i.e. BCR and TCR)",
          "*" = "The per-chain data for one of the provided columns has been
                 previously filtered"
        ))
      }
    }
  }

  # Get types to use for coercing columns
  typs <- tidyr::separate_rows(typs, all_of(sep_cols), sep = sep)
  typs <- purrr::map(typs, readr::guess_parser)
  typs <- typs[typs != "character"]
  typs <- purrr::map_chr(typs, ~ paste0("as.", .x))

  typ_cols <- names(typs)

  # Split sep_cols and convert types
  # tidyr::separate_rows with convert = TRUE is slower than determining types
  # and then converting
  # strsplit is faster than str_split_n used by separate_rows
  .str_convert      <- function(x, fn) do.call(fn, list(x = x))
  .str_convert_list <- function(l, fn) map(l, .str_convert, fn = fn)

  res <- purrr::modify_at(
    df_in, sep_cols,
    ~ strsplit(as.character(.x), split = sep, fixed = TRUE)
  )

  if (!unnest) {
    res <- dplyr::mutate(
      res,
      across(all_of(typ_cols), ~ .str_convert_list(.x, typs[[cur_column()]]))
    )

    return(res)
  }

  res <- tidyr::unchop(res, all_of(sep_cols))

  res <- mutate(
    res, across(all_of(typ_cols), ~ .str_convert(.x, typs[[cur_column()]]))
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
                          cell_col = global$cell_col) {

  # Check clone_col
  no_clone_col <- is.null(cols_in) && !is.null(clone_col) &&
    !clone_col %in% colnames(df_in)

  if (no_clone_col) {
    clone_col <- NULL

    cli::cli_warn(
      "The {clone_col} column is not present in the input data, provide a
       column containing clonotype IDs to increase performance"
    )
  }

  # If clone_col and cols_in are both NULL, use all columns
  # columns should not include cell_col
  if (is.null(clone_col) && is.null(cols_in)) cols_in <- colnames(df_in)

  # If cols_in is NULL, identify columns with V(D)J data based on NAs in
  # clone_col
  if (is.null(cols_in)) {
    cols_in <- dplyr::mutate(df_in, across(dplyr::everything(), is.na))

    cols_in <- purrr::keep(
      cols_in, ~ identical(.x, dplyr::pull(cols_in, clone_col))
    )

    cols_in <- colnames(cols_in)
  }

  cols_in <- cols_in[cols_in != cell_col]

  # Identify columns to unnest based on sep
  # check first 1000 non-NA rows of every column
  # detect is faster for columns containing sep but slower when column does not
  # contain sep
  # ~ !is.null(detect(x, ~ grepl(sep, .x)))
  sep_cols <- NULL

  if (!is.null(sep)) sep_cols <- .detect_sep(df_in, cols_in, sep, n_rows = 1000)

  # Return list of vectors
  res <- list(
    "vdj" = cols_in,
    "sep" = sep_cols
  )

  res
}


#' Check that columns are present in object
#'
#' This also checks that the provided chain is found in chain_col
#'
#' @param input input object
#' @param ... Columns to check for in input
#' @param chain Chain
#' @param chain_col Column containing chains
#' @param list_avail List available columns in error message
#' @return Error when column(s) are missing
#' @noRd
.check_obj_cols <- function(input, ..., chain = NULL,
                            chain_col = global$chain_col, list_avail = FALSE) {

  dat      <- .get_meta(input, row_col = global$cell_col)
  dat_cols <- colnames(dat)

  cols <- c(...)

  if (!is.null(chain)) {
    if (is.null(chain_col)) {
      cli::cli_abort("`chain_col` must be provided when `chain` is specified")
    }

    if (chain_col %in% dat_cols) {
      chns <- dat[[chain_col]]
      chns <- chns[!is.na(chns)]
      chns <- any(grepl(paste0(chain, collapse = "|"), chns))

      if (!chns) {
        cli::cli_abort(
          "The specified chain ({chain}) is not present in {chain_col}"
        )
      }
    }

    cols <- c(cols, chain_col)
  }

  cols <- unique(cols)
  chk  <- cols %in% dat_cols

  if (!all(chk)) {
    bad <- cols[!chk]

    msg <- "Some columns are not present in input data: {bad}"

    if (list_avail) {
      dat_cols <- paste0(dat_cols, collapse = ", ")

      msg <- c(
        "!" = msg,
        "i" = "Available columns include: {dat_cols}"
      )
    }

    cli::cli_abort(msg)
  }
}


#' Check type and length of argument
#'
#' @param arg Name of argument to check
#' @param Class Must be one of the specified classes, if arg is a list, each
#' item will get checked for the provided classes. To only check that arg is a
#' list, just set Class to 'list'.
#' @param len_one Should argument be length 1
#' @param allow_null Can argument be NULL
#' @param envir Environment to use
#' @return Error when argument is not expected class and length
#' @noRd
.check_arg <- function(arg, Class = "character", len_one = TRUE,
                       allow_null = FALSE, envir) {

  val <- eval(sym(arg), envir = envir)
  len <- length(val)

  if (!allow_null && is.null(val)) cli::cli_abort("`{arg}` cannot be `NULL`")

  if (is.null(val)) return(invisible())

  if (len_one && len != 1) cli::cli_abort("`{arg}` must be length 1")

  if (is.list(val) && identical(Class, "list")) return(invisible())

  # Check each item of list
  if (!is.list(val)) val <- list(val)

  purrr::walk(val, ~ {
    val <- .x

    # val must be one of the classes specified to Class
    chk <- purrr::map_lgl(Class, ~ methods::is(val, .x))

    if (!any(chk)) cli::cli_abort("`{arg}` must be class {.or {Class}}")
  })
}

.check_args <- function(..., arg_classes = global$arg_classes) {
  envir   <- parent.frame()
  new_cls <- list(...)

  arg_classes[names(new_cls)] <- new_cls

  arg_classes <- purrr::imap(arg_classes, ~ {
    .x$arg <- .y
    .x
  })

  arg_classes <- arg_classes[names(arg_classes) %in% ls(envir = envir)]

  purrr::walk(arg_classes, ~ purrr::pmap(.x, .check_arg, envir = envir))
}


#' Check possible values for argument
#'
#' @param ... Name value pairs providing argument name and possible values
#' @param .internal Internal error
#' @return Error if argument does not match one of the provided values
#' @noRd
.check_possible_values <- function(..., .internal = FALSE) {
  envir <- parent.frame()
  vals  <- list(...)

  purrr::iwalk(vals, ~ {
    val <- eval(sym(.y), envir = envir)

    if (!all(val %in% .x)) {
      cli::cli_abort("`{.y}` must be {.or {.x}}", .internal = .internal)
    }
  })
}


#' Check that package is installed
#'
#' @importFrom rlang is_installed
#' @noRd
.check_packages <- function(pkgs, db = "CRAN") {
  chks <- purrr::map_lgl(pkgs, rlang::is_installed)
  pkgs <- paste0("\'", pkgs, "\'")

  missing <- pkgs[!chks]

  if (any(!chks)) {
    cli::cli_abort(paste0(
      "Package{?s} {cli::qty(missing)} must be installed to use this function.
       Th{?is/ese} package{?s} {?is/are} available on ", db, "."
    ))
  }
}
