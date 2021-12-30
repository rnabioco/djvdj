#' Modify V(D)J data in object
#'
#' @param input Single cell object or data.frame containing V(D)J data. If a
#' data.frame is provided, cell barcodes should be stored as row names.
#' @param ... Name-value pairs to use for creating or modifying meta.data
#' columns
#' @param clonotype_col meta.data column containing clonotype IDs. This is used
#' to identify columns containing V(D)J data.
#' @param vdj_cols meta.data columns containing V(D)J data to modify. If NULL,
#' data are automatically selected by identifying columns that have NAs in the
#' same rows as clonotype_col.
#' @param return_df Return results as a data.frame. If FALSE, results will be
#' added to the input object.
#' @param sep Separator used for storing per cell V(D)J data
#' @return Object with mutated meta.data
#' @export
mutate_vdj <- function(input, ..., clonotype_col = "clonotype_id", vdj_cols = NULL,
                       return_df = FALSE, sep = ";") {

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

  # Allow sep to be NULL so user can skip .unnest_vdj
  if (purrr::is_empty(sep_cols) && !is.null(sep)) {
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
  }

  # Mutate meta.data
  vdj <- dplyr::mutate(vdj, ...)
  vdj <- dplyr::ungroup(vdj)

  # Re-nest list-cols
  res <- .nest_vdj(vdj, sep_cols = sep_cols, sep = sep)

  if (return_df) {
    input <- res
  }

  res <- .add_meta(input, meta = res)

  res
}


#' Summarize V(D)J data for each cell
#'
#' Summarize values for each cell using a function or purrr-style lambda. This
#' is useful since some V(D)J metrics include data for each chain. Per-chain
#' data will be summarized, which is helpful for plotting these metrics per
#' cell.
#'
#' @param input Single cell object or data.frame containing V(D)J data. If a
#' data.frame is provided, the cell barcodes should be stored as row names.
#' @param vdj_cols meta.data column(s) containing V(D)J data to summarize for
#' each cell
#' @param fn Function to apply to each selected column, possible values can be
#' either a function, e.g. mean, or a purrr-style lambda, e.g. ~ mean(.x,
#' na.rm = TRUE). If NULL, the mean will be calculated for numeric values,
#' non-numeric columns will be combined into a single string.
#' @param ... Additional arguments to pass to fn
#' @param chain Chain to use for summarizing V(D)J data
#' @param chain_col meta.data column(s) containing chains for each cell
#' @param col_names A glue specification that describes how to name the output
#' columns, this can use {.col} to stand for the selected column name
#' @param return_df Return results as a data.frame. If FALSE, results will be
#' added to the input object.
#' @param sep Separator used for storing per cell V(D)J data
#' @return data.frame containing V(D)J data summarized for each cell
#' @importFrom glue glue
#' @export
summarize_vdj <- function(input, vdj_cols, fn = NULL, ..., chain = NULL, chain_col = "chains",
                          sep = ";", col_names = "{.col}", return_df = FALSE) {

  # Names of new columns
  new_cols <- gsub("\\{.col\\}", "{vdj_cols}", col_names)
  new_cols <- glue::glue(new_cols)

  # Fetch V(D)J data
  # * this is a key performance bottleneck as the length of vdj_cols increases
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

  # Set default fn
  # * default is mean if all vdj_col are numeric
  # * if !is.null(chain), collapse vdj_cols into string
  # * if is.null(chain), return input unchanged
  if (is.null(fn)) {
    is_num <- purrr::map_lgl(
      res[, vdj_cols],
      ~ all(purrr::map_lgl(.x, is.numeric))
    )

    is_num <- all(is_num)

    fn <- mean

    if (!is_num) {
      if (is.null(chain)) {
        return(input)
      }

      fn <- ~ paste0(.x, collapse = sep)
    }
  }

  # Filter chains
  if (!is.null(chain)) {
    prfx <- "FILT_"

    res <- .filter_chains(
      res,
      vdj_cols  = vdj_cols,
      chain     = chain,
      chain_col = chain_col,
      col_names = paste0(prfx, "{.col}"),
      empty_val = NA
    )

    # Add prefix to vdj_cols so temporary columns are used
    vdj_cols <- paste0(prfx, vdj_cols)

    # Set col_names so prefix is removed from columns
    col_names <- gsub(
      "\\{.col\\}",
      paste0('{sub(\\"^', prfx, '\\", "", .col)}'),
      col_names
    )
  }

  # Summarize columns with user fn
  # * normally if the per-row result from user fn is greater than length 1,
  #   mutate will throw an error unless result is wrapped in list
  # * want per-row results > length 1 to return list-col
  # * want results == 1 to return vector
  # * skip row if is.na
  # * not clear how to do this with across, instead do:
  #     1) use walk to iterate through vdj_cols
  #     2) use map to apply fn to each row in column
  #     3) within map check length of result and set (<<-) variable if length > 1
  fn <- purrr::as_mapper(fn, ...)

  purrr::walk2(vdj_cols, new_cols, ~ {
    LENGTH_ONE <- TRUE

    x <- res[[.x]]

    x <- purrr::map(x, ~ {
      if (all(is.na(.x))) {
        return(NA)
      }

      r <- fn(.x)

      if (length(r) > 1) {
        LENGTH_ONE <<- FALSE
      }

      r
    })

    if (LENGTH_ONE) {
      x <- unlist(x)
    }

    res <<- dplyr::mutate(res, !!sym(.y) := x)
  })

  # If chain provided remove temporary columns
  if (!is.null(chain)) {
    res <- dplyr::select(res, -all_of(vdj_cols))
  }

  # Re-nest vdj_cols
  res <- .nest_vdj(res, sep_cols = NULL, sep = sep)

  # Add results back to object
  if (return_df) {
    input <- res
  }

  res <- .add_meta(input, meta = res)

  res
}


#' Modify object meta.data
#'
#' @param input Single cell object or data.frame containing V(D)J data. If a
#' data.frame is provided, the cell barcodes should be stored as row names.
#' @param fn Function or formula to use for modifying the meta.data. If a
#' formula is provided, use .x to refer to the meta.data table.
#' @param ... Arguments to pass to the provided function
#' @return Object with mutated meta.data
#' @export
mutate_meta <- function(input, fn, ...) {

  if (!purrr::is_function(fn) && !is_formula(fn)) {
    stop("fn must be either a function or a formula")
  }

  meta <- .get_meta(input)

  if (purrr::is_formula(fn)) {
    fn <- as_mapper(fn, ...)
  }

  res <- fn(meta, ...)

  res <- .add_meta(input, meta = res)

  res
}

