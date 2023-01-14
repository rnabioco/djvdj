#' Fetch V(D)J data from object
#'
#' Fetch per-chain V(D)J data from object. Within the object meta.data, each
#' row represents a single cell and can include information for multiple
#' chains. This function can return a data.frame where each row
#' represents a single chain. This is useful for plotting per-chain metrics
#' such as CDR3 length or the number of insertions/deletions.
#'
#' @param input Single cell object or data.frame containing V(D)J data. If a
#' data.frame is provided, the cell barcodes should be stored as row names.
#' @param data_cols meta.data columns containing per-chain V(D)J data to unnest.
#' If NULL, V(D)J data are automatically selected by identifying columns that
#' have NAs in the same rows as clonotype_col.
#' @param clonotype_col meta.data column containing clonotype IDs. This column
#' is used to determine which columns have V(D)J data. If both clonotype_col
#' and data_cols are NULL, all columns are included.
#' @param filter_cells Remove cells that do not have V(D)J data, clonotype_col
#' must be provided to determine which cells to filter.
#' @param per_cell Return per-cell data instead of per-chain data.
#' @param unnest If FALSE, a nested data.frame is returned where each row
#' represents a cell and V(D)J data is stored as list-cols. If TRUE, columns
#' are unnested so each row represents a single chain.
#' @param sep Separator used for storing per cell V(D)J data. This is used to
#' identify columns containing per-chain data that can be unnested.
#' @return data.frame containing V(D)J data
#'
#' @examples
#' # Fetch per-chain V(D)J data
#' fetch_vdj(vdj_so)
#'
#' # To increase performance, specify which columns to return per-chain data,
#' # per-cell data will be returned for all other columns
#' fetch_vdj(
#'   vdj_sce,
#'   data_cols = c("chains", "reads")
#' )
#'
#' # Only include cells that have V(D)J data
#' # clonotype_col must be specified to identify cells with V(D)J data
#' fetch_vdj(
#'   vdj_so,
#'   filter_cells = TRUE,
#'   clonotype_col = "clonotype_id"
#' )
#'
#' @export
fetch_vdj <- function(input, data_cols = NULL, clonotype_col = NULL,
                      filter_cells = FALSE, per_cell = FALSE, unnest = TRUE,
                      sep = global$sep) {

  # Check that columns are present in object
  .check_obj_cols(input, data_cols, clonotype_col)

  # Check input classes
  .check_args(data_cols = list(len_one = FALSE, allow_null = TRUE))

  # Format input data
  meta <- .get_meta(input)

  # Filter cells
  if (filter_cells) {
    if (is.null(data_cols) && is.null(clonotype_col)) {
      cli::cli_abort(
        "When `filter_cells` is `TRUE`, `data_cols` or `clonotype_col` must be
         provided to determine which cells to filter"
      )
    }

    filt_cols <- clonotype_col %||% data_cols

    meta <- dplyr::filter(meta, dplyr::if_all(filt_cols, ~ !is.na(.x)))
  }

  # If NULL sep or per_cell, return meta.data
  if (is.null(sep) || per_cell) return(meta)

  # Identify columns with V(D)J data
  col_list <- .get_vdj_cols(
    df_in     = meta,
    clone_col = clonotype_col,
    cols_in   = data_cols,
    sep       = sep
  )

  sep_cols <- col_list$sep

  if (purrr::is_empty(sep_cols)) {
    cli::cli_warn("`sep` ({sep}) was not found in {.or {data_cols}}")

    return(meta)
  }

  # Unnest V(D)J data
  res <- .unnest_vdj(
    meta,
    sep_cols = sep_cols,
    sep      = sep,
    unnest   = unnest
  )

  res
}


#' Modify V(D)J data in object
#'
#' Modify per-chain V(D)J data for each cell. This function offers greater
#' flexibility than [summarize_vdj()], but is less user-friendly.
#'
#' @param input Single cell object or data.frame containing V(D)J data. If a
#' data.frame is provided, cell barcodes should be stored as row names.
#' @param ... Name-value pairs to use for creating or modifying per-chain V(D)J
#' meta.data, e.g. mean_umis = mean(umis).
#'
#' To allow for modification of per-chain V(D)J data, the data for each cell is
#' converted into a vector, e.g. 'IGH;IGK' is equivalent to c('IGH', 'IGK').
#' This allows R vector operations to be performed on the per-chain values.
#' Any operations that produce a result greater than length 1 need to be
#' returned as a list(), e.g. new_col = umis + 1 will return a new value for
#' each chain, to prevent an error this must be written as
#' new_col = list(umis + 1).
#'
#' @param clonotype_col meta.data column containing clonotype IDs. This is used
#' to identify columns containing V(D)J data.
#' @param data_cols meta.data columns containing V(D)J data to modify. If NULL,
#' data are automatically selected by identifying columns that have NAs in the
#' same rows as clonotype_col.
#' @param return_df Return results as a data.frame. If FALSE, results will be
#' added to the input object.
#' @param sep Separator used for storing per cell V(D)J data. If NULL, columns
#' containing V(D)J data will not be converted to vectors for filtering.
#' @return Object with modified meta.data
#'
#' @examples
#' # Calculate mean reads and UMIs per cell
#' res <- mutate_vdj(
#'   vdj_so,
#'   mean_umis  = mean(umis),
#'   mean_reads = mean(reads)
#' )
#'
#' head(slot(res, "meta.data"), 3)
#'
#' # Calculate the total number of insertions + deletions for each chain
#' # we have to wrap our expression in list() since a value is returned for
#' # each chain
#' res <- mutate_vdj(
#'   vdj_sce,
#'   indels = list(all_ins + all_del)
#' )
#'
#' head(slot(res, "colData"), 3)
#'
#' # Create a new column showing the unique chains for each cell
#' res <- mutate_vdj(
#'   vdj_so,
#'   unique_chains = stringr::str_c(unique(chains), collapse = "_")
#' )
#'
#' head(slot(res, "meta.data"), 3)
#'
#' # Determine which cells have both an IGK and IGL chain
#' res <- mutate_vdj(
#'   vdj_sce,
#'   both_light = all(c("IGK", "IGL") %in% chains)
#' )
#'
#' head(slot(res, "colData"), 1)
#'
#' # Determine which cells have multiple light chains
#' res <- mutate_vdj(
#'   vdj_so,
#'   multi_light = sum(chains %in% c("IGK", "IGL")) > 1
#' )
#'
#' head(slot(res, "meta.data"), 3)
#'
#' @export
mutate_vdj <- function(input, ..., clonotype_col = global$clonotype_col,
                       data_cols = NULL, return_df = FALSE, sep = global$sep) {

  # Check that columns are present in object
  .check_obj_cols(input, data_cols, clonotype_col)

  # Check input classes
  .check_args(data_cols = list(len_one = FALSE, allow_null = TRUE))

  # Format input data
  meta <- .get_meta(input)
  vdj  <- meta

  # Identify columns with V(D)J data
  col_list <- .get_vdj_cols(
    df_in     = vdj,
    clone_col = clonotype_col,
    cols_in   = data_cols,
    sep       = sep
  )

  vdj_cols <- col_list$vdj
  sep_cols <- col_list$sep

  # Allow sep to be NULL so user can skip .unnest_vdj
  if (!is.null(sep) && purrr::is_empty(sep_cols)) {
    cli::cli_warn("`sep` ('{sep}') is not present in the data")
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
  res <- .nest_vdj(vdj, sep_cols = NULL, sep = sep)

  if (return_df) input <- res

  res <- .add_meta(input, meta = res)

  res
}


#' Summarize V(D)J data for each cell
#'
#' Summarize per-chain values for each cell using a function or purrr-style
#' lambda. This is useful for plotting or filtering cells based on the V(D)J
#' meta.data.
#'
#' @param input Single cell object or data.frame containing V(D)J data. If a
#' data.frame is provided, the cell barcodes should be stored as row names.
#' @param data_cols meta.data column(s) containing V(D)J data to summarize for
#' each cell
#' @param fn Function to apply to each selected column, possible values can be
#' either a function, e.g. mean, or a purrr-style lambda, e.g. ~ mean(.x,
#' na.rm = TRUE). If NULL, the mean will be calculated for numeric values,
#' non-numeric columns will be combined into a single string.
#' @param ... Additional arguments to pass to fn
#' @param chain Chain to use for summarizing V(D)J data
#' @param chain_col meta.data column(s) containing chains for each cell
#' @param col_names A glue specification that describes how to name the output
#' columns, use \{.col\} to refer to the original column name. If col_names is
#' NULL, the original column names will be used.
#' @param return_df Return results as a data.frame. If FALSE, results will be
#' added to the input object.
#' @param sep Separator used for storing per cell V(D)J data
#' @return Object containing V(D)J data summarized for each cell
#' @importFrom glue glue
#'
#' @examples
#' # Summarize numeric columns
#' # by default the mean will be calculated for numeric columns
#' res <- summarize_vdj(
#'   vdj_so,
#'   data_cols = c("all_del", "all_ins")
#' )
#'
#' head(slot(res, "meta.data"), 3)
#'
#' # Specifying a different summary function
#' # this calculates the median number of insertions and deletions for each
#' # cell
#' res <- summarize_vdj(
#'   vdj_sce,
#'   data_cols = c("all_del", "all_ins"),
#'   fn = stats::median
#' )
#'
#' head(slot(res, "colData"), 3)
#'
#' # Summarize values for a specific chain
#' res <- summarize_vdj(
#'   vdj_so,
#'   data_cols = c("all_del", "all_ins"),
#'   chain = "IGK"
#' )
#'
#' head(slot(res, "meta.data"), 3)
#'
#' # Specifying new names for summarized columns
#' # use {.col} to refer to the original column name
#' res <- summarize_vdj(
#'   vdj_sce,
#'   data_cols = c("all_del", "all_ins"),
#'   fn = stats::median,
#'   col_names = "median_{.col}"
#' )
#'
#' head(slot(res, "colData"), 1)
#'
#' # Return a data.frame instead of adding the results to the input object
#' res <- summarize_vdj(
#'   vdj_so,
#'   data_cols = c("all_del", "all_ins"),
#'   return_df = TRUE
#' )
#'
#' head(res, 1)
#'
#' # Using a lambda function to summarize values
#' # use '.x' to refer to values in the column
#' # this creates a new column showing the unique chains for each cell
#' res <- summarize_vdj(
#'   vdj_sce,
#'   data_cols = "chains",
#'   fn = ~ paste0(unique(.x), collapse = "_"),
#'   col_names = "unique_chains"
#' )
#'
#' head(slot(res, "colData"), 3)
#'
#' # Creating an index column to use for filtering/plotting
#' # this creates a column indicating which cells have no insertions
#' # the V(D)J data can be filtered based on this new column
#' res <- summarize_vdj(
#'   vdj_so,
#'   data_cols = "all_ins",
#'   fn = ~ all(.x == 0),
#'   col_names = "no_insertions"
#' )
#'
#' res <- filter_vdj(
#'   res,
#'   filt = no_insertions
#' )
#'
#' head(slot(res, "meta.data"), 3)
#'
#' @export
summarize_vdj <- function(input, data_cols, fn = NULL, ..., chain = NULL,
                          chain_col = global$chain_col, col_names = "{.col}",
                          return_df = FALSE, sep = global$sep) {

  # Check that columns are present in object
  .check_obj_cols(input, data_cols, chain = chain, chain_col = chain_col)

  # Check input classes
  .check_args()

  # Names of new columns
  new_cols <- gsub("\\{.col\\}", "{data_cols}", col_names)
  new_cols <- glue::glue(new_cols)

  # Fetch V(D)J data
  # * this is a key performance bottleneck as the length of data_cols increases
  fetch_cols <- data_cols

  if (!is.null(chain)) fetch_cols <- c(data_cols, chain_col)

  res <- fetch_vdj(
    input,
    clonotype_col = NULL,
    data_cols     = fetch_cols,
    sep           = sep,
    unnest        = FALSE
  )

  # Set default fn
  # * default is mean if all vdj_col are numeric
  # * if !is.null(chain), collapse data_cols into string
  # * if is.null(chain), return input unchanged
  if (is.null(fn)) {
    is_num <- purrr::map_lgl(
      res[, data_cols],
      ~ all(purrr::map_lgl(.x, is.numeric))
    )

    is_num <- all(is_num)

    fn <- mean

    if (!is_num) {
      if (is.null(chain)) return(input)

      fn <- ~ paste0(.x, collapse = sep)
    }
  }

  # Filter chains
  # check that all data_cols contain per-chain data
  if (!is.null(chain)) {
    is_lst      <- purrr::map_lgl(res[, data_cols], is.list)
    not_lst     <- names(is_lst[!is_lst])
    filt_chains <- purrr::is_empty(not_lst)

    if (filt_chains) {
      prfx <- "FILT_"

      res <- .filter_chains(
        res,
        data_cols = data_cols,
        chain     = chain,
        chain_col = chain_col,
        col_names = paste0(prfx, "{.col}"),
        empty_val = NA
      )

      # Add prefix to data_cols so temporary columns are used
      data_cols <- paste0(prfx, data_cols)

      # Set col_names so prefix is removed from columns
      col_names <- gsub(
        "\\{.col\\}",
        paste0('{sub(\\"^', prfx, '\\", "", .col)}'),
        col_names
      )

    } else {
      cli::cli_warn(
        "Some columns do not contain per-chain V(D)J data, can only filter
         based on `chain` if all `data_cols` contain per-chain data: {not_lst}"
      )
    }
  }

  # Summarize columns with user fn
  # * normally if the per-row result from user fn is greater than length 1,
  #   mutate will throw an error unless result is wrapped in list
  # * want per-row results > length 1 to return list-col
  # * want results == 1 to return vector
  # * skip row if is.na
  # * use for loops instead of map to avoid '<<-' operator
  # * not clear how to do this with across, instead do:
  #     1) iterate through data_cols
  #     2) apply fn to each row in column
  #     3) check length of fn result and set variable if length > 1
  #
  # With this current approach, user cannot refer to other columns in provided
  # fn, e.g. ~ .x[productive] does not work :(
  fn <- purrr::as_mapper(fn, ...)

  for (i in seq_along(data_cols)) {
    clmn       <- data_cols[i]
    new_clmn   <- new_cols[i]
    length_one <- TRUE

    x <- res[[clmn]]

    for (j in seq_along(x)) {
      val <- x[[j]]

      if (all(is.na(val))) {
        r <- NA

      } else {
        r <- fn(val)

        if (purrr::is_empty(r)) {
          r <- NA

        } else if (length(r) > 1) {
          length_one <- FALSE
        }
      }

      x[[j]] <- r
    }

    if (length_one) x <- unlist(x)

    res[[new_clmn]] <- x
  }

  # If chain provided remove temporary columns
  if (!is.null(chain) && filt_chains) {
    res <- dplyr::select(res, -all_of(data_cols))
  }

  # Re-nest data_cols
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
#'
#' @param fn Function to use for modifying object meta.data. This can be either
#' a function, e.g. mean, or a purrr-style lambda, e.g. ~ mean(.x,
#' na.rm = TRUE) where ".x" refers to the meta.data table.
#' @param ... Additional arguments to pass to the provided function
#' @return Object with mutated meta.data
#'
#' @examples
#' # Sum two meta.data columns
#' # all additional arguments provided to mutate_meta() are passed directly to
#' # the function (in this case, dplyr::mutate())
#' res <- mutate_meta(
#'   tiny_so,
#'   dplyr::mutate,
#'   NEW = nCount_RNA + nFeature_RNA
#' )
#'
#' head(slot(res, "meta.data"), 1)
#'
#' # Pass a purrr-style lambda
#' # this produces the same result as the previous example
#' res <- mutate_meta(
#'   tiny_so,
#'   ~ dplyr::mutate(.x, NEW = nCount_RNA + nFeature_RNA)
#' )
#'
#' head(slot(res, "meta.data"), 1)
#'
#' # Modify multiple meta.data columns
#' res <- mutate_meta(
#'   tiny_sce,
#'   dplyr::mutate,
#'   NEW_1 = nCount_RNA + nFeature_RNA,
#'   NEW_2 = stringr::str_c(orig.ident, seurat_clusters)
#' )
#'
#' head(slot(res, "colData"), 1)
#'
#' # Remove meta.data columns
#' # any function can be passed to mutate_meta(), in this example
#' # dplyr::select() is used to remove columns
#' res <- mutate_meta(
#'   tiny_so,
#'   dplyr::select,
#'   -UMAP_1
#' )
#'
#' head(slot(res, "meta.data"), 1)
#'
#' # Perform grouped operations using dplyr
#' # multi-line commands can be passed using brackets, just refer to the
#' # meta.data with '.x'
#' # this calculates the mean number of features for each group in the
#' # orig.ident meta.data column
#' res <- mutate_meta(tiny_sce, ~ {
#'   y <- dplyr::group_by(.x, orig.ident)
#'   y <- dplyr::mutate(y, mean_genes = mean(nFeature_RNA))
#'   y
#' })
#'
#' head(slot(res, "colData"), 1)
#'
#' @export
mutate_meta <- function(input, fn, ...) {

  psbl_fns <- c("function", "formula")

  if (!purrr::is_function(fn) && !purrr::is_formula(fn)) {
    cli::cli_abort("`fn` must be a {.or {psbl_fns}}")
  }

  meta <- .get_meta(input)

  if (purrr::is_formula(fn)) {
    fn <- as_mapper(fn, ...)
  }

  res <- fn(meta, ...)

  res <- .add_meta(input, meta = res)

  res
}

