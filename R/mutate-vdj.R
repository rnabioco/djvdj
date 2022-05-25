#' Fetch V(D)J data from object
#'
#' Fetch per-chain V(D)J data from object. Within the object meta.data, each
#' row represents a single cell and can include information for multiple
#' chains. This function returns an unnested data.frame where each row
#' represents a single chain. This is useful for plotting per-chain metrics
#' such as CDR3 length or the number of insertions/deletions.
#'
#' @param input Single cell object or data.frame containing V(D)J data. If a
#' data.frame is provided, the cell barcodes should be stored as row names.
#' @param vdj_cols meta.data columns containing per-chain V(D)J data to unnest.
#' If NULL, V(D)J data are automatically selected by identifying columns that
#' have NAs in the same rows as clonotype_col.
#' @param clonotype_col meta.data column containing clonotype IDs. This column
#' is used to determine which columns have V(D)J data. If both clonotype_col
#' and vdj_cols are NULL, all columns are included.
#' @param filter_cells Remove cells that do not have V(D)J data, clonotype_col
#' must be provided to determine which cells to filter.
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
#'   vdj_cols = c("chains", "n_insertion")
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
fetch_vdj <- function(input, vdj_cols = NULL, clonotype_col = NULL, filter_cells = FALSE,
                      unnest = TRUE, sep = ";") {

  # Format input data
  meta <- .get_meta(input)

  if (is.null(sep)) {
    return(meta)
  }

  # Identify columns with V(D)J data
  col_list <- .get_vdj_cols(
    df_in     = meta,
    clone_col = clonotype_col,
    cols_in   = vdj_cols,
    sep       = sep
  )

  sep_cols <- col_list$sep

  if (purrr::is_empty(sep_cols)) {
    warning(
      "The separator '", sep, "' was not identified in any columns specified ",
      "by vdj_cols, the unmodified meta.data will be returned"
    )

    return(meta)
  }

  # Filter cells
  if (filter_cells) {
    if (is.null(clonotype_col)) {
      stop("clonotype_col must be provided to determine which cells to filter.")
    }

    meta <- dplyr::filter(meta, !is.na(!!sym(clonotype_col)))
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
#' @param vdj_cols meta.data columns containing V(D)J data to modify. If NULL,
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
#' head(res@meta.data, 3)
#'
#' # Calculate the total number of insertions + deletions for each chain
#' # we have to wrap our expression in list() since a value is returned for
#' # each chain
#' res <- mutate_vdj(
#'   vdj_sce,
#'   indels = list(n_insertion + n_deletion)
#' )
#'
#' head(res@colData, 3)
#'
#' # Create a new column showing the unique chains for each cell
#' res <- mutate_vdj(
#'   vdj_so,
#'   unique_chains = stringr::str_c(unique(chains), collapse = "_")
#' )
#'
#' head(res@meta.data, 3)
#'
#' # Determine which cells have both an IGK and IGL chain
#' res <- mutate_vdj(
#'   vdj_sce,
#'   both_light = all(c("IGK", "IGL") %in% chains)
#' )
#'
#' head(res@colData, 1)
#'
#' # Determine which cells have multiple light chains
#' res <- mutate_vdj(
#'   vdj_so,
#'   multi_light = sum(chains %in% c("IGK", "IGL")) > 1
#' )
#'
#' head(res@meta.data, 3)
#'
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
  if (!is.null(sep) && purrr::is_empty(sep_cols)) {
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
  res <- .nest_vdj(vdj, sep_cols = NULL, sep = sep)

  if (return_df) {
    input <- res
  }

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
#'   vdj_cols = c("n_deletion", "n_insertion")
#' )
#'
#' head(res@meta.data, 3)
#'
#' # Specifying a different summary function
#' # this calculates the median number of insertions and deletions for each
#' # cell
#' res <- summarize_vdj(
#'   vdj_sce,
#'   vdj_cols = c("n_deletion", "n_insertion"),
#'   fn = stats::median
#' )
#'
#' head(res@colData, 3)
#'
#' # Summarize values for a specific chain
#' res <- summarize_vdj(
#'   vdj_so,
#'   vdj_cols = c("n_deletion", "n_insertion"),
#'   chain = "IGK"
#' )
#'
#' head(res@meta.data, 3)
#'
#' # Specifying new names for summarized columns
#' # use {.col} to refer to the original column name
#' res <- summarize_vdj(
#'   vdj_sce,
#'   vdj_cols = c("n_deletion", "n_insertion"),
#'   fn = stats::median,
#'   col_names = "median_{.col}"
#' )
#'
#' head(res@colData, 1)
#'
#' # Return a data.frame instead of adding the results to the input object
#' res <- summarize_vdj(
#'   vdj_so,
#'   vdj_cols = c("n_deletion", "n_insertion"),
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
#'   vdj_cols = "chains",
#'   fn = ~ paste0(unique(.x), collapse = "_"),
#'   col_names = "unique_chains"
#' )
#'
#' head(res@colData, 3)
#'
#' # Creating an index column to use for filtering/plotting
#' # this creates a column indicating which cells have no insertions
#' # the V(D)J data can be filtered based on this new column
#' res <- summarize_vdj(
#'   vdj_so,
#'   vdj_cols = "n_insertion",
#'   fn = ~ all(.x == 0),
#'   col_names = "no_insertions"
#' )
#'
#' res <- filter_vdj(
#'   res,
#'   filt = no_insertions
#' )
#'
#' head(res@meta.data, 3)
#'
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
  # check that all vdj_cols contain per-chain data
  if (!is.null(chain)) {
    is_lst      <- purrr::map_lgl(res[, vdj_cols], is.list)
    not_lst     <- names(is_lst[!is_lst])
    filt_chains <- purrr::is_empty(not_lst)

    if (filt_chains) {
      prfx <- "FILT_"

      res <- .filter_chains(
        res,
        vdj_cols  = vdj_cols,
        chain     = chain,
        chain_col = chain_col,
        col_names = paste0(prfx, "{.col}"),
        empty_val = NA
      )

      # Add prefix to vdj_cols so temporary columns are usedvdj %>% hea
      vdj_cols <- paste0(prfx, vdj_cols)

      # Set col_names so prefix is removed from columns
      col_names <- gsub(
        "\\{.col\\}",
        paste0('{sub(\\"^', prfx, '\\", "", .col)}'),
        col_names
      )

    } else {
      not_lst <- paste0(not_lst, collapse = ", ")

      warning(
        not_lst, " does not contain per-chain V(D)J data, can only filter ",
        "based on chain if all vdj_cols contain per-chain data."
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
  #     1) iterate through vdj_cols
  #     2) apply fn to each row in column
  #     3) check length of fn result and set variable if length > 1
  #
  # With this current approach, user cannot refer to other columns in provided
  # fn, e.g. ~ .x[productive] does not work :(
  fn <- purrr::as_mapper(fn, ...)

  for (i in seq_along(vdj_cols)) {
    clmn       <- vdj_cols[i]
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

    if (length_one) {
      x <- unlist(x)
    }

    res[[new_clmn]] <- x
  }

  # If chain provided remove temporary columns
  if (!is.null(chain) && filt_chains) {
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
#' head(res@meta.data, 1)
#'
#' # Pass a purrr-style lambda
#' # this produces the same result as the previous example
#' res <- mutate_meta(
#'   tiny_so,
#'   ~ dplyr::mutate(.x, NEW = nCount_RNA + nFeature_RNA)
#' )
#'
#' head(res@meta.data, 1)
#'
#' # Modify multiple meta.data columns
#' res <- mutate_meta(
#'   tiny_sce,
#'   dplyr::mutate,
#'   NEW_1 = nCount_RNA + nFeature_RNA,
#'   NEW_2 = stringr::str_c(orig.ident, seurat_clusters)
#' )
#'
#' head(res@colData, 1)
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
#' head(res@meta.data, 1)
#'
#' # Perform grouped operations using dplyr
#' # multi-line commands can be passed using brackets, just refer to the
#' # meta.data with ".x"
#' # this calculates the mean number of features for each group in the
#' # orig.ident meta.data column
#' res <- mutate_meta(tiny_sce, ~ {
#'   y <- dplyr::group_by(.x, orig.ident)
#'   y <- dplyr::mutate(y, mean_genes = mean(nFeature_RNA))
#'   y
#' })
#'
#' head(res@colData, 1)
#'
#' @export
mutate_meta <- function(input, fn, ...) {

  if (!purrr::is_function(fn) && !purrr::is_formula(fn)) {
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

