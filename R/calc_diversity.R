#' Calculate repertoire diversity
#'
#' @param input Single cell object or data.frame containing V(D)J data. If a
#' data.frame is provided, the cell barcodes should be stored as row names.
#' @param clonotype_col meta.data column containing clonotype IDs to use for
#' calculating diversity
#' @param cluster_col meta.data column containing cluster IDs to use for
#' calculating diversity. If cluster_col is omitted, diversity index will be
#' calculated for all clonotypes.
#' @param method Method to use for calculating diversity. A named list can also
#' be passed to use multiple methods. The names should specify names for the
#' output columns.
#' @param prefix Prefix to add to new columns
#' @param return_df Return results as a data.frame. If set to FALSE, results
#' will be added to the input object.
#' @return Single cell object or data.frame with diversity metrics
#' @export
calc_diversity <- function(input, clonotype_col = "cdr3_nt", cluster_col = NULL,
                           method = abdiv::simpson, prefix = "", return_df = FALSE) {

  if (length(method) > 1 && is.null(names(method))) {
    stop("Must include names if using a list of methods.")
  }

  if (length(method) == 1 && is.null(names(method))) {
    nm <- as.character(substitute(method))
    nm <- last(nm)

    method <- purrr::set_names(list(method), nm)
  }

  # Format input data
  meta <- .get_meta(input)
  vdj  <- dplyr::filter(meta, !is.na(!!sym(clonotype_col)))

  # Count clonotypes
  vdj_cols <- clonotype_col

  if (!is.null(cluster_col)) {
    vdj_cols <- c(cluster_col, vdj_cols)
    vdj   <- dplyr::group_by(vdj, !!sym(cluster_col))
  }

  vdj <- dplyr::group_by(vdj, !!sym(clonotype_col), .add = TRUE)

  vdj <- dplyr::summarize(
    vdj,
    .n      = dplyr::n_distinct(.data$.cell_id),
    .groups = "drop"
  )

  # Calculate diversity
  if (!is.null(cluster_col)) {
    vdj <- dplyr::group_by(vdj, !!sym(cluster_col))
  }

  if (length(method) == 1 && is.null(names(method))) {
    nm <- as.character(substitute(method))
    nm <- dplyr::last(nm)

    method        <- list(method)
    names(method) <- nm
  }

  div_cols      <- paste0(prefix, names(method))
  names(method) <- div_cols

  vdj <- purrr::imap_dfr(method, ~ {
    dplyr::mutate(
      vdj,
      diversity = .x(.data$.n),
      met       = .y
    )
  })

  vdj <- tidyr::pivot_wider(
    vdj,
    names_from  = "met",
    values_from = "diversity"
  )

  vdj <- dplyr::ungroup(vdj)
  vdj <- dplyr::select(vdj, all_of(c(vdj_cols, div_cols)))

  # Format results
  res <- dplyr::left_join(meta, vdj, by = vdj_cols)

  if (return_df) {
    input <- meta
  }

  res <- .add_meta(input, meta = res)

  res
}
