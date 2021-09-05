#' Calculate clonotype abundance
#'
#' @param input Single cell object or data.frame containing V(D)J data. If a
#' data.frame is provided, the cell barcodes should be stored as row names.
#' @param clonotype_col meta.data column containing clonotype IDs
#' @param cluster_col meta.data column containing cluster IDs to use for
#' grouping cells when calculating clonotype abundance
#' @param prefix Prefix to add to new columns
#' @param return_df Return results as a data.frame. If set to FALSE, results
#' will be added to the input object.
#' @return Single cell object or data.frame with clonotype abundance metrics
#' @export
calc_abundance <- function(input, clonotype_col = "cdr3_nt", cluster_col = NULL,
                           prefix = "", return_df = FALSE) {

  # Format input data
  meta <- .get_meta(input)
  vdj  <- dplyr::filter(meta, !is.na(!!sym(clonotype_col)))

  vdj <- dplyr::select(
    vdj,
    .data$.cell_id, all_of(c(cluster_col, clonotype_col))
  )

  # Calculate clonotype abundance
  vdj <- .calc_abund(
    df_in     = vdj,
    cell_col  = ".cell_id",
    clone_col = clonotype_col,
    clust_col = cluster_col
  )

  new_cols <- c("freq", "pct")

  if (!is.null(cluster_col)) {
    new_cols <- c(new_cols, "shared")
  }

  new_cols <- purrr::set_names(
    paste0(".", new_cols),
    paste0(prefix, "clone_", new_cols)
  )

  vdj <- select(vdj, .data$.cell_id, !!!syms(new_cols))

  # Format results
  res <- dplyr::left_join(meta, vdj, by = ".cell_id")

  if (return_df) {
    input <- meta
  }

  res <- .add_meta(input, meta = res)

  res
}

#' Calculate clonotype abundance
#'
#' @param df_in Input data.frame
#' @param cell_col Column containing cell IDs
#' @param clone_col Column containing clonotype IDs
#' @param clust_col Column containing cluster IDs
#' @return data.frame
.calc_abund <- function(df_in, cell_col, clone_col, clust_col = NULL) {

  # Count number of cells in each group
  if (!is.null(clust_col)) {
    df_in <- dplyr::group_by(df_in, !!sym(clust_col))
  }

  df_in <- dplyr::mutate(
    df_in,
    .n_cells = dplyr::n_distinct(!!sym(cell_col))
  )

  # Calculate frequency
  res <- dplyr::group_by(df_in, !!sym(clone_col), .add = TRUE)

  res <- dplyr::mutate(
    res,
    .freq = dplyr::n_distinct(!!sym(cell_col)),
    .pct  = (.data$.freq / .data$.n_cells) * 100
  )

  # Identify shared clonotypes
  if (!is.null(clust_col)) {
    res <- dplyr::group_by(res, !!sym(clone_col))

    res <- dplyr::mutate(
      res,
      .shared = dplyr::n_distinct(!!sym(clust_col)) > 1
    )
  }

  res <- dplyr::ungroup(res)

  res
}
