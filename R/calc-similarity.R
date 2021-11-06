#' Calculate repertoire overlap
#'
#' @param input Object containing V(D)J data. If a data.frame is provided, the
#' cell barcodes should be stored as row names.
#' @param cluster_col meta.data column containing cluster IDs to use for
#' calculating overlap
#' @param method Method to use for calculating similarity between clusters
#' @param clonotype_col meta.data column containing clonotype IDs to use for
#' calculating overlap
#' @param prefix Prefix to add to new columns
#' @param return_mat Return a matrix with similarity values. If set to
#' FALSE, results will be added to the input object.
#' @return Single cell object or data.frame with similarity values
#' @export
calc_similarity <- function(input, cluster_col, method = abdiv::jaccard, clonotype_col = "clonotype_id",
                            prefix = NULL, return_mat = FALSE) {

  # If no prefix provided, use method name
  if (is.null(prefix)) {
    prefix <- as.character(substitute(method))
    prefix <- dplyr::last(prefix)
    prefix <- paste0(prefix, "_")
  }

  # Format input data
  meta <- .get_meta(input)
  vdj  <- dplyr::filter(meta, !is.na(!!sym(clonotype_col)))

  vdj <- dplyr::select(
    vdj,
    all_of(c(".cell_id", clonotype_col, cluster_col))
  )

  vdj <- dplyr::group_by(
    vdj,
    !!!syms(c(cluster_col, clonotype_col))
  )

  vdj <- dplyr::summarize(
    vdj,
    n       = dplyr::n_distinct(.data$.cell_id),
    .groups = "drop"
  )

  vdj <- tidyr::pivot_wider(
    vdj,
    names_from  = all_of(cluster_col),
    values_from = .data$n
  )

  # Calculate similarity index
  clsts <- colnames(vdj)
  clsts <- clsts[clsts != clonotype_col]

  vdj <- dplyr::mutate(vdj, across(
    all_of(clsts), ~ tidyr::replace_na(.x, 0)
  ))

  combs <- utils::combn(clsts, 2, simplify = FALSE)

  res <- map_dfr(combs, ~ {
    x <- pull(vdj, .x[1])
    y <- pull(vdj, .x[2])

    tibble::tibble(
      Var1 = .x[1],
      Var2 = .x[2],
      sim  = method(x, y)
    )
  })

  # Return matrix
  if (return_mat) {
    res <- tidyr::pivot_wider(
      res,
      names_from  = .data$Var1,
      values_from = .data$sim
    )

    res <- tibble::column_to_rownames(res, "Var2")
    res <- as.matrix(res)

    return(res)
  }

  # Add inverse combinations
  res_i <- dplyr::rename(res, Var1 = .data$Var2, Var2 = .data$Var1)

  clsts <- tibble::tibble(Var1 = clsts, Var2 = clsts, sim = 1)
  res   <- dplyr::bind_rows(res, res_i, clsts)
  res   <- dplyr::mutate(res, Var1 = paste0(prefix, .data$Var1))

  res <- tidyr::pivot_wider(
    res,
    names_from  = .data$Var1,
    values_from = .data$sim
  )

  # Add results to input
  j_cols <- purrr::set_names("Var2", cluster_col)

  res <- dplyr::left_join(meta, res, by = j_cols)

  # Format results
  # join will cause factor levels to be lost, add these back
  clst_col <- pull(meta, cluster_col)

  if (is.factor(clst_col)) {
    res <- dplyr::mutate(
      res,
      !!sym(cluster_col) := factor(!!sym(cluster_col), levels = levels(clst_col))
    )
  }

  res <- .add_meta(input, meta = res)

  res
}