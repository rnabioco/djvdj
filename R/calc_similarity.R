#' Calculate repertoire overlap
#'
#' @export
calc_similarity <- function(input, ...) {
  UseMethod("calc_similarity", input)
}

#' @rdname calc_similarity
#' @param input Single cell object or data.frame containing V(D)J data. If a
#' data.frame is provided, the cell barcodes should be stored as row names.
#' @param clonotype_col meta.data column containing clonotype IDs to use for
#' calculating overlap
#' @param cluster_col meta.data column containing cluster IDs to use for
#' calculating overlap
#' @param method Method to use for calculating similarity between clusters
#' @param prefix Prefix to add to new columns
#' @param return_mat Return a matrix with similarity values. If set to
#' FALSE, results will be added to the input object.
#' @param ... Arguments passed to other methods
#' @return Single cell object or data.frame with similarity values
#' @export
calc_similarity.default <- function(input, clonotype_col = "cdr3_nt", cluster_col, method = abdiv::jaccard,
                                    prefix = NULL, return_mat = FALSE, ...) {

  if (is.null(prefix)) {
    prefix <- as.character(substitute(method))
    prefix <- dplyr::last(prefix)
    prefix <- paste0(prefix, "_")
  }

  # Format meta.data
  input   <- tibble::as_tibble(input, rownames = ".cell_id")
  meta_df <- dplyr::filter(input, !is.na(!!sym(clonotype_col)))

  meta_df <- dplyr::select(
    meta_df,
    all_of(c(".cell_id", clonotype_col, cluster_col))
  )

  vdj_df <- dplyr::group_by(
    meta_df,
    !!!syms(c(cluster_col, clonotype_col))
  )

  vdj_df <- dplyr::summarize(
    vdj_df,
    n       = dplyr::n_distinct(.data$.cell_id),
    .groups = "drop"
  )

  vdj_df <- tidyr::pivot_wider(
    vdj_df,
    names_from  = all_of(cluster_col),
    values_from = .data$n
  )

  # Calculate similarity index
  clusts <- colnames(vdj_df)
  clusts <- clusts[clusts != clonotype_col]

  vdj_df <- dplyr::mutate(vdj_df, across(
    all_of(clusts), ~ tidyr::replace_na(.x, 0)
  ))

  combs <- utils::combn(clusts, 2, simplify = FALSE)

  res <- map_dfr(combs, ~ {
    ins <- paste0("vdj_df$", .x)

    x <- pull(vdj_df, .x[1])
    y <- pull(vdj_df, .x[2])

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

  clusts <- tibble::tibble(Var1 = clusts, Var2 = clusts, sim = 1)
  res    <- dplyr::bind_rows(res, res_i, clusts)
  res    <- dplyr::mutate(res, Var1 = paste0(prefix, .data$Var1))

  res <- tidyr::pivot_wider(
    res,
    names_from  = .data$Var1,
    values_from = .data$sim
  )

  # Add results to input
  j_cols <- purrr::set_names("Var2", cluster_col)

  res <- dplyr::left_join(input, res, by = j_cols)

  # Format results
  # join will cause factor levels to be lost, add these back
  clsts <- pull(input, cluster_col)

  if (is.factor(clsts)) {
    res <- dplyr::mutate(
      res,
      !!sym(cluster_col) := factor(!!sym(cluster_col), levels = levels(clsts))
    )
  }

  res <- tibble::column_to_rownames(res, ".cell_id")

  res
}

#' @rdname calc_similarity
#' @export
calc_similarity.Seurat <- function(input, clonotype_col = "cdr3_nt", cluster_col, method = abdiv::jaccard,
                                   prefix = NULL, return_mat = FALSE, ...) {

  if (is.null(prefix)) {
    prefix <- as.character(substitute(method))
    prefix <- dplyr::last(prefix)
    prefix <- paste0(prefix, "_")
  }

  res <- calc_similarity(
    input         = input@meta.data,
    clonotype_col = clonotype_col,
    cluster_col   = cluster_col,
    method        = method,
    prefix        = prefix,
    return_mat    = return_mat
  )

  if (!return_mat) {
    res <- Seurat::AddMetaData(input, metadata = res)
  }

  res
}
