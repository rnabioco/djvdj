#' Calculate gene usage
#'
#' Gene usage is calculated as the percent of cells that express each V(D)J
#' gene present in gene_cols. Cells that lack V(D)J data and have NAs present
#' in gene_cols are excluded from this calculation.
#'
#' @export
calc_usage <- function(input, ...) {
  UseMethod("calc_usage", input)
}

#' @rdname calc_usage
#' @param input Seurat object containing V(D)J data
#' @param gene_cols meta.data column containing genes used for each clonotype
#' @param cluster_col meta.data column containing cell clusters to use when
#' calculating gene usage
#' @param chain Chain to use for calculating gene usage. Set to NULL to include
#' all chains.
#' @param chain_col meta.data column containing chains for each cell
#' @param sep Separator to use for expanding gene_cols
#' @return data.frame containing gene usage summary
#' @export
calc_usage.default <- function(input, gene_cols, cluster_col = NULL, chain = NULL,
                               chain_col = "chains", sep = ";") {

  # data.frame to calculate usage
  sep_cols <- c(gene_cols, chain_col)
  vdj_cols   <- c(".cell_id", cluster_col, sep_cols)

  input   <- tibble::as_tibble(input, rownames = ".cell_id")
  meta_df <- dplyr::select(input, all_of(vdj_cols))

  meta_df <- dplyr::filter(meta_df, across(
    all_of(gene_cols),
    ~ !is.na(.x)
  ))

  res <- dplyr::mutate(meta_df, across(
    all_of(sep_cols),
    ~ strsplit(as.character(.x), sep)
  ))

  # Filter chains
  if (!is.null(chain)) {
    if (is.null(chain_col)) {
      stop("Must specify chain_col.")
    }

    res <- dplyr::mutate(res, across(all_of(gene_cols), ~ {
      g_col <- .x

      purrr::map2(g_col, !!sym(chain_col), ~ {
        .x <- dplyr::if_else(
          any(.y %in% chain),
          list(.x[.y %in% chain]),
          list("None")
        )

        unlist(.x)
      })
    }))

    res <- dplyr::select(res, -all_of(chain_col))
  }

  res <- tidyr::unnest(res, cols = all_of(gene_cols))
  res <- dplyr::distinct(res)

  # Count genes used
  res <- dplyr::group_by(res, !!!syms(gene_cols))

  if (!is.null(cluster_col)) {
    res <- dplyr::group_by(res, !!sym(cluster_col), .add = TRUE)
  }

  res <- dplyr::summarize(
    res,
    n_cells = n_distinct(meta_df$.cell_id),
    freq    = dplyr::n_distinct(.data$.cell_id),
    .groups = "drop"
  )

  # Calculate percentage used
  if (!is.null(cluster_col)) {
    clusts       <- pull(meta_df, cluster_col)
    clust_counts <- table(clusts)
    clusts       <- unique(clusts)

    res <- tidyr::pivot_wider(
      res,
      names_from  = all_of(cluster_col),
      values_from = .data$freq,
      values_fill = 0
    )

    res <- tidyr::pivot_longer(
      res,
      cols      = all_of(clusts),
      names_to  = cluster_col,
      values_to = "freq"
    )

    res <- dplyr::mutate(  # why doesn't .before work here??
      res,
      n_cells = as.numeric(clust_counts[!!sym(cluster_col)])
    )

    res <- dplyr::relocate(res, .data$n_cells, .before = .data$freq)
  }

  res <- dplyr::mutate(res, pct = (.data$freq / .data$n_cells) * 100)

  res
}

#' @rdname calc_usage
#' @export
calc_usage.Seurat <- function(input, gene_cols, cluster_col = NULL, chain = NULL,
                              chain_col = "chains", sep = ";") {

  res <- calc_usage(
    input       = input@meta.data,
    gene_cols   = gene_cols,
    cluster_col = cluster_col,
    chain       = chain,
    chain_col   = chain_col,
    sep         = sep
  )

  res
}


