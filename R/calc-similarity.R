#' Calculate repertoire overlap
#'
#' @param input Object containing V(D)J data. If a data.frame is provided, the
#' cell barcodes should be stored as row names.
#' @param cluster_col meta.data column containing cluster IDs to use for
#' calculating repertoire overlap
#' @param method Method to use for calculating similarity between clusters
#' @param clonotype_col meta.data column containing clonotype IDs to use for
#' calculating overlap
#' @param prefix Prefix to add to new columns
#' @param return_mat Return a matrix with similarity values. If set to
#' FALSE, results will be added to the input object.
#' @return Single cell object or data.frame with similarity values
#' @importFrom abdiv jaccard
#'
#' @examples
#' # Calculate repertoire overlap
#' res <- calc_similarity(
#'   vdj_so,
#'   clonotype_col = "clonotype_id",
#'   cluster_col = "orig.ident",
#'   method = abdiv::jaccard
#' )
#'
#' head(res@meta.data, 1)
#'
#' # Add a prefix to the new columns
#' # this is useful if multiple calculations are stored in the meta.data
#' res <- calc_similarity(
#'   vdj_sce,
#'   cluster_col = "orig.ident",
#'   prefix = "bcr_"
#' )
#'
#' head(res@colData, 1)
#'
#' # Return a matrix instead of adding the results to the input object
#' calc_similarity(
#'   vdj_so,
#'   cluster_col = "orig.ident",
#'   return_mat = TRUE
#' )
#'
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
    all_of(c(CELL_COL, clonotype_col, cluster_col))
  )

  if (dplyr::n_distinct(vdj[[cluster_col]]) < 2) {
    stop("cluster_col must contain at least two unique groups.")
  }

  vdj <- dplyr::group_by(
    vdj,
    !!!syms(c(cluster_col, clonotype_col))
  )

  vdj <- dplyr::summarize(
    vdj,
    n       = dplyr::n_distinct(!!sym(CELL_COL)),
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

  res <- purrr::map_dfr(combs, ~ {
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

  # Calculate self similarity
  res_s <- purrr::map_dfr(clsts, ~ {
    d <- pull(vdj, .x)

    tibble::tibble(
      Var1 = .x,
      Var2 = .x,
      sim  = method(d, d)
    )
  })

  # Add inverse combinations
  res_i <- dplyr::rename(res, Var1 = .data$Var2, Var2 = .data$Var1)

  res <- dplyr::bind_rows(res, res_i, res_s)
  res <- dplyr::mutate(res, Var1 = paste0(prefix, .data$Var1))

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


#' Plot repertoire overlap
#'
#' @param input Single cell object or data.frame containing V(D)J data. If a
#' data.frame is provided, the cell barcodes should be stored as row names.
#' @param cluster_col meta.data column containing cluster IDs to use for
#' calculating overlap
#' @param method Method to use for calculating similarity between clusters
#' @param clonotype_col meta.data column containing clonotype IDs to use for
#' calculating overlap
#' @param plot_colors Character vector containing colors for plotting
#' @param ... Additional arguments to pass to ggplot2, e.g. color, fill, size,
#' linetype, etc.
#' @return ggplot object
#' @importFrom abdiv jaccard
#'
#' @examples
#' # Plot repertoire overlap
#' plot_similarity(
#'   vdj_so,
#'   cluster_col = "orig.ident"
#' )
#'
#' # Specify method to use for calculating repertoire overlap
#' plot_similarity(
#'   vdj_sce,
#'   cluster_col = "orig.ident",
#'   method = abdiv::jaccard
#' )
#'
#' # Specify colors to use for heatmap
#' plot_similarity(
#'   vdj_so,
#'   cluster_col = "orig.ident",
#'   plot_color = c("white", "red")
#' )
#'
#' # Pass additional aesthetic parameters to ggplot2
#' plot_similarity(
#'   vdj_sce,
#'   cluster_col = "orig.ident",
#'   color = "black",
#'   size = 2
#' )
#'
#' @export
plot_similarity <- function(input, cluster_col, method = abdiv::jaccard, clonotype_col = "clonotype_id",
                            plot_colors = NULL, ...) {

  # Calculate similarity
  plt_dat <- calc_similarity(
    input         = input,
    cluster_col   = cluster_col,
    method        = method,
    clonotype_col = clonotype_col,
    prefix        = "",
    return_mat    = TRUE
  )

  sim_col <- as.character(substitute(method))
  sim_col <- dplyr::last(sim_col)

  var_lvls <- unique(c(rownames(plt_dat), colnames(plt_dat)))
  var_lvls <- sort(var_lvls)

  plt_dat <- tibble::as_tibble(plt_dat, rownames = "Var1")

  plt_dat <- tidyr::pivot_longer(
    plt_dat,
    cols      = -.data$Var1,
    names_to  = "Var2",
    values_to = sim_col
  )

  # Set Var levels
  plt_dat <- dplyr::mutate(
    plt_dat,
    Var1 = factor(.data$Var1, levels = rev(var_lvls)),
    Var2 = factor(.data$Var2, levels = var_lvls)
  )

  # Create heatmap
  res <- .create_heatmap(
    plt_dat,
    x     = "Var1",
    y     = "Var2",
    .fill = sim_col,
    clrs  = plot_colors,
    ...
  )

  res
}

