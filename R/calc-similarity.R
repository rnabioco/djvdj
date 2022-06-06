#' Calculate repertoire overlap
#'
#' @param input Object containing V(D)J data. If a data.frame is provided, the
#' cell barcodes should be stored as row names.
#' @param cluster_col meta.data column containing cluster IDs to use for
#' calculating repertoire overlap
#' @param method Method to use for comparing clusters, possible values are:
#'
#' - 'mds', perform multidimensional scaling
#' - 'count', count the number of clonotypes overlapping between each cluster
#' - A function that takes two numeric vectors containing counts for each
#' unique value in the column provided to the data_col column, e.g.
#' abdiv::jaccard()
#'
#' @param data_col meta.data column containing clonotype IDs to use for
#' calculating overlap
#' @param prefix Prefix to add to new columns
#' @param return_mat Return a matrix with similarity values. If set to
#' FALSE, results will be added to the input object.
#' @param sep Separator used for storing per-chain V(D)J data for each cell
#' @return Single cell object or data.frame with similarity values
#' @importFrom abdiv jaccard
#'
#' @examples
#' # Calculate repertoire overlap
#' res <- calc_similarity(
#'   vdj_so,
#'   data_col = "clonotype_id",
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
calc_similarity <- function(input, data_col, cluster_col, method = abdiv::jaccard,
                            chain = NULL, chain_col = "chains", prefix = NULL,
                            return_mat = FALSE, sep = ";") {

  # Check inputs
  is_mds    <- identical(method, "mds")
  is_counts <- identical(method, "count")

  if (!is.function(method) && !is_mds && !is_counts) {
    stop(
      "method must be 'mds' or a function to use for comparing ",
      "sample repertoires, e.g. method = abdiv::jaccard."
    )
  }

  if (is_mds) method <- abdiv::jaccard

  if (is_counts) {
    method <- function(x, y) length(x[x > 0 & y > 0])

    if (is.null(prefix)) prefix <- "count_"
  }

  # If no prefix provided, use method name
  if (is.null(prefix) && !is_mds) {
    prefix <- as.character(substitute(method))
    prefix <- dplyr::last(prefix)
    prefix <- paste0(prefix, "_")
  }

  # Format input data
  # filter chains if provided
  if (is.null(chain)) {
    meta <- .get_meta(input)
    vdj  <- dplyr::filter(meta, !is.na(!!sym(data_col)))

  } else {
    vdj <- fetch_vdj(
      input,
      vdj_cols      = c(data_col, chain_col),
      clonotype_col = data_col,
      filter_cells  = TRUE,
      unnest        = FALSE
    )

    vdj <- .filter_chains(
      vdj,
      vdj_cols   = data_col,
      chain      = chain,
      chain_col  = chain_col,
      col_names  = "{.col}",
      allow_dups = FALSE
    )

    vdj <- dplyr::filter(vdj, !is.na(!!sym(data_col)))
  }

  vdj <- dplyr::select(
    vdj,
    all_of(c(CELL_COL, data_col, cluster_col))
  )

  if (dplyr::n_distinct(vdj[[cluster_col]]) < 2) {
    stop("cluster_col must contain at least two unique groups.")
  }

  vdj <- dplyr::group_by(
    vdj,
    !!!syms(c(cluster_col, data_col))
  )

  vdj <- dplyr::summarize(
    vdj,
    n       = dplyr::n_distinct(!!sym(CELL_COL)),
    .groups = "drop"
  )

  clsts <- unique(vdj[[cluster_col]])

  vdj <- tidyr::pivot_wider(
    vdj,
    names_from  = all_of(cluster_col),
    values_from = .data$n,
    values_fill = 0
  )

  # Calculate similarity index
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

  # Calculate self similarity
  res_s <- purrr::map_dfr(clsts, ~ {
    d <- pull(vdj, .x)

    tibble::tibble(
      Var1 = .x,
      Var2 = .x,
      sim  = method(d, d)
    )
  })

  # Combine with inverse combinations
  res_i <- dplyr::rename(res, Var1 = .data$Var2, Var2 = .data$Var1)

  res <- dplyr::bind_rows(res, res_i, res_s)

  # Format data.frame
  clmns <- sort(unique(res$Var2))

  res <- dplyr::arrange(res, Var2)

  res <- tidyr::pivot_wider(
    res,
    names_from  = .data$Var1,
    values_from = .data$sim
  )

  res <- dplyr::select(res, !!sym(cluster_col) := .data$Var2, all_of(clmns))

  if (is_mds || return_mat) {
    res <- tibble::column_to_rownames(res, cluster_col)
    res <- as.matrix(res)
  }

  # Calculate MDS
  if (is_mds) {
    clmns <- c("MDS_1", "MDS_2")

    mds_fn <- purrr::quietly(MASS::isoMDS)
    res    <- mds_fn(as.dist(res))
    res    <- res$result$points

    colnames(res) <- clmns

    res <- tibble::as_tibble(res, rownames = cluster_col)

  # Return matrix if method not mds
  } else if (return_mat) return(res)

  # Return data.frame if method mds
  if (return_mat) return(res)

  # Add column prefixes
  res <- dplyr::rename_with(
    res,
    ~ paste0(prefix, .x),
    all_of(clmns)
  )

  # Add results to input
  res <- dplyr::left_join(meta, res, by = cluster_col)

  # Format results
  # join will cause factor levels to be lost, add these back
  clst_col <- pull(meta, cluster_col)

  if (is.factor(clst_col)) {
    res <- dplyr::mutate(
      res,
      !!sym(cluster_col) := factor(
        !!sym(cluster_col),
        levels = levels(clst_col)
      )
    )
  }

  res <- .add_meta(input, meta = res)

  res
}


#' Plot repertoire overlap
#'
#' @param input Single cell object or data.frame containing V(D)J data. If a
#' data.frame is provided, the cell barcodes should be stored as row names.
#' @param data_col meta.data column containing clonotype IDs to use for
#' calculating overlap
#' @param cluster_col meta.data column containing cluster IDs to use for
#' calculating overlap
#' @param group_col meta.data column to use for grouping cluster IDs present in
#' cluster_col. This is useful when there are multiple replicates or patients
#' for each treatment condition.
#' @param method Method to use for comparing clusters, possible values are:
#'
#' - 'mds', perform multidimensional scaling, this will generate a scatter plot
#' - 'circos', create circos plot implemented with circlize::chordDiagram()
#' - 'count', count the number of clonotypes overlapping between each cluster,
#' this will generate a heatmap.
#' - A function that takes two numeric vectors containing counts for each
#' clonotype in the object, such as the beta diversity functions provided by
#' the abdiv package. This will generate a heatmap.
#'
#' @param chain Chain to use for calculating gene usage. Set to NULL to include
#' all chains.
#' @param chain_col meta.data column containing chains for each cell
#' @param plot_colors Character vector containing colors for plotting
#' @param plot_lvls Levels to use for ordering clusters
#' @param include_upper_triangle If FALSE, upper triangle for heatmap will not
#' be shown and rows/columns will not be clustered.
#' @param include_diagonal If FALSE, diagonal for heatmap will not be shown and
#' rows/columns will not be clustered.
#' @param sep Separator used for storing per-chain V(D)J data for each cell
#' @param ... Additional arguments to pass to plotting function
#' @importFrom abdiv jaccard
#' @importFrom ComplexHeatmap Heatmap
#' @importFrom circlize chordDiagram
#' @importFrom grid gpar
#' @seealso [calc_similarity()], [circlize::chordDiagram()], [ComplexHeatmap::Heatmap()],
#' [ggplot2::geom_point()]
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
plot_similarity <- function(input, data_col, cluster_col, group_col = NULL,
                            method = abdiv::jaccard, chain = NULL,
                            chain_col = "chains", plot_colors = NULL,
                            plot_lvls = names(plot_colors),
                            include_upper_triangle = TRUE,
                            include_diagonal = TRUE, sep = ";", ...) {

  # Check inputs
  is_mds  <- identical(method, "mds")
  is_circ <- identical(method, "circos")

  if (is_circ) method <- "count"

  .chk_group_cols(cluster_col, group_col)

  # Calculate similarity
  plt_dat <- calc_similarity(
    input       = input,
    data_col    = data_col,
    cluster_col = cluster_col,
    method      = method,
    chain       = chain,
    chain_col   = chain_col,
    prefix      = "",
    return_mat  = TRUE,
    sep         = sep
  )

  # Create circos plot
  plt_args <- list(...)

  if (is_circ) {
    grps <- NULL

    if (!is.null(group_col)) {
      meta <- .get_meta(input)
      grps <- dplyr::distinct(meta, !!!syms(c(cluster_col, group_col)))
      grps <- purrr::set_names(grps[[group_col]], grps[[cluster_col]])
    }

    .create_circos(
      plt_dat,
      clrs = plot_colors,
      lvls = plot_lvls,
      grps = grps
    )

    return(invisible())
  }

  # Create MDS plot
  if (is_mds) {
    res <- plot_features(
      plt_dat,
      x = "MDS_1",
      y = "MDS_2",
      feature     = cluster_col,
      plot_colors = plot_colors,
      plot_lvls   = plot_lvls,
      ...
    ) +
      ggrepel::geom_text_repel(ggplot2::aes(label = !!sym(cluster_col)))

    return(res)
  }

  # Similarity column
  sim_col <- as.character(substitute(method))
  sim_col <- dplyr::last(sim_col)




  # Plot colors and levels
  plt_args$col <- plot_colors %||% c("white", "#6A51A3")

  if (is.null(plot_lvls)) {
    plot_lvls <- unique(c(rownames(plt_dat), colnames(plt_dat)))
    plot_lvls <- sort(plot_lvls)

  } else {
    plt_args$cluster_rows    <- FALSE
    plt_args$cluster_columns <- FALSE
  }

  plt_dat <- plt_dat[plot_lvls, plot_lvls]

  # Remove upper triangle and/or diagonal
  if (!include_upper_triangle || !include_diagonal) {
    lvls_key <- purrr::set_names(plot_lvls)

    lvls_key <- purrr::imap(lvls_key, ~ {
      idx <- grep(.x, plot_lvls)

      v <- c()

      if (!include_upper_triangle) {
        v <- plot_lvls[idx:length(plot_lvls)]
        v <- v[v != .y]
      }

      if (!include_diagonal) {
        v <- c(v, .y)
      }

      v
    })

    for (i in seq_along(lvls_key)) {
      plt_dat[names(lvls_key[i]), lvls_key[[i]]] <- NA
    }

    # Remove rows/columns with all NAs
    na_idx  <- is.na(plt_dat)
    r_idx   <- rowSums(na_idx) != ncol(plt_dat)
    c_idx   <- colSums(na_idx) != nrow(plt_dat)
    plt_dat <- plt_dat[r_idx, c_idx]

    # Set plot arguments
    plt_args$cluster_rows    <- FALSE
    plt_args$cluster_columns <- FALSE
    plt_args$row_names_side  <- plt_args$row_names_side %||% "left"
    plt_args$na_col          <- plt_args$na_col %||% NA
  }

  # Set final heatmap parameters
  # use computed similarities for clustering, so use a distance function that
  # just converts a matrix to a dist object
  # can't directly use as.dist as the function since it accepts too many
  # arguments
  dist_fn <- function(input) stats::as.dist(input)

  lgd_params <- list(
    title_gp      = grid::gpar(fontface = "plain"),
    legend_height = ggplot2::unit(80, "pt"),
    title         = sim_col
  )

  plt_args$matrix                      <- plt_dat
  plt_args$heatmap_legend_param        <- plt_args$heatmap_legend_param %||% lgd_params
  plt_args$clustering_distance_rows    <- plt_args$clustering_distance_rows %||% dist_fn
  plt_args$clustering_distance_columns <- plt_args$clustering_distance_columns %||% dist_fn

  # Create heatmap
  res <- purrr::lift_dl(ComplexHeatmap::Heatmap)(plt_args)

  res
}

