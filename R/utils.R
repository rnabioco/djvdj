#' Add V(D)J data to Seurat object
#'
#' @param sobj_in Seurat object
#' @param vdj_dir cellranger vdj output directories. If a vector of multiple
#' paths is provided, an equal number of cell prefixes must also be provided.
#' If a named vector is given, the names will be used to prefix each cell
#' barcode.
#' @param prefix Prefix to add to new meta.data columns
#' @param cell_prefix Prefix to add to cell barcodes
#' @param filter_contigs Only include chains with at least one productive
#' contig
#' @param sep Separator to use for storing per cell clonotype information in
#' the meta.data
#' @return Seurat object with V(D)J data added to meta.data
#' @export
import_vdj <- function(sobj_in, vdj_dir, prefix = "", cell_prefix = "",
                       filter_contigs = TRUE, sep = ";") {

  vdj_dir <- purrr::map_chr(vdj_dir, ~ file.path(.x, "outs"))

  # VDJ columns
  count_cols <- c("reads", "umis")

  sep_cols <- c(
    "v_gene",     "d_gene",
    "j_gene",     "c_gene",
    "chains",     "cdr3",
    "cdr3_nt",    count_cols,
    "productive", "full_length"
  )

  vdj_cols <- c(
    "barcode", "clonotype_id",
    sep_cols
  )

  # Check path names
  if (is.null(names(vdj_dir))) {
    if (length(vdj_dir) != length(cell_prefix)) {
      stop("Must provide a cell prefix for each path passed to vdj_dir.")
    }

    names(vdj_dir) <- cell_prefix
  }

  if (any(is.na(names(vdj_dir)))) {
    stop("Cell prefixes must not include NAs.")
  }

  nms <- !grepl("_$", names(vdj_dir))
  names(vdj_dir)[nms] <- paste0(names(vdj_dir)[nms], "_")

  # Load contigs
  contigs <- "filtered_contig_annotations.csv"
  contigs <- purrr::map(vdj_dir, ~ file.path(.x, contigs))
  contigs <- purrr::map(contigs, readr::read_csv, col_types = readr::cols())

  # Add cell prefixes
  contigs <- purrr::imap(contigs, ~ {
    .x <- dplyr::mutate(
      .x,
      barcode      = paste0(.y, .data$barcode),
      clonotype_id = paste0(.y, raw_clonotype_id)
    )

    .x <- dplyr::rename(.x, chains = .data$chain)
  })

  contigs <- dplyr::bind_rows(contigs)

  # Filter for productive contigs
  if (filter_contigs) {
    contigs <- dplyr::filter(contigs, .data$productive, .data$full_length)
  }

  contigs <- dplyr::select(contigs, all_of(vdj_cols))

  # Sum contig reads and UMIs for chains
  grp_cols <- vdj_cols[!vdj_cols %in% count_cols]
  contigs  <- dplyr::group_by(contigs, !!!syms(grp_cols))

  contigs  <- dplyr::summarize(
    contigs,
    across(all_of(count_cols), sum),
    .groups = "drop"
  )

  # Merge rows for each cell
  contigs <- dplyr::arrange(
    contigs,
    .data$barcode, .data$clonotype_id, .data$chains
  )

  contigs <- dplyr::group_by(contigs, .data$barcode, .data$clonotype_id)

  meta_df <- summarize(
    contigs,
    n_chains = n(),
    dplyr::across(
      all_of(sep_cols),
      ~ paste0(as.character(.x), collapse = sep)
    ),
    .groups = "drop"
  )

  # Filter for cells present in sobj_in
  cells <- Seurat::Cells(sobj_in)

  if (!any(cells %in% meta_df$barcode)) {
    stop("No VDJ cell barcodes are present in the Seurat object. Are you sure you are using the correct cell prefixes?")
  }

  meta_df <- dplyr::filter(meta_df, .data$barcode %in% cells)

  # Add meta.data to Seurat object
  meta_df <- tibble::column_to_rownames(meta_df, "barcode")
  meta_df <- dplyr::rename_with(meta_df, ~ paste0(prefix, .x))

  res <- Seurat::AddMetaData(sobj_in, metadata = meta_df)

  res
}


#' Filter V(D)J meta.data
#'
#' @param sobj_in Seurat object containing V(D)J data
#' @param filt Condition to use for filtering meta.data
#' @param clonotype_col meta.data column containing clonotype IDs. This column
#' is used to determine which cells have V(D)J data. If clonotype_col is set to
#' NULL, filtering is performed regardless of whether V(D)J data is present for
#' the cell
#' @param filter_cells Should cells be removed from object? If set FALSE
#' (default) V(D)J data will be removed from meta.data but no cells will be
#' removed from the object
#' @param sep Separator to use for expanding meta.data columns
#' @param vdj_cols meta.data columns containing VDJ data to use for filtering.
#' If set to NULL (recommended) columns are automatically selected.
#' @return Seurat object
#' @export
filter_vdj <- function(sobj_in, filt, clonotype_col = "clonotype_id", filter_cells = FALSE,
                       sep = ";", vdj_cols = NULL) {

  if (!filter_cells && is.null(clonotype_col)) {
    stop("clonotype_col must be provided when filter_cells is set to FALSE.")
  }

  # Identify columns with VDJ data
  meta_df <- sobj_in@meta.data
  meta_df <- tibble::as_tibble(meta_df, rownames = ".cell_id")

  col_list <- .get_vdj_cols(
    df_in     = meta_df,
    clone_col = clonotype_col,
    cols_in   = vdj_cols,
    sep       = sep
  )

  vdj_cols <- col_list$vdj
  sep_cols <- col_list$sep

  # Create list-cols for VDJ columns that contain sep
  if (!purrr::is_empty(sep_cols)) {
    sep_cols <- set_names(
      x  = sep_cols,
      nm = paste0(".", sep_cols)
    )

    meta_df <- .split_vdj(
      df_in    = meta_df,
      sep_cols = sep_cols,
      sep      = sep
    )

    meta_df <- dplyr::rowwise(meta_df)
  }

  # Store results for filtering
  meta_df <- dplyr::mutate(meta_df, .KEEP = {{filt}})
  meta_df <- dplyr::ungroup(meta_df)

  # Remove VDJ data from meta.data (without filtering cells)
  if (!filter_cells) {
    meta_df <- dplyr::mutate(
      meta_df,
      dplyr::across(
        all_of(c(vdj_cols, names(sep_cols))),
        ~ ifelse(.KEEP, .x, NA)
      )
    )

  # Filter cells from object
  # In clonotype_col != NULL only filter cells with VDJ data
  } else {
    if (!is.null(clonotype_col)) {
      meta_df <- dplyr::mutate(
        meta_df,
        .KEEP = dplyr::if_else(
          is.na(!!sym(clonotype_col)),
          TRUE,
          .data$.KEEP
        )
      )
    }

    meta_df <- dplyr::filter(meta_df, .data$.KEEP)
  }

  # Remove columns created for filtering
  if (!purrr::is_empty(sep_cols)) {
    sep_cols <- purrr::set_names(
      x  = names(sep_cols),
      nm = unname(sep_cols)
    )

    meta_df <- dplyr::select(
      meta_df,
      !dplyr::all_of(c(names(sep_cols), ".KEEP"))
    )

    meta_df <- dplyr::rename(meta_df, !!!syms(sep_cols))
  }

  # Add meta.data to Seurat object
  meta_df <- tibble::column_to_rownames(meta_df, ".cell_id")

  cells <- rownames(meta_df)
  res   <- subset(sobj_in, cells = cells)
  res   <- Seurat::AddMetaData(res, meta_df)

  res
}


#' Mutate V(D)J meta.data
#'
#' @param sobj_in Seurat object containing V(D)J data
#' @param ... Name-value pairs to use for creating or modifying meta.data
#' columns
#' @param clonotype_col meta.data column containing clonotype IDs. This column
#' is used to determine which cells have V(D)J data
#' @param sep Separator to use for expanding meta.data columns
#' @param vdj_cols meta.data columns containing VDJ data to use for filtering.
#' If set to NULL (recommended) columns are automatically selected.
#' @return Seurat object
#' @export
mutate_vdj <- function(sobj_in, ..., clonotype_col = "clonotype_id", sep = ";", vdj_cols = NULL) {

  # Identify columns with VDJ data
  meta_df <- sobj_in@meta.data
  meta_df <- tibble::as_tibble(meta_df, rownames = ".cell_id")

  col_list <- .get_vdj_cols(
    df_in     = meta_df,
    clone_col = clonotype_col,
    cols_in   = vdj_cols,
    sep       = sep
  )

  vdj_cols <- col_list$vdj
  sep_cols <- col_list$sep

  # Create list-cols for VDJ columns that contain sep
  if (!purrr::is_empty(sep_cols)) {
    sep_cols <- set_names(
      x  = sep_cols,
      nm = paste0(".", sep_cols)
    )

    meta_df <- .split_vdj(
      df_in    = meta_df,
      sep_cols = sep_cols,
      sep      = sep
    )

    meta_df <- dplyr::rowwise(meta_df)
  }

  # Mutate meta.data
  meta_df <- dplyr::mutate(meta_df, ...)
  meta_df <- dplyr::ungroup(meta_df)

  # Remove columns created for mutate
  if (!purrr::is_empty(sep_cols)) {
    sep_cols <- purrr::set_names(
      x  = names(sep_cols),
      nm = unname(sep_cols)
    )

    meta_df <- dplyr::select(
      meta_df,
      !dplyr::all_of(names(sep_cols))
    )

    meta_df <- dplyr::rename(meta_df, !!!syms(sep_cols))
  }

  # Add meta.data to Seurat object
  meta_df <- tibble::column_to_rownames(meta_df, ".cell_id")

  cells <- rownames(meta_df)
  res   <- subset(sobj_in, cells = cells)
  res   <- Seurat::AddMetaData(res, meta_df)

  res
}


#' Identify clonotypes
#'
#' Identify clonotypes using enclone from 10x Genomics
#' (https://10xgenomics.github.io/enclone)
#'
#' @param sobj_in Seurat object
#' @param vdj_dir cellranger vdj output directories. If a vector of multiple
#' paths is provided, an equal number of cell prefixes must also be provided.
#' If a named vector is given, the names will be used to prefix each cell
#' barcode.
#' @param csv_file Output file to write enclone results
#' @param prefix Prefix to add to new meta.data columns
#' @param cell_prefix Prefix to add to cell barcodes
#' @param overwrite If csv_file already exists, should the file be overwritten?
#' If csv_file already exists and overwrite = FALSE, the existing file will be
#' loaded. To re-run enclone, set overwrite = TRUE.
#' @return Seurat object
#' @export
identify_clonotypes <- function(sobj_in, vdj_dir, csv_file = "enclone_output.csv", prefix = "",
                                cell_prefix = "", overwrite = FALSE) {

  # Check path names
  if (is.null(names(vdj_dir))) {
    if (length(vdj_dir) != length(cell_prefix)) {
      stop("Must provide a cell prefix for each path passed to vdj_dir.")
    }

    names(vdj_dir) <- cell_prefix
  }

  if (any(is.na(names(vdj_dir)))) {
    stop("Cell prefixes must not include NAs.")
  }

  nms <- !grepl("_$", names(vdj_dir))
  names(vdj_dir)[nms] <- paste0(names(vdj_dir)[nms], "_")

  # Run enclone
  en_cols <- c(
    "datasets", "barcodes",
    "group_id", "group_ncells"
    # "clonotype_id", "clonotype_ncells"
  )

  if (!file.exists(csv_file) || overwrite) {
    out_cols <- paste0(en_cols, collapse = ",")
    en_dir   <- paste0(vdj_dir, collapse = ";")
    en_cmd   <- paste0("enclone BCR=\"", en_dir, "\" POUT=", csv_file, " PCOLS=", out_cols)

    system(en_cmd, ignore.stdout = TRUE)

  } else {
    warning(paste0(csv_file, " already exists, using previous results. Set overwrite = TRUE to re-run enclone."))
  }

  samples <- purrr::set_names(
    names(vdj_dir),
    basename(vdj_dir)
  )

  # Import enclone results
  meta_df <- readr::read_csv(csv_file, col_types = readr::cols())

  meta_df <- dplyr::mutate(
    meta_df,
    barcodes = strsplit(as.character(.data$barcodes), ",")
  )

  meta_df <- tidyr::unnest(meta_df, cols = "barcodes")

  meta_df <- dplyr::mutate(
    meta_df,
    barcodes = paste0(samples[.data$datasets], .data$barcodes)
  )

  meta_df <- dplyr::select(meta_df, -.data$datasets)

  # Filter for cells present in sobj_in
  cells <- Seurat::Cells(sobj_in)

  if (!any(cells %in% meta_df$barcodes)) {
    stop("No VDJ cell barcodes are present in the Seurat object. Are you sure you are using the correct cell prefixes?")
  }

  meta_df <- dplyr::filter(meta_df, .data$barcodes %in% cells)

  # Add meta.data to Seurat object
  meta_df <- tibble::column_to_rownames(meta_df, "barcodes")
  meta_df <- dplyr::rename_with(meta_df, ~ paste0(prefix, .x))

  res <- Seurat::AddMetaData(sobj_in, metadata = meta_df)

  res
}


#' Calculate clonotype abundance
#'
#' @param sobj_in Seurat object containing V(D)J data
#' @param clonotype_col meta.data column containing clonotype IDs
#' @param cluster_col meta.data column containing cluster IDs to use for
#' grouping cells when calculating clonotype abundance
#' @param prefix Prefix to add to new meta.data columns
#' @param return_seurat Return a Seurat object. If set to FALSE, a tibble
#' summarizing the results is returned.
#' @return Seurat object with clonotype abundance added to meta.data
#' @export
calc_abundance <- function(sobj_in, clonotype_col = NULL, cluster_col = NULL,
                           prefix = "", return_seurat = TRUE) {

  # Format meta.data
  meta_df <- sobj_in@meta.data
  meta_df <- tibble::as_tibble(meta_df, rownames = ".cell_id")
  meta_df <- dplyr::filter(meta_df, !is.na(!!sym(clonotype_col)))

  meta_df <- dplyr::select(
    meta_df,
    .data$.cell_id, all_of(c(cluster_col, clonotype_col))
  )

  # Calculate abundance
  if (!is.null(cluster_col)) {
    meta_df <- dplyr::group_by(meta_df, !!sym(cluster_col))
  }

  meta_df <- dplyr::mutate(
    meta_df,
    n_cells = dplyr::n_distinct(.data$.cell_id)
  )

  meta_df <- dplyr::group_by(
    meta_df,
    !!sym(clonotype_col),
    .add = TRUE
  )

  meta_df <- dplyr::mutate(
    meta_df,
    freq = dplyr::n_distinct(.data$.cell_id),
    pct  = (freq / .data$n_cells) * 100
  )

  meta_df <- dplyr::ungroup(meta_df)

  # Add results to meta.data
  if (!return_seurat) {
    res <- dplyr::select(meta_df, -.cell_id)
    res <- unique(res)

    return(res)
  }

  freq_col <- paste0(prefix, "clone_freq")
  pct_col  <- paste0(prefix, "clone_pct")

  meta_df <- dplyr::select(
    meta_df,
    -all_of(c("n_cells", cluster_col)),
    !!sym(freq_col) := freq,
    !!sym(pct_col)  := pct
  )

  meta_df <- tibble::column_to_rownames(meta_df, ".cell_id")

  res <- Seurat::AddMetaData(sobj_in, metadata = meta_df)

  res
}


#' Calculate repertoire diversity
#'
#' @param sobj_in Seurat object
#' @param clonotype_col meta.data column containing clonotype IDs to use for
#' calculating diversity
#' @param cluster_col meta.data column containing cluster IDs to use for
#' calculating diversity. If cluster_col is omitted, diversity index will be
#' calculated for all clonotypes.
#' @param method Method to use for calculating diversity. Can pass also pass
#' a named list to use multiple methods.
#' @param prefix Prefix to add to new meta.data columns
#' @param return_seurat Return a Seurat object. If set to FALSE, the meta.data
#' table is returned.
#' @return Seurat object with diversity index added to meta.data
#' @export
calc_diversity <- function(sobj_in, clonotype_col = NULL, cluster_col = NULL,
                           method = abdiv::simpson, prefix = "", return_seurat = TRUE) {

  # Format meta.data
  vdj_cols <- clonotype_col
  meta_df  <- tibble::as_tibble(sobj_in@meta.data, rownames = ".cell_id")
  vdj_df   <- dplyr::filter(meta_df, !is.na(!!sym(clonotype_col)))

  # Count clonotypes
  if (!is.null(cluster_col)) {
    vdj_cols <- c(cluster_col, vdj_cols)
    vdj_df   <- dplyr::group_by(vdj_df, !!sym(cluster_col))
  }

  vdj_df <- dplyr::group_by(vdj_df, !!sym(clonotype_col), .add = TRUE)

  vdj_df <- dplyr::summarize(
    vdj_df,
    .n      = dplyr::n_distinct(.data$.cell_id),
    .groups = "drop"
  )

  # Calculate diversity
  if (!is.null(cluster_col)) {
    vdj_df <- dplyr::group_by(vdj_df, !!sym(cluster_col))
  }

  if (length(method) == 1 && is.null(names(method))) {
    nm <- as.character(substitute(method))
    nm <- dplyr::last(nm)

    method        <- list(method)
    names(method) <- nm
  }

  div_cols      <- paste0(prefix, names(method))
  names(method) <- div_cols

  vdj_df <- purrr::imap_dfr(method, ~ {
    dplyr::mutate(
      vdj_df,
      diversity = .x(.data$.n),
      met       = .y
    )
  })

  vdj_df <- tidyr::pivot_wider(
    vdj_df,
    names_from  = "met",
    values_from = "diversity"
  )

  vdj_df <- dplyr::ungroup(vdj_df)
  vdj_df <- dplyr::select(vdj_df, all_of(c(vdj_cols, div_cols)))

  # Return data.frame
  if (!return_seurat) {
    res <- dplyr::select(vdj_df, dplyr::all_of(c(cluster_col, div_cols)))
    res <- dplyr::distinct(res)

    if (is.null(cluster_col)) {
      res <- unlist(res)
    }

    return(res)
  }

  # Add results to meta.data
  meta_df <- dplyr::left_join(meta_df, vdj_df, by = vdj_cols)
  meta_df <- tibble::column_to_rownames(meta_df, ".cell_id")
  meta_df <- as.data.frame(meta_df)

  res <- Seurat::AddMetaData(sobj_in, metadata = meta_df)

  res
}


#' Calculate repertoire similarity between clusters
#'
#' @param sobj_in Seurat object
#' @param clonotype_col meta.data column containing clonotype IDs to use for
#' calculating overlap
#' @param cluster_col meta.data column containing cluster IDs to use for
#' calculating overlap
#' @param method Method to use for calculating similarity between clusters
#' @param prefix Prefix to add to new meta.data columns
#' @param return_seurat Return a Seurat object. If set to FALSE, a matrix is
#' returned
#' @return Seurat object with similarity index added to meta.data
#' @export
calc_similarity <- function(sobj_in, clonotype_col = NULL, cluster_col,
                            method = abdiv::jaccard, prefix = "sim_", return_seurat = TRUE) {

  # Format meta.data
  meta_df <- sobj_in@meta.data
  meta_df <- tibble::as_tibble(meta_df, rownames = ".cell_id")
  meta_df <- filter(meta_df, !is.na(!!sym(clonotype_col)))

  meta_df <- dplyr::select(
    meta_df,
    dplyr::all_of(c(".cell_id", clonotype_col, cluster_col))
  )

  vdj_df <- dplyr::group_by(
    meta_df,
    !!!syms(c(cluster_col, clonotype_col))
  )

  vdj_df <- dplyr::summarize(
    vdj_df,
    n       = n_distinct(.data$.cell_id),
    .groups = "drop"
  )

  vdj_df <- tidyr::pivot_wider(
    vdj_df,
    names_from  = dplyr::all_of(cluster_col),
    values_from = .data$n
  )

  # Calculate similarity index
  clusts <- colnames(vdj_df)
  clusts <- clusts[clusts != clonotype_col]

  vdj_df <- dplyr::mutate(vdj_df, dplyr::across(
    dplyr::all_of(clusts), ~ replace_na(.x, 0)
  ))

  combs <- utils::combn(clusts, 2, simplify = FALSE)

  res <- map_dfr(combs, ~ {
    ins <- paste0("vdj_df$", .x)

    tibble::tibble(
      Var1 = .x[1],
      Var2 = .x[2],
      sim  = method(dplyr::pull(vdj_df, .x[1]), dplyr::pull(vdj_df, .x[2]))
    )
  })

  res_i <- dplyr::rename(res, Var1 = .data$Var2, Var2 = .data$Var1)

  # Return matrix
  if (!return_seurat) {
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
  clusts <- tibble::tibble(Var1 = clusts, Var2 = clusts, sim = 1)
  res    <- dplyr::bind_rows(res, res_i, clusts)
  res    <- dplyr::mutate(res, Var1 = paste0(prefix, .data$Var1))

  res <- tidyr::pivot_wider(
    res,
    names_from  = .data$Var1,
    values_from = .data$sim
  )

  # Add similarity index to meta.data
  j_cols <- purrr::set_names("Var2", cluster_col)

  meta_df <- dplyr::left_join(meta_df, res, by = j_cols)
  meta_df <- dplyr::select(meta_df, -all_of(cluster_col))
  meta_df <- tibble::column_to_rownames(meta_df, ".cell_id")
  meta_df <- as.data.frame(meta_df)

  res <- Seurat::AddMetaData(sobj_in, meta_df)

  res
}


#' Calculate gene usage
#'
#' Gene usage is calculated as the percent of cells that express each V(D)J
#' gene present in gene_cols. Cells that lack V(D)J data and have NAs present
#' in gene_cols are excluded from this calculation.
#'
#' @param sobj_in Seurat object containing V(D)J data
#' @param gene_cols meta.data column containing genes used for each clonotype
#' @param cluster_col meta.data column containing cell clusters to use when
#' calculating gene usage
#' @param chain Chain to use for calculating gene usage. Set to NULL to include
#' all chains.
#' @param chain_col meta.data column containing chains for each cell
#' @param sep Separator to use for expanding gene_cols
#' @return data.frame containing gene usage summary
#' @export
calc_usage <- function(sobj_in, gene_cols, cluster_col = NULL, chain = NULL,
                       chain_col = NULL, sep = ";") {

  # data.frame to calculate usage
  sep_cols <- c(gene_cols, chain_col)
  vdj_cols   <- c(".cell_id", cluster_col, sep_cols)

  meta_df <- sobj_in@meta.data
  meta_df <- tibble::as_tibble(meta_df, rownames = ".cell_id")
  meta_df <- dplyr::select(meta_df, dplyr::all_of(vdj_cols))

  meta_df <- dplyr::filter(meta_df, dplyr::across(
    dplyr::all_of(gene_cols),
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

    res <- dplyr::mutate(res, dplyr::across(dplyr::all_of(gene_cols), ~ {
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
    clusts       <- dplyr::pull(meta_df, cluster_col)
    clust_counts <- table(clusts)
    clusts       <- unique(clusts)

    res <- tidyr::pivot_wider(
      res,
      names_from  = dplyr::all_of(cluster_col),
      values_from = .data$freq,
      values_fill = 0
    )

    res <- tidyr::pivot_longer(
      res,
      cols      = dplyr::all_of(clusts),
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


#' Summarize values for chains
#'
#' Summarize values present for each column provided to the data_cols argument.
#' For each cell, the function(s) provided will be applied to each unique label
#' in chain_col.
#'
#' @param sobj_in Seurat object containing V(D)J data
#' @param data_cols meta.data columns to summarize
#' @param fn Function to use for summarizing data_cols
#' @param chain_col meta.data column(s) containing labels for each chain
#' expressed in the cell. These labels are used for grouping the summary
#' output. Set chain_col to NULL to group solely based on the cell ID.
#' @param include_cols Additional columns to include in the output data.frame
#' @param sep Separator to use for expanding data_cols and chain_col
#' @return data.frame containing summary results
#' @export
summarize_chains <- function(sobj_in, data_cols = c("umis", "reads"), fn,
                             chain_col = NULL, include_cols = NULL, sep = ";") {

  # Fetch meta.data
  fetch_cols <- c(data_cols, chain_col, include_cols)

  meta_df <- Seurat::FetchData(sobj_in, fetch_cols)
  meta_df <- tibble::as_tibble(meta_df, rownames = ".cell_id")

  meta_df <- dplyr::filter(
    meta_df,
    dplyr::across(dplyr::all_of(data_cols), ~ !is.na(.x))
  )

  # Expand meta.data
  res <- dplyr::mutate(meta_df, dplyr::across(
    dplyr::all_of(c(data_cols, chain_col)),
    ~ strsplit(as.character(.x), sep)
  ))

  res <- tidyr::unnest(res, cols = dplyr::all_of(c(data_cols, chain_col)))

  # Summarize data_cols for each chain present for the cell
  res <- dplyr::mutate(res, dplyr::across(
    dplyr::all_of(data_cols),
    ~ .convert_char(.x, as.numeric)
  ))

  grp_cols <- c(".cell_id", chain_col, include_cols)
  res      <- dplyr::group_by(res, !!!syms(grp_cols))

  res <- dplyr::summarize(
    res,
    dplyr::across(dplyr::all_of(data_cols), fn),
    .groups = "drop"
  )

  res
}


#' Split V(D)J meta.data columns into vectors
#'
#' @param df_in data.frame to modify
#' @param sep Separator to use for expanding meta.data columns
#' @param sep_cols Columns to split based on sep
#' @param num_cols Columns to convert to numeric
#' @param lgl_cols Columns to convert to logical
#' @return data.frame
.split_vdj <- function(df_in, sep = ";", sep_cols, num_cols = c("umis", "reads"),
                       lgl_cols = c("productive", "full_length")) {

  # Split columns into vectors
  res <- dplyr::mutate(df_in, !!!syms(sep_cols))

  res <- dplyr::mutate(res, dplyr::across(
    dplyr::all_of(unname(sep_cols)),
    ~ strsplit(as.character(.x), sep)
  ))

  # Convert to numeric or logical
  res <- dplyr::mutate(
    res,
    across(
      dplyr::all_of(num_cols),
      map, ~ .convert_char(.x, as.numeric)
    ),
    across(
      dplyr::all_of(lgl_cols),
      map, ~ .convert_char(.x, as.logical)
    )
  )

  res
}


#' Identify columns with V(D)J data
#'
#' @param df_in data.frame
#' @param clone_col Column containing clonotype IDs to use for identifying
#' columns with V(D)J data
#' @param cols_in meta.data columns containing VDJ data to use for filtering.
#' If set to NULL (recommended) columns are automatically selected.
#' @param sep Separator to search for in columns
#' @return data.frame
.get_vdj_cols <- function(df_in, clone_col, cols_in, sep) {

  # Identify columns with VDJ data based on NAs in clonotype_col
  # If no clonotype_col is given use all columns
  if (is.null(clone_col)) {
    cols_in <- colnames(df_in)
    cols_in <- cols_in[cols_in != ".cell_id"]
  }

  if (is.null(cols_in)) {
    cols_in <- dplyr::mutate(
      df_in,
      dplyr::across(dplyr::everything(), is.na)
    )

    cols_in <- purrr::keep(
      cols_in,
      ~ identical(.x, dplyr::pull(cols_in, clone_col))
    )

    cols_in <- colnames(cols_in)
  }

  # Identify columns to split based on sep
  sep_cols <- NULL

  if (!is.null(sep)) {
    sep_cols <- dplyr::select(df_in, all_of(cols_in))

    sep_cols <- purrr::keep(
      sep_cols,
      ~ any(purrr::map_lgl(na.omit(.x), grepl, pattern = sep))
    )

    sep_cols <- colnames(sep_cols)
  }

  # Return list of vectors
  res <- list(
    "vdj" = cols_in,
    "sep" = sep_cols
  )

  res
}


#' Attempt to convert character vector using provided function
#'
#' @param x Character vector to convert
#' @param fn Function to try
#' @return Value converted using fn
.convert_char <- function(x, fn) {
  if (!is.character(x)) {
    return(x)
  }

  suppressWarnings(ifelse(!is.na(fn(x)) | is.na(x), fn(x), x))
}







