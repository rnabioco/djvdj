#' Add V(D)J data to Seurat object
#'
#' @param sobj_in Seurat object, if set to NULL a tibble containing the V(D)J
#' data will be returned
#' @param vdj_dir Directory containing the output from cellranger vdj. A vector
#' or named vector can be given to load data from several runs. If a named
#' vector is given, the cell barcodes will be prefixed with the provided names.
#' This mimics the behavior of the Read10X function found in the Seurat
#' package. Cell barcode prefixes can also be provided using the cell_prefix
#' argument.
#' @param prefix Prefix to add to new meta.data columns
#' @param cell_prefix Prefix to add to cell barcodes, this is helpful when
#' loading data from multiple runs into a single object. If set to NULL, cell
#' barcode prefixes will be automatically generated in a similar way as the
#' Read10X function found in the Seurat package.
#' @param filter_contigs Only include chains with at least one productive
#' contig
#' @param sep Separator to use for storing per cell clonotype information in
#' the meta.data
#' @return Seurat object with V(D)J data added to meta.data
#' @export
import_vdj <- function(sobj_in = NULL, vdj_dir, prefix = "", cell_prefix = NULL,
                       filter_contigs = TRUE, sep = ";") {

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
  # if no cell prefixes are provided, auto generate prefixes
  # add an underscore if not included in the provided cell prefixes
  if (is.null(names(vdj_dir))) {

    if (is.null(cell_prefix)) {
      cell_prefix    <- paste0(seq_along(vdj_dir), "_")
      cell_prefix[1] <- ""
    }

    if (length(vdj_dir) != length(cell_prefix)) {
      stop("cell_prefix must be the same length as vdj_dir (", length(vdj_dir),").")
    }

    names(vdj_dir) <- cell_prefix
  }

  if (any(is.na(names(vdj_dir)))) {
    stop("Cell prefixes cannot include NAs.")
  }

  if (any(duplicated(names(vdj_dir)))) {
    dups <- duplicated(names(vdj_dir))
    dups <- names(vdj_dir)[dups]
    dups <- paste0(dups, collapse = ", ")

    warning("Some cell barcode prefixes are duplicated: ", dups)
  }

  nms                 <- names(vdj_dir) != "" & !grepl("_$", names(vdj_dir))
  names(vdj_dir)[nms] <- paste0(names(vdj_dir)[nms], "_")

  # Load contigs
  # check given dir before adding "outs" to path
  contigs <- "filtered_contig_annotations.csv"

  contigs <- purrr::map_chr(vdj_dir, ~ {
    path <- case_when(
      file.exists(file.path(.x, contigs))         ~ file.path(.x, contigs),
      file.exists(file.path(.x, "outs", contigs)) ~ file.path(.x, "outs", contigs)
    )

    if (is.na(path)) {
      stop(contigs, " not found in ", vdj_dir, ".")
    }

    path
  })

  contigs <- purrr::map(
    contigs,
    readr::read_csv,
    col_types = readr::cols()
  )

  # Add cell prefixes and bind rows
  contigs <- purrr::imap_dfr(contigs, ~ {
    .x <- dplyr::mutate(
      .x,
      barcode      = paste0(.y, .data$barcode),
      clonotype_id = paste0(.y, raw_clonotype_id)
    )

    dplyr::rename(.x, chains = .data$chain)
  })



  # CHECK BARCODE OVERLAP BEFORE BINDING ROWS



  # Filter for productive contigs
  if (filter_contigs) {
    contigs <- dplyr::filter(contigs, .data$productive, .data$full_length)
  }

  # Remove contigs that do not have an assigned clonotype_id
  n_remove <- contigs$raw_clonotype_id
  n_remove <- n_remove[is.na(n_remove)]
  n_remove <- length(n_remove)

  if (n_remove > 0) {
    warning(n_remove, " contigs do not have an assigned clonotype_id, these contigs will be removed.")

    contigs <- dplyr::filter(contigs, !is.na(.data$raw_clonotype_id))
  }

  # Select V(D)J columns to keep
  contigs <- dplyr::select(contigs, all_of(vdj_cols))

  # Check if sep is already in sep_cols
  sep <- gsub("[[:space:]]", "", sep)

  if (any(grepl(sep, contigs[, sep_cols]))) {
    stop("The string '", sep, "' is already present in the V(D)J data, select a different value for sep.")
  }

  # Sum contig reads and UMIs for chains
  # some chains are supported by multiple contigs
  grp_cols <- vdj_cols[!vdj_cols %in% count_cols]
  contigs  <- dplyr::group_by(contigs, !!!syms(grp_cols))

  contigs <- dplyr::summarize(
    contigs,
    across(all_of(count_cols), sum),
    .groups = "drop"
  )

  # Order chains and CDR3 sequences
  # when the rows are collapsed, the cdr3 sequences must be in the same order
  # for every cell. This is required so the cdr3 columns can be used directly
  # as the clonotype ID
  contigs <- dplyr::arrange(
    contigs,
    .data$barcode, .data$chains, .data$cdr3_nt
  )

  # Extract isotypes from c_gene for IGH chain
  iso_pat <- "^IGH[ADEGM]"

  if (any(grepl(iso_pat, contigs$c_gene))) {
    contigs <- dplyr::group_by(contigs, .data$barcode)

    contigs <- dplyr::mutate(
      contigs,
      isotype = list(
        substr(
          .data$c_gene, 1,
          attr(regexpr(iso_pat, .data$c_gene), "match.length", exact = TRUE)
        )
      ),
      isotype = map_chr(.data$isotype, ~ {
        isos <- unique(.x)
        isos <- tidyr::replace_na(isos, "")
        isos <- isos[isos != ""]

        if (length(isos) == 0) {
          isos <- ""
        }

        isos <- dplyr::case_when(
          length(isos) > 1 ~ "Multi",
          isos == ""       ~ "None",
          TRUE             ~ isos
        )

        unique(isos)
      })
    )

    contigs <- dplyr::group_by(contigs, .data$isotype)
  }

  # Collapse chains into a single row for each cell
  # Include isotype and clonotype_id as groups so that they are included in the
  # summarized results
  contigs <- dplyr::group_by(
    contigs,
    .data$barcode, .data$clonotype_id,
    .add = TRUE
  )

  meta_df <- summarize(
    contigs,
    across(
      all_of(sep_cols),
      ~ paste0(as.character(.x), collapse = sep)
    ),
    n_chains = n(),
    .groups = "drop"
  )

  # Check for duplicated cell barcodes
  if (any(duplicated(meta_df$barcode))) {
    stop("Malformed inport data, multiple clonotype_ids are associated with the same cell barcode.")
  }

  # Return tibble if sobj_in is NULL
  if (is.null(sobj_in)) {
    return(meta_df)
  }

  # Check overlap between V(D)J data and the Seurat object
  so_bcs  <- Seurat::Cells(sobj_in)
  vdj_bcs <- meta_df$barcode

  n_overlap <- length(so_bcs[so_bcs %in% vdj_bcs])
  pct_so    <- round(n_overlap / length(so_bcs), 2) * 100
  pct_vdj   <- round(n_overlap / length(vdj_bcs), 2) * 100

  if (n_overlap == 0) {
    stop("
      Cell barcodes from the V(D)J data were not found in the Seurat object, are you using the
      correct cell barcode prefixes? Cell barcode prefixes can be provided to the cell_prefix
      argument or by passing a named vector to the vdj_dir argument.
    ")
  }

  if (pct_so < 25) {
    warning("Only ", pct_so, "% (", n_overlap, ") of cell barcodes present in the Seurat object overlap with the V(D)J data.")
  }

  if (pct_vdj < 25) {
    warning("Only ", pct_vdj, "% (", n_overlap, ") of cell barcodes present in the V(D)J data overlap with the Seurat object.")
  }

  # Filter for cells present in sobj_in
  meta_df <- dplyr::filter(meta_df, .data$barcode %in% so_bcs)

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
filter_vdj <- function(sobj_in, filt, clonotype_col = "cdr3_nt", filter_cells = FALSE,
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
      across(
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
      !all_of(c(names(sep_cols), ".KEEP"))
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


#' Manipulate Seurat object meta.data
#'
#' @param sobj_in Seurat object
#' @param .fun Function or formula to use for modifying the meta.data. If a
#' formula is provided, use .x to refer to the meta.data table.
#' @param ... Arguments to pass to the provided function
#' @return Seurat object
#' @export
mutate_meta <- function(sobj_in, .fun, ...) {

  if (!is_function(.fun) && !is_formula(.fun)) {
    stop(".fun must be either a function or a formula.")
  }

  .x <- sobj_in@meta.data
  .x <- tibble::rownames_to_column(.x, ".cell_id")

  if (is_formula(.fun)) {
    .fun <- as_mapper(.fun, ...)
  }

  res <- .fun(.x, ...)
  res <- tibble::column_to_rownames(res, ".cell_id")

  sobj_in@meta.data <- res

  sobj_in
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
mutate_vdj <- function(sobj_in, ..., clonotype_col = "cdr3_nt", sep = ";", vdj_cols = NULL) {

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
      !all_of(names(sep_cols))
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
calc_abundance <- function(sobj_in, clonotype_col = "cdr3_nt", cluster_col = NULL,
                           prefix = "", return_seurat = TRUE) {

  # Format meta.data
  meta_df <- sobj_in@meta.data
  meta_df <- tibble::as_tibble(meta_df, rownames = ".cell_id")
  meta_df <- dplyr::filter(meta_df, !is.na(!!sym(clonotype_col)))

  meta_df <- dplyr::select(
    meta_df,
    .data$.cell_id, all_of(c(cluster_col, clonotype_col))
  )

  # Calculate clonotype abundance
  meta_df <- .calc_abund(
    df_in     = meta_df,
    cell_col  = ".cell_id",
    clone_col = clonotype_col,
    clust_col = cluster_col
  )

  # Add results to meta.data
  if (!return_seurat) {
    res <- dplyr::select(meta_df, -.data$.cell_id)
    res <- dplyr::distinct(res)

    return(res)
  }

  new_cols <- c("freq", "pct")

  if (!is.null(cluster_col)) {
    new_cols <- c(new_cols, "shared")
  }

  new_cols <- purrr::set_names(
    new_cols,
    paste0(prefix, "clone_", new_cols)
  )

  meta_df <- select(
    meta_df,
    -all_of(c("n_cells", cluster_col)),
    !!!syms(new_cols)
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
calc_diversity <- function(sobj_in, clonotype_col = "cdr3_nt", cluster_col = NULL,
                           method = abdiv::simpson, prefix = "", return_seurat = TRUE) {

  if (length(method) > 1 && is.null(names(method))) {
    stop("Must include names if using a list of methods.")
  }

  if (length(method) == 1 && is.null(names(method))) {
    nm <- as.character(substitute(method))
    nm <- last(nm)

    method <- purrr::set_names(list(method), nm)
  }

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
    res <- dplyr::select(vdj_df, all_of(c(cluster_col, div_cols)))
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
calc_similarity <- function(sobj_in, clonotype_col = "cdr3_nt", cluster_col, method = abdiv::jaccard,
                            prefix = NULL, return_seurat = TRUE) {

  if (is.null(prefix)) {
    prefix <- as.character(substitute(method))
    prefix <- dplyr::last(prefix)
    prefix <- paste0(prefix, "_")
  }

  # Format meta.data
  meta_df <- sobj_in@meta.data
  meta_df <- tibble::as_tibble(meta_df, rownames = ".cell_id")
  meta_df <- filter(meta_df, !is.na(!!sym(clonotype_col)))

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
    n       = n_distinct(.data$.cell_id),
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
                       chain_col = "chains", sep = ";") {

  # data.frame to calculate usage
  sep_cols <- c(gene_cols, chain_col)
  vdj_cols   <- c(".cell_id", cluster_col, sep_cols)

  meta_df <- sobj_in@meta.data
  meta_df <- tibble::as_tibble(meta_df, rownames = ".cell_id")
  meta_df <- dplyr::select(meta_df, all_of(vdj_cols))

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
                             chain_col = "chains", include_cols = NULL, sep = ";") {

  # Fetch meta.data
  fetch_cols <- c(data_cols, chain_col, include_cols)

  meta_df <- Seurat::FetchData(sobj_in, fetch_cols)
  meta_df <- tibble::as_tibble(meta_df, rownames = ".cell_id")

  meta_df <- dplyr::filter(
    meta_df,
    across(all_of(data_cols), ~ !is.na(.x))
  )

  # Expand meta.data
  res <- dplyr::mutate(meta_df, across(
    all_of(c(data_cols, chain_col)),
    ~ strsplit(as.character(.x), sep)
  ))

  res <- tidyr::unnest(res, cols = all_of(c(data_cols, chain_col)))

  # Summarize data_cols for each chain present for the cell
  res <- dplyr::mutate(res, across(
    all_of(data_cols),
    ~ .convert_char(.x, as.numeric)
  ))

  grp_cols <- c(".cell_id", chain_col, include_cols)
  res      <- dplyr::group_by(res, !!!syms(grp_cols))

  res <- dplyr::summarize(
    res,
    across(all_of(data_cols), fn),
    .groups = "drop"
  )

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
    n_cells = dplyr::n_distinct(!!sym(cell_col))
  )

  # Calculate frequency
  res <- dplyr::group_by(df_in, !!sym(clone_col), .add = TRUE)

  res <- dplyr::mutate(
    res,
    freq = dplyr::n_distinct(!!sym(cell_col)),
    pct  = (.data$freq / .data$n_cells) * 100
  )

  # Identify shared clonotypes
  if (!is.null(clust_col)) {
    res <- dplyr::group_by(res, !!sym(clone_col))

    res <- dplyr::mutate(
      res,
      shared = dplyr::n_distinct(!!sym(clust_col)) > 1
    )
  }

  res <- dplyr::ungroup(res)

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

  res <- dplyr::mutate(res, across(
    all_of(unname(sep_cols)),
    ~ strsplit(as.character(.x), sep)
  ))

  # Convert to numeric or logical
  res <- dplyr::mutate(
    res,
    across(
      all_of(num_cols),
      map, ~ .convert_char(.x, as.numeric)
    ),
    across(
      all_of(lgl_cols),
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
      across(dplyr::everything(), is.na)
    )

    cols_in <- purrr::keep(
      cols_in,
      ~ identical(.x, pull(cols_in, clone_col))
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
