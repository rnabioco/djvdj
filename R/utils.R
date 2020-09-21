#' Add VDJ data to Seurat object
#'
#' @param sobj_in Seurat object
#' @param vdj_dir cellranger vdj output directory
#' @param prefix Prefix to add to new meta.data columns
#' @param cell_prefix Prefix to add to cell barcodes
#' @return Seurat object with VDJ data added to meta.data
#' @export
import_vdj <- function(sobj_in, vdj_dir, prefix = "", cell_prefix = "") {

  # Load contigs
  grp_cols <- c(
    "reads",  "umis",
    "v_gene", "d_gene",
    "j_gene", "c_gene",
    "chain",  "cdr3",
    "cdr3_nt"
  )

  vdj_cols <- c(
    "barcode", "raw_clonotype_id",
    grp_cols
  )

  vdj_dir <- file.path(vdj_dir, "filtered_contig_annotations.csv")
  contigs <- readr::read_csv(vdj_dir)

  # Extract chain from v_gene, add cell barcode prefix
  contigs <- dplyr::mutate(
    contigs,
    chain   = stringr::str_extract(.data$v_gene, "^[A-Z]{3}"),
    barcode = stringr::str_c(cell_prefix, .data$barcode)
  )

  # Filter for productive full length contigs
  contigs <- dplyr::filter(contigs, .data$productive, .data$full_length)
  contigs <- dplyr::select(contigs, all_of(vdj_cols))

  # Merge rows for each cell
  contigs <- dplyr::arrange(contigs, .data$barcode, .data$raw_clonotype_id, .data$chain)
  contigs <- dplyr::group_by(contigs, .data$barcode, .data$raw_clonotype_id)

  meta_df <- summarize(
    contigs,
    n_chains = n(),
    dplyr::across(
      all_of(grp_cols),
      ~ purrr::reduce(as.character(.x), stringr::str_c, sep = ";")
    ),
    .groups = "drop"
  )

  # Filter for cells present in sobj_in
  cells   <- Seurat::Cells(sobj_in)
  meta_df <- dplyr::filter(meta_df, barcode %in% cells)
  meta_df <- dplyr::rename(meta_df, clonotype_id = .data$raw_clonotype_id)

  # Calculate stats
  meta_df <- dplyr::group_by(meta_df, .data$clonotype_id)
  meta_df <- dplyr::mutate(meta_df, clone_freq = dplyr::n_distinct(.data$barcode))
  meta_df <- dplyr::ungroup(meta_df)
  meta_df <- dplyr::mutate(meta_df, clone_prop = .data$clone_freq / nrow(meta_df))

  # Add meta.data to Seurat object
  meta_df <- tibble::column_to_rownames(meta_df, "barcode")
  meta_df <- dplyr::rename_with(meta_df, ~ stringr::str_c(prefix, .x))

  res <- Seurat::AddMetaData(sobj_in, metadata = meta_df)

  res
}


#' Calculate receptor diversity (inverse Simpson Index)
#'
#' @param sobj_in Seurat object
#' @param clonotype_col meta.data column containing clonotype ids
#' @param cluster_col meta.data column containing cluster ids used for
#' calculating receptor diversity. If cluster_col is omitted, diversity index
#' will be calculated for all clonotypes.
#' @param prefix Prefix to add to new meta.data columns
#' @return Seurat object with inverse Simpson index added to meta.data
#' @export
calc_diversity <- function(sobj_in, clonotype_col = "clonotype_id",
                           cluster_col = NULL, prefix = "") {

  # meta.data
  vdj_cols      <- clonotype_col
  clonotype_col <- dplyr::sym(clonotype_col)
  meta_df       <- tibble::as_tibble(sobj_in@meta.data, rownames = "cell_id")
  vdj_df        <- dplyr::filter(meta_df, !is.na(!!clonotype_col))

  # Count clonotypes
  if (!is.null(cluster_col)) {
    vdj_cols    <- c(cluster_col, vdj_cols)
    cluster_col <- dplyr::sym(cluster_col)
    vdj_df      <- dplyr::group_by(vdj_df, !!cluster_col)
  }

  vdj_df <- dplyr::group_by(vdj_df, !!clonotype_col, .add = T)
  vdj_df <- dplyr::summarize(
    .data   = vdj_df,
    num     = dplyr::n_distinct(.data$cell_id),
    .groups = "drop"
  )

  # Calculate diversity
  if (!is.null(cluster_col)) {
    vdj_df <- dplyr::group_by(vdj_df, !!cluster_col)
  }

  div_col <- stringr::str_c(prefix, "diversity")

  vdj_df <- dplyr::mutate(
    vdj_df,
    frac     = .data$num / sum(.data$num),
    sum_frac = sum(.data$frac ^ 2),

    !!dplyr::sym(div_col) := 1 - .data$sum_frac
    # !!dplyr::sym(div_col) := 1 / sum_frac
  )

  vdj_df <- dplyr::select(vdj_df, all_of(c(vdj_cols, div_col)))

  # Add resuts to meta.data
  meta_df <- dplyr::left_join(meta_df, vdj_df, by = vdj_cols)
  meta_df <- tibble::column_to_rownames(meta_df, "cell_id")
  meta_df <- as.data.frame(meta_df)

  res <- Seurat::AddMetaData(sobj_in, metadata = meta_df)

  res
}


#' Calculate Jaccard index for cell identities
#'
#' @param sobj_in Seurat object
#' @param clonotype_col meta.data column containing clonotype ids
#' @param cluster_col meta.data column containing cell clusters
#' @param ref_cluster Cluster id to use as a reference for calculating Jaccard
#' index. If ref_cluster is omitted, Jaccard index will be calculated for all
#' combinations of clusters.
#' @param return_matrix Return matrix instead of Seurat object
#' @param prefix Prefix to add to new meta.data columns
#' @return Seurat object with Jaccard index added to meta.data
#' @export
calc_jaccard <- function(sobj_in, clonotype_col = "clonotype_id", cluster_col,
                         ref_cluster = NULL, prefix = "") {

  # Helper to calculate jaccard index
  calc_jidx <- function(df_in, comparison, clonotype_col) {

    if (length(comparison) != 2) {
      stop("comparison must be a character vector containing two elements")
    }

    uniq_vars     <- unique(comparison)
    clonotype_col <- dplyr::sym(clonotype_col)

    if (length(uniq_vars) == 1) {
      dup_var       <- stringr::str_c(uniq_vars, "_dup")
      comparison[3] <- dup_var
      dup_var       <- dplyr::sym(dup_var)
      uniq_vars     <- dplyr::sym(uniq_vars)

      df_in <- dplyr::mutate(df_in, !!dup_var := !!uniq_vars)
    }

    row_sums <- dplyr::select(df_in, !!clonotype_col, dplyr::all_of(comparison))
    row_sums <- tidyr::pivot_longer(row_sums, cols = -!!clonotype_col)
    row_sums <- dplyr::group_by(row_sums, !!clonotype_col)
    row_sums <- dplyr::summarize(row_sums, row_sum = sum(.data$value), .groups = "drop")
    row_sums <- row_sums$row_sum

    a     <- length(row_sums[row_sums == 2])      # Intersection
    a_b_c <- length(row_sums[row_sums == 1]) + a  # Union

    jaccard <- a / a_b_c

    res <- tibble::tibble(
      Var1 = comparison[1],
      Var2 = comparison[2],
      jaccard
    )

    res
  }

  # Fetch clonotypes and cell identities
  # so_idents   <- Seurat::Idents(sobj_in)
  so_idents   <- Seurat::FetchData(sobj_in, cluster_col)
  so_idents   <- so_idents[, cluster_col]
  so_idents   <- as.character(so_idents)
  uniq_idents <- unique(so_idents)
  uniq_idents <- stats::na.omit(uniq_idents)

  ctypes   <- Seurat::FetchData(sobj_in, clonotype_col)
  vdj_meta <- dplyr::bind_cols(ctypes, idents = so_idents)
  vdj_meta <- filter(vdj_meta, !is.na(!!dplyr::sym(clonotype_col)))

  # Create data.frame for calculating Jaccard index
  j_df <- dplyr::mutate(vdj_meta, num = 1)
  j_df <- dplyr::mutate(j_df, num = 1)

  j_df <- dplyr::group_by(j_df, .data$idents)
  j_df <- dplyr::group_split(j_df)
  j_df <- purrr::map(
    .x          = j_df,
    .f          = tidyr::pivot_wider,
    names_from  = .data$idents,
    values_from = .data$num,
    values_fn   = list
  )

  j_df <- purrr::map(j_df, tidyr::unnest, cols = -!!dplyr::sym(clonotype_col))
  j_df <- purrr::map(j_df, unique)
  j_df <- purrr::reduce(j_df, dplyr::full_join, by = clonotype_col)
  j_df <- dplyr::mutate_all(j_df, tidyr::replace_na, replace = 0)

  # Create data.frame of comparisons
  # Currently combinations include duplicates
  # comps <- combinations(
  #   n = length(uniq_idents),
  #   r = 2,
  #   v = uniq_idents
  # ) %>%
  #   as.data.frame(stringsAsFactors = F) %>%
  #   dplyr::bind_rows(data.frame(V1 = uniq_idents, V2 = uniq_idents))

  comps <- expand.grid(
    uniq_idents, uniq_idents,
    stringsAsFactors = F
  )

  if (!is.null(ref_cluster)) {
    comps <- data.frame(
      Var1 = ref_cluster,
      Var2 = uniq_idents,
      stringsAsFactors = F
    )
  }

  # Calculate Jaccard index for comparisons
  res <- purrr::map2_dfr(
    .x = comps$Var1,
    .y = comps$Var2,
    .f = ~ calc_jidx(j_df, c(.x, .y), clonotype_col = clonotype_col)
  )

  res <- dplyr::mutate(res, Var1 = stringr::str_c(prefix, .data$Var1, "_jaccard"))
  res <- tidyr::pivot_wider(res, names_from = .data$Var1, values_from = .data$jaccard)

  # Add jaccard index to meta.data
  vdj_meta <- tibble::as_tibble(vdj_meta, rownames = "cell_id")

  vdj_meta <- dplyr::left_join(vdj_meta, res, by = c("idents" = "Var2"))
  vdj_meta <- dplyr::select(vdj_meta, -.data$idents)
  vdj_meta <- tibble::column_to_rownames(vdj_meta, "cell_id")
  vdj_meta <- as.data.frame(vdj_meta)

  res <- Seurat::AddMetaData(sobj_in, vdj_meta)

  res
}


#' Cluster cells based on receptor sequence
#'
#' @param sobj_in Seurat object
#' @param cdr3_col meta.data column containing CDR3 sequences to use for
#' calculating Levenshtein distance
#' @param resolution Clustering resolution to pass to FindClusters
#' @param use_chains Chains to use for calculating Levenshtein distance. If
#' multiple sequences are present for a chain, the first sequence is used.
#' @param prefix Prefix to add to graph name
#' @param ... Additional parameters to pass to FindClusters
#' @return Seurat object with an added shared nearest neighbors graph (vdj_snn)
#' and a meta.data column containing cluster ids
#' @export
cluster_vdj <- function(sobj_in, cdr3_col = "cdr3", resolution = 0.1,
                        use_chains = NULL, prefix = "vdj_", ...) {

  # Extract sequences
  # Only include cells with VDJ data
  orig_idents <- Idents(sobj_in)

  seqs <- Seurat::FetchData(sobj_in, cdr3_col)
  seqs <- tibble::rownames_to_column(seqs, "cell_id")
  seqs <- stats::na.omit(seqs)

  # Select chains to used for calculating distance
  if (!is.null(use_chains)) {
    re       <- stringr::str_c("(?<=", use_chains, ":)[A-Z]+")
    seq_cols <- stringr::str_c("V", seq_along(use_chains))

    seqs <- purrr::map2(re, seq_cols, ~ {
      n_col <- dplyr::sym(.y)
      c_col <- dplyr::sym(cdr3_col)

      dplyr::mutate(seqs, !!n_col := stringr::str_extract(!!c_col, .x))
    })

    seqs <- purrr::reduce(seqs, dplyr::left_join, by = c("cell_id", cdr3_col))
    seqs <- stats::na.omit(seqs)
    seqs <- dplyr::mutate(seqs, seq_c = stringr::str_c(!!!dplyr::syms(seq_cols)))

  } else {
    seqs <- mutate(
      seqs,
      seq_c = stringr::str_extract_all(!!dplyr::sym(cdr3_col), "(?<=:)[A-Z]+"),
      seq_c = purrr::map(.data$seq_c, purrr::reduce, stringr::str_c)
    )
  }

  # Create Levenshtein distance matrix
  seqs     <- purrr::set_names(seqs$seq_c, seqs$cell_id)
  vdj_dist <- utils::adist(seqs)

  # Create nearest neighbors graph
  # Add graph this way or error thrown due to differing number of cells
  vdj_snn  <- Seurat::FindNeighbors(vdj_dist, distance.matrix = T)
  snn_name <- stringr::str_c(prefix, "snn")

  sobj_in@graphs[[snn_name]] <- vdj_snn$snn

  # Find clusters
  res <- Seurat::FindClusters(
    object     = sobj_in,
    resolution = resolution,
    graph.name = snn_name,
    ...
  )

  Idents(res) <- orig_idents

  res
}


#' Run UMAP using Seurat object containing VDJ nearest neighbors graph
#'
#' @param sobj_in Seurat object containing shared nearest neighbors graph for
#' VDJ data
#' @param umap_key Key to use for UMAP columns in meta.data
#' @param vdj_graph Name of shared nearest neighbors graph stored in Seurat
#' object
#' @return Seurat object containing UMAP coordinates in meta.data
#' @export
run_umap_vdj <- function(sobj_in, umap_key = "vdjUMAP_", vdj_graph = "vdj_snn") {

  # Subset sobj_in to only include VDJ cells and add vdj_snn graph
  # RunUMAP does not like running with a graph that does not include results
  # for all cells in the object
  vdj_cells <- rownames(sobj_in[[vdj_graph]])

  vdj_so <- subset(sobj_in, cells = vdj_cells)
  vdj_so[[vdj_graph]] <- sobj_in[[vdj_graph]]

  # Run UMAP and add reduction object back to original object
  vdj_so <- Seurat::RunUMAP(
    object         = vdj_so,
    reduction.name = "vdj_umap",
    reduction.key  = umap_key,
    graph          = vdj_graph
  )

  umap_coords <- Seurat::Embeddings(vdj_so, reduction = "vdj_umap")
  umap_cols   <- stringr::str_c(umap_key, c("1", "2"))

  res <- Seurat::AddMetaData(
    object   = sobj_in,
    metadata = umap_coords,
    col.name = umap_cols
  )

  res
}


#' Subset Seurat object based on VDJ meta.data
#'
#' @param sobj_in Seurat object containing CDR3 sequences
#' @param filt Expression to use for filtering object
#' @param new_col Instead of filtering object create a new column with values
#' based on the filtering expression
#' @param true Value to include in new_col when the filtering expression
#' evaluates to TRUE
#' @param false Value to include in new_col when the filtering expression
#' evaluates to FALSE
#' @param clonotype_col meta.data column containing clonotype IDs. This column
#' is used to determine which cells have VDJ data. If the clonotype_col is set
#' to NULL, filtering is performed regardless of whether VDJ data is present
#' for the cell.
#' @param split_cols The names of meta.data columns to split into vectors. This
#' allows for filtering based on multiple terms present in the column.
#' @param split_sep Separator to use for splitting columns provided by the
#' split_cols argument
#' @param return_seurat Return a Seurat object. If set to FALSE, the meta.data
#' table is returned
#' @return Filtered Seurat object
#' @export
filter_vdj <- function(sobj_in, filt, new_col = NULL, true = TRUE, false = FALSE, clonotype_col = "clonotype_id",
                       split_cols = c("chain", "cdr3"), split_sep = ";", return_seurat = T) {

  meta_df <- tibble::as_tibble(sobj_in@meta.data, rownames = ".cell_id")

  # Split columns into vectors
  if (!is.null(split_cols)) {

    # Save original columns
    split_names <- purrr::set_names(
      x  = split_cols,
      nm = stringr::str_c(".", split_cols)
    )

    meta_df <- dplyr::mutate(meta_df, !!!dplyr::syms(split_names))

    # Split into vectors
    meta_df <- dplyr::mutate(meta_df, dplyr::across(
      dplyr::all_of(split_cols),
      ~ str_split(.x, split_sep)
    ))

    meta_df <- tidyr::unnest(meta_df, cols = dplyr::all_of(split_cols))

    # Convert to numeric if possible
    meta_df <- dplyr::mutate(meta_df, across(
      dplyr::all_of(split_names),
      ~ ifelse(!is.na(suppressWarnings(as.numeric(.x))), as.numeric(.x), .x)
    ))

    meta_df <- dplyr::group_by(meta_df, .data$.cell_id)
  }

  # Store results from filtering expression
  meta_df <- dplyr::mutate(meta_df, .KEEP = {{filt}})

  # Add new column with values based on filtering expression
  if (!is.null(new_col)) {
    meta_df <- dplyr::mutate(
      meta_df,
      !!dplyr::sym(new_col) := dplyr::if_else(.data$.KEEP, true = true, false = false)
    )

    if (!is.null(clonotype_col)) {
      meta_df <- dplyr::mutate(
        meta_df,
        !!dplyr::sym(new_col) := ifelse(
          is.na(!!dplyr::sym(clonotype_col)),
          NA,
          !!dplyr::sym(new_col)
        )
      )
    }

  } else {
    if (!is.null(clonotype_col)) {
      meta_df <- dplyr::mutate(
        meta_df,
        .KEEP = dplyr::if_else(is.na(!!dplyr::sym(clonotype_col)), TRUE, .data$.KEEP)
      )
    }

    meta_df <- dplyr::filter(meta_df, .data$.KEEP)
  }

  # Remove columns created for filtering
  if (!is.null(split_cols)) {
    split_names <- purrr::set_names(
      x  = names(split_names),
      nm = unname(split_names)
    )

    meta_df <- dplyr::select(meta_df, !dplyr::all_of(c(split_cols, ".KEEP")))
    meta_df <- dplyr::rename(meta_df, !!!dplyr::syms(split_names))
    meta_df <- dplyr::distinct(meta_df)
    meta_df <- dplyr::ungroup(meta_df)
  }

  # Add meta.data to Seurat object
  meta_df <- tibble::column_to_rownames(meta_df, ".cell_id")

  if (!return_seurat) {
    return(meta_df)
  }

  cells <- rownames(meta_df)
  res   <- subset(sobj_in, cells = cells)
  res   <- Seurat::AddMetaData(res, meta_df)

  res
}

#' Calculate gene usage
#'
#' @param sobj_in Seurat object containing VDJ data
#' @param gene_col meta.data column containing genes used for each clonotype
#' @param cluster_col meta.data column containing cell clusters to use when
#' calculating gene usage
#' @param n_genes Number of top genes to include in output
#' @return tibble containing gene usage summary
#' @export
calc_usage <- function(sobj_in, gene_col, cluster_col = "orig.ident",
                       n_genes = NULL) {

  # data.frame to calculate usage
  vdj_cols <- c("cell_id", gene_col, cluster_col)

  meta_data <- sobj_in@meta.data
  meta_data <- as_tibble(meta_data, rownames = "cell_id")
  meta_data <- dplyr::select(meta_data, dplyr::all_of(vdj_cols))
  meta_data <- dplyr::filter(meta_data, !is.na(!!dplyr::sym(gene_col)))
  meta_data <- dplyr::group_by(meta_data, !!dplyr::sym(cluster_col))
  meta_data <- dplyr::mutate(
    meta_data,
    !!dplyr::sym(gene_col) := stringr::str_split(!!dplyr::sym(gene_col), ";"),
    n_cells = n_distinct(.data$cell_id)
  )
  meta_data <- dplyr::ungroup(meta_data)

  # All genes used
  vdj_genes <- pull(meta_data, gene_col)
  vdj_genes <- unlist(vdj_genes)
  vdj_genes <- unique(vdj_genes)
  vdj_genes <- sort(vdj_genes)

  # Create data.frame with gene usage
  res <- map(vdj_genes, ~ {
    gene_name <- .x

    res <- mutate(
      meta_data,
      used = purrr::map_lgl(!!dplyr::sym(gene_col), ~ gene_name %in% .x)
    )

    res <- dplyr::group_by(res, !!dplyr::sym(cluster_col), n_cells)
    res <- dplyr::tally(res, used)
    res <- dplyr::ungroup(res)
    res <- dplyr::mutate(res, !!dplyr::sym(gene_name) := n / n_cells)
    res <- dplyr::select(res, -n, -n_cells)

    res
  })

  res <- purrr::reduce(res, dplyr::left_join, by = cluster_col)
  res <- tidyr::pivot_longer(
    res,
    cols      = c(-!!dplyr::sym(cluster_col)),
    names_to  = gene_col,
    values_to = "usage"
  )
  res <- dplyr::group_by(res, !!dplyr::sym(gene_col))
  res <- dplyr::mutate(res, ave_usage = mean(.data$usage))
  res <- dplyr::ungroup(res)

  # Filter for top used genes
  if (!is.null(n_genes)) {
    res <- dplyr::top_n(res, n = n_genes, wt = .data$ave_usage)
  }

  res
}



