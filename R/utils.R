#' Add VDJ data to Seurat object
#'
#' @param sobj_in Seurat object
#' @param vdj_dir cellranger vdj output directory
#' @param prefix Prefix to add to new meta.data columns
#' @param cell_prefix Prefix to add to cell barcodes
#' @param filter_contigs Only include chains with at least one productive full
#' length contig
#' @return Seurat object with VDJ data added to meta.data
#' @export
import_vdj <- function(sobj_in, vdj_dir, prefix = "", cell_prefix = "",
                       filter_contigs = TRUE) {

  # Load contigs
  count_cols <- c("reads", "umis")

  split_cols <- c(
    "v_gene",     "d_gene",
    "j_gene",     "c_gene",
    "chain",      "cdr3",
    "cdr3_nt",    count_cols,
    "productive", "full_length"
  )

  vdj_cols <- c(
    "barcode", "raw_clonotype_id",
    split_cols
  )

  vdj_dir <- file.path(vdj_dir, "filtered_contig_annotations.csv")
  contigs <- readr::read_csv(vdj_dir)

  # Extract chain, add cell barcode prefix
  pat <- "^[A-Z]{3}"

  contigs <- dplyr::mutate(
    contigs,
    barcode = stringr::str_c(cell_prefix, .data$barcode),
    chain   = case_when(
      str_detect(.data$v_gene, pat) ~ str_extract(.data$v_gene, pat),
      str_detect(.data$d_gene, pat) ~ str_extract(.data$d_gene, pat),
      str_detect(.data$j_gene, pat) ~ str_extract(.data$j_gene, pat),
      str_detect(.data$c_gene, pat) ~ str_extract(.data$c_gene, pat),
      TRUE ~ "None"
    )
  )

  # Filter for productive full length contigs
  if (filter_contigs) {
    contigs <- dplyr::filter(contigs, .data$productive, .data$full_length)
  }

  contigs <- dplyr::select(contigs, all_of(vdj_cols))

  # Sum contig reads and UMIs for clonotype chains
  grp_cols <- vdj_cols[!vdj_cols %in% count_cols]
  contigs  <- dplyr::group_by(contigs, !!!syms(grp_cols))
  contigs  <- dplyr::summarize(
    contigs,
    across(all_of(count_cols), sum),
    .groups = "drop"
  )

  # Merge rows for each cell
  contigs <- dplyr::arrange(contigs, .data$barcode, .data$raw_clonotype_id, .data$chain)
  contigs <- dplyr::group_by(contigs, .data$barcode, .data$raw_clonotype_id)

  meta_df <- summarize(
    contigs,
    n_chains = n(),
    dplyr::across(
      all_of(split_cols),
      ~ stringr::str_c(as.character(.x), collapse = ";")
    ),
    .groups = "drop"
  )

  # Filter for cells present in sobj_in
  cells   <- Seurat::Cells(sobj_in)
  meta_df <- dplyr::filter(meta_df, .data$barcode %in% cells)
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


#' Subset Seurat object based on VDJ meta.data
#'
#' @param sobj_in Seurat object containing VDJ data
#' @param filt Condition to use for filtering object
#' @param new_col Instead of filtering object create a new column with values
#' based on filtering condition
#' @param true Value to use for new_col when filtering condition evaluates to
#' TRUE
#' @param false Value to use for new_col when filtering condition evaluates to
#' FALSE
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
filter_vdj <- function(sobj_in, filt, new_col = NULL, true = TRUE, false = FALSE,
                       clonotype_col = "clonotype_id", split_cols = c("chain", "cdr3"),
                       split_sep = ";", return_seurat = TRUE) {

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

    # Convert to numeric or logical if possible
    convert_type <- function(x, fun) {
      suppressWarnings(ifelse(!is.na(fun(x)) | is.na(x), fun(x), x))
    }

    meta_df <- dplyr::mutate(
      meta_df,
      across(
        dplyr::all_of(split_cols),
        ~ convert_type(.x, as.numeric)
      ),
      across(
        dplyr::all_of(split_cols),
        ~ convert_type(.x, as.logical)
      )
    )

    meta_df <- dplyr::group_by(meta_df, .data$.cell_id)
  }

  # Store results from filtering expression
  meta_df <- dplyr::mutate(meta_df, .KEEP = {{filt}})

  # Add new column with values based on filtering expression
  if (!is.null(new_col)) {
    meta_df <- dplyr::mutate(
      meta_df,
      !!dplyr::sym(new_col) := ifelse(
        .data$.KEEP,
        yes = {{true}},
        no  = {{false}})
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


#' Calculate clonotype abundance
#'
#' @param sobj_in Seurat object containing VDJ data
#' @param clonotype_col meta.data column containing clonotype IDs
#' @param cluster_col meta.data column containing cluster IDs to use for
#' grouping cells when calculating clonotype abundance
#' @param prefix Prefix to add to new meta.data columns
#' @return ggplot object
#' @export
calc_abundance <- function(sobj_in, clonotype_col = "clonotype_id",
                           cluster_col = NULL, prefix = "") {

  meta_df <- sobj_in@meta.data
  meta_df <- tibble::as_tibble(meta_df, rownames = ".cell_id")
  meta_df <- dplyr::filter(meta_df, !is.na(!!dplyr::sym(clonotype_col)))

  if (!is.null(cluster_col)) {
    meta_df <- dplyr::group_by(meta_df, !!sym(cluster_col))
  }

  freq_col  <- stringr::str_c(prefix, "clone_freq")
  abund_col <- stringr::str_c(prefix, "clone_abund")

  meta_df <- dplyr::mutate(
    meta_df,
    .n_cells = dplyr::n_distinct(.data$.cell_id)
  )

  meta_df <- dplyr::group_by(
    meta_df,
    !!dplyr::sym(clonotype_col),
    add = TRUE
  )

  meta_df <- dplyr::mutate(
    meta_df,
    !!sym(freq_col)  := dplyr::n_distinct(.data$.cell_id),
    !!sym(abund_col) := !!dplyr::sym(freq_col) / .n_cells
  )

  meta_df <- dplyr::select(meta_df, -.n_cells)
  meta_df <- tibble::column_to_rownames(meta_df, ".cell_id")

  res <- Seurat::AddMetaData(sobj_in, metadata = meta_df)

  res
}


#' Calculate receptor diversity (inverse Simpson Index)
#'
#' @param sobj_in Seurat object
#' @param clonotype_col meta.data column containing clonotype ids
#' @param cluster_col meta.data column containing cluster IDs used for
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
  meta_df       <- tibble::as_tibble(sobj_in@meta.data, rownames = ".cell_id")
  vdj_df        <- dplyr::filter(meta_df, !is.na(!!clonotype_col))

  # Count clonotypes
  if (!is.null(cluster_col)) {
    vdj_cols    <- c(cluster_col, vdj_cols)
    cluster_col <- dplyr::sym(cluster_col)
    vdj_df      <- dplyr::group_by(vdj_df, !!cluster_col)
  }

  vdj_df <- dplyr::group_by(vdj_df, !!clonotype_col, .add = TRUE)

  vdj_df <- dplyr::summarize(
    vdj_df,
    .n      = dplyr::n_distinct(.data$.cell_id),
    .groups = "drop"
  )

  # Calculate diversity
  if (!is.null(cluster_col)) {
    vdj_df <- dplyr::group_by(vdj_df, !!cluster_col)
  }

  div_col <- stringr::str_c(prefix, "diversity")

  vdj_df <- dplyr::mutate(
    vdj_df,
    frac     = .data$.n / sum(.data$.n),
    sum_frac = sum(.data$frac ^ 2),

    !!dplyr::sym(div_col) := 1 - .data$sum_frac
    # !!dplyr::sym(div_col) := 1 / sum_frac
  )

  vdj_df <- dplyr::select(vdj_df, all_of(c(vdj_cols, div_col)))

  # Add resuts to meta.data
  meta_df <- dplyr::left_join(meta_df, vdj_df, by = vdj_cols)
  meta_df <- tibble::column_to_rownames(meta_df, ".cell_id")
  meta_df <- as.data.frame(meta_df)

  res <- Seurat::AddMetaData(sobj_in, metadata = meta_df)

  res
}


#' Calculate Jaccard index for cell identities
#'
#' @param sobj_in Seurat object
#' @param clonotype_col meta.data column containing clonotype ids
#' @param cluster_col meta.data column containing cell clusters
#' @param prefix Prefix to add to new meta.data columns
#' @param return_seurat Return a Seurat object. If set to FALSE, a matrix is
#' returned
#' @return Seurat object with Jaccard index added to meta.data
#' @export
calc_jaccard <- function(sobj_in, clonotype_col = "clonotype_id", cluster_col,
                         prefix = "", return_seurat = TRUE) {

  # Helper to calculate jaccard index
  calc_jidx <- function(df_in, comparison, ctype_col) {

    if (length(comparison) != 2) {
      stop("comparison must be a character vector with two column names")
    }

    if (length(unique(df_in[, ctype_col])) != length(df_in[, ctype_col])) {
      stop("duplicate clonotype IDs are present")
    }

    uniq_vars <- unique(comparison)

    if (length(uniq_vars) == 1) {
      dup_var       <- stringr::str_c(uniq_vars, "_dup")
      comparison[3] <- dup_var
      dup_var       <- dplyr::sym(dup_var)
      uniq_vars     <- dplyr::sym(uniq_vars)

      df_in <- dplyr::mutate(df_in, !!dup_var := !!uniq_vars)
    }

    row_sums  <- dplyr::select(df_in, dplyr::all_of(c(ctype_col, comparison)))
    ctype_col <- dplyr::sym(ctype_col)
    row_sums  <- tidyr::pivot_longer(row_sums, cols = -!!ctype_col)
    row_sums  <- dplyr::group_by(row_sums, !!ctype_col)

    row_sums <- dplyr::summarize(
      row_sums,
      row_sum = sum(.data$value),
      .groups = "drop"
    )

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
  vdj_meta        <- Seurat::FetchData(sobj_in, c(clonotype_col, cluster_col))
  vdj_meta        <- dplyr::rename(vdj_meta, idents = cluster_col)
  vdj_meta$idents <- as.character(vdj_meta$idents)

  uniq_idents <- unique(vdj_meta$idents)
  uniq_idents <- stats::na.omit(uniq_idents)
  uniq_idents <- sort(uniq_idents)

  # Create data.frame for calculating Jaccard index
  # Use 1's and 0's to indicate if a clonotype is present in a cluster
  vdj_meta <- dplyr::filter(vdj_meta, !is.na(!!dplyr::sym(clonotype_col)))
  j_df     <- dplyr::mutate(vdj_meta, .n = 1)

  j_df <- tidyr::pivot_wider(
    j_df,
    names_from  = "idents",
    values_from = ".n",
    values_fn   = list
  )

  j_df <- dplyr::mutate(j_df, dplyr::across(
    dplyr::all_of(uniq_idents),
    ~ purrr::map(.x, unique)
  ))

  j_df <- tidyr::unnest(j_df, cols = dplyr::all_of(uniq_idents))

  j_df <- dplyr::mutate(j_df, dplyr::across(
    dplyr::all_of(uniq_idents),
    ~ replace_na(.x, 0)
  ))

  j_df <- dplyr::arrange(j_df, !!sym(clonotype_col))
  j_df <- select(j_df, all_of(c(clonotype_col, uniq_idents)))

  # Create data.frame of comparisons
  comps <- expand.grid(
    uniq_idents, uniq_idents,
    stringsAsFactors = F
  )

  comps <- dplyr::mutate(comps, comp = map2(Var1, Var2, ~ sort(c(.x, .y))))
  comps <- dplyr::select(comps, .data$comp)
  comps <- unique(comps)

  Var1 <- purrr::map_chr(comps$comp, ~ .x[1])
  Var2 <- purrr::map_chr(comps$comp, ~ .x[2])

  # Calculate Jaccard index for comparisons
  res <- purrr::map2_dfr(
    Var1, Var2,
    ~ calc_jidx(j_df, c(.x, .y), ctype_col = clonotype_col)
  )

  res_i <- dplyr::rename(res, Var1 = Var2, Var2 = Var1)
  res   <- dplyr::bind_rows(res, res_i)
  res   <- unique(res)

  res <- dplyr::mutate(
    res,
    Var1 = stringr::str_c(prefix, .data$Var1, "_jaccard")
  )

  res <- tidyr::pivot_wider(
    res,
    names_from  = .data$Var1,
    values_from = .data$jaccard
  )

  # Return matrix
  if (!return_seurat) {
    res <- dplyr::mutate(
      res,
      Var2 = stringr::str_c(prefix, .data$Var2, "_jaccard")
    )

    res <- tibble::column_to_rownames(res, "Var2")
    res <- as.matrix(res)

    return(res)
  }

  # Add jaccard index to meta.data
  vdj_meta <- tibble::as_tibble(vdj_meta, rownames = ".cell_id")
  vdj_meta <- dplyr::left_join(vdj_meta, res, by = c("idents" = "Var2"))
  vdj_meta <- dplyr::select(vdj_meta, -.data$idents)
  vdj_meta <- tibble::column_to_rownames(vdj_meta, ".cell_id")
  vdj_meta <- as.data.frame(vdj_meta)

  res <- Seurat::AddMetaData(sobj_in, vdj_meta)

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
  vdj_cols <- c(".cell_id", gene_col, cluster_col)

  meta_df <- sobj_in@meta.data
  meta_df <- as_tibble(meta_df, rownames = ".cell_id")
  meta_df <- dplyr::select(meta_df, dplyr::all_of(vdj_cols))
  meta_df <- dplyr::filter(meta_df, !is.na(!!dplyr::sym(gene_col)))
  meta_df <- dplyr::group_by(meta_df, !!dplyr::sym(cluster_col))
  meta_df <- dplyr::mutate(
    meta_df,
    !!dplyr::sym(gene_col) := stringr::str_split(!!dplyr::sym(gene_col), ";"),
    n_cells = n_distinct(.data$.cell_id)
  )
  meta_df <- dplyr::ungroup(meta_df)

  # All genes used
  vdj_genes <- pull(meta_df, gene_col)
  vdj_genes <- unlist(vdj_genes)
  vdj_genes <- unique(vdj_genes)
  vdj_genes <- sort(vdj_genes)

  # Create data.frame with gene usage
  res <- map(vdj_genes, ~ {
    gene_name <- .x

    res <- mutate(
      meta_df,
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



