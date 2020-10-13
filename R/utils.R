#' Add V(D)J data to Seurat object
#'
#' @param sobj_in Seurat object
#' @param vdj_dir cellranger vdj output directory
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
      ~ stringr::str_c(as.character(.x), collapse = sep)
    ),
    .groups = "drop"
  )

  # Filter for cells present in sobj_in
  cells   <- Seurat::Cells(sobj_in)
  meta_df <- dplyr::filter(meta_df, .data$barcode %in% cells)

  meta_df <- dplyr::rename(
    meta_df,
    clonotype_id = .data$raw_clonotype_id,
    chains = .data$chain
  )

  # Add meta.data to Seurat object
  meta_df <- tibble::column_to_rownames(meta_df, "barcode")
  meta_df <- dplyr::rename_with(meta_df, ~ stringr::str_c(prefix, .x))

  res <- Seurat::AddMetaData(sobj_in, metadata = meta_df)

  res
}


#' Subset Seurat object based on V(D)J meta.data
#'
#' @param sobj_in Seurat object containing V(D)J data
#' @param filt Condition to use for filtering object. If set to NULL, all cells
#' will be returned. An argument must be passed to filt and/or new_col.
#' @param new_col Instead of filtering object create a new column with values
#' based on filtering condition. If filt is set to NULL, the values in new_col
#' will be based on the expression passed to true.
#' @param true Expression to generate values for new_col when filtering
#' condition evaluates to TRUE
#' @param false Expression to generate values for new_col when filtering
#' condition evaluates to FALSE
#' @param clonotype_col meta.data column containing clonotype IDs. This column
#' is used to determine which cells have V(D)J data. If clonotype_col is set to
#' NULL, filtering is performed regardless of whether V(D)J data is present for
#' the cell.
#' @param vdj_cols The names of meta.data columns to expand for filtering
#' @param sep Separator to use for expanding meta.data columns
#' @param return_seurat Return a Seurat object. If set to FALSE, the meta.data
#' table is returned.
#' @return Filtered Seurat object
#' @export
filter_vdj <- function(sobj_in, filt = NULL, new_col = NULL, true = TRUE, false = FALSE,
                       clonotype_col = "clonotype_id", vdj_cols = c("chains", "cdr3"),
                       sep = ";", return_seurat = TRUE) {

  if (is.null(dplyr::enexpr(filt)) && is.null(new_col)) {
    stop("filt and/or new_col are required")
  }

  meta_df <- tibble::as_tibble(sobj_in@meta.data, rownames = ".cell_id")

  # Split columns into vectors
  if (!is.null(vdj_cols)) {

    # Save original columns
    split_names <- purrr::set_names(
      x  = vdj_cols,
      nm = stringr::str_c(".", vdj_cols)
    )

    meta_df <- dplyr::mutate(meta_df, !!!dplyr::syms(split_names))

    # Split into vectors
    meta_df <- dplyr::mutate(meta_df, dplyr::across(
      dplyr::all_of(vdj_cols),
      ~ stringr::str_split(.x, sep)
    ))

    meta_df <- tidyr::unnest(meta_df, cols = dplyr::all_of(vdj_cols))

    # Convert to numeric or logical if possible
    meta_df <- dplyr::mutate(
      meta_df,
      across(
        dplyr::all_of(vdj_cols),
        ~ convert_char(.x, as.numeric)
      ),
      across(
        dplyr::all_of(vdj_cols),
        ~ convert_char(.x, as.logical)
      )
    )

    meta_df <- dplyr::group_by(meta_df, .data$.cell_id)
  }

  # Store results from filtering expression
  meta_df <- dplyr::mutate(meta_df, .KEEP = TRUE)

  if (!is.null(dplyr::enexpr(filt))) {
    meta_df <- dplyr::mutate(meta_df, .KEEP = {{filt}})
  }

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
        .KEEP = dplyr::if_else(
          is.na(!!dplyr::sym(clonotype_col)),
          TRUE,
          .data$.KEEP
        )
      )
    }

    meta_df <- dplyr::filter(meta_df, .data$.KEEP)
  }

  # Remove columns created for filtering
  if (!is.null(vdj_cols)) {
    split_names <- purrr::set_names(
      x  = names(split_names),
      nm = unname(split_names)
    )

    meta_df <- dplyr::select(meta_df, !dplyr::all_of(c(vdj_cols, ".KEEP")))
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
#' @param sobj_in Seurat object containing V(D)J data
#' @param clonotype_col meta.data column containing clonotype IDs
#' @param cluster_col meta.data column containing cluster IDs to use for
#' grouping cells when calculating clonotype abundance
#' @param prefix Prefix to add to new meta.data columns
#' @param return_seurat Return a Seurat object. If set to FALSE, the meta.data
#' table is returned.
#' @return Seurat object with clonotype abundance added to meta.data
#' @export
calc_abundance <- function(sobj_in, clonotype_col = "clonotype_id", cluster_col = NULL,
                           prefix = "", return_seurat = TRUE) {

  meta_df <- sobj_in@meta.data
  meta_df <- tibble::as_tibble(meta_df, rownames = ".cell_id")
  meta_df <- dplyr::filter(meta_df, !is.na(!!dplyr::sym(clonotype_col)))

  if (!is.null(cluster_col)) {
    meta_df <- dplyr::group_by(meta_df, !!sym(cluster_col))
  }

  freq_col  <- stringr::str_c(prefix, "clone_freq")
  pct_col <- stringr::str_c(prefix, "clone_pct")

  meta_df <- dplyr::mutate(
    meta_df,
    .n_cells = dplyr::n_distinct(.data$.cell_id)
  )

  meta_df <- dplyr::group_by(
    meta_df,
    !!dplyr::sym(clonotype_col),
    .add = TRUE
  )

  meta_df <- dplyr::mutate(
    meta_df,
    !!sym(freq_col) := dplyr::n_distinct(.data$.cell_id),
    !!sym(pct_col)  := (!!dplyr::sym(freq_col) / .data$.n_cells) * 100
  )

  meta_df <- dplyr::select(meta_df, -.data$.n_cells)
  meta_df <- tibble::column_to_rownames(meta_df, ".cell_id")

  if (!return_seurat) {
    return(meta_df)
  }

  res <- Seurat::AddMetaData(sobj_in, metadata = meta_df)

  res
}


#' Calculate receptor diversity (inverse Simpson Index)
#'
#' @param sobj_in Seurat object
#' @param clonotype_col meta.data column containing clonotype IDs to use for
#' calculating diversity
#' @param cluster_col meta.data column containing cluster IDs to use for
#' calculating diversity. If cluster_col is omitted, diversity index will be
#' calculated for all clonotypes.
#' @param prefix Prefix to add to new meta.data columns
#' @param return_seurat Return a Seurat object. If set to FALSE, the meta.data
#' table is returned.
#' @return Seurat object with inverse Simpson index added to meta.data
#' @export
calc_diversity <- function(sobj_in, clonotype_col = "clonotype_id", cluster_col = NULL,
                           prefix = "", return_seurat = TRUE) {

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
    # !!dplyr::sym(div_col) := 1 / .data$sum_frac
  )

  vdj_df <- dplyr::select(vdj_df, all_of(c(vdj_cols, div_col)))

  # Add results to meta.data
  meta_df <- dplyr::left_join(meta_df, vdj_df, by = vdj_cols)
  meta_df <- tibble::column_to_rownames(meta_df, ".cell_id")
  meta_df <- as.data.frame(meta_df)

  if (!return_seurat) {
    return(meta_df)
  }

  res <- Seurat::AddMetaData(sobj_in, metadata = meta_df)

  res
}


#' Calculate repertoire overlap between clusters
#'
#' @param sobj_in Seurat object
#' @param clonotype_col meta.data column containing clonotype IDs to use for
#' calculating overlap
#' @param cluster_col meta.data column containing cluster IDs to use for
#' calculating overlap
#' @param prefix Prefix to add to new meta.data columns
#' @param return_seurat Return a Seurat object. If set to FALSE, a matrix is
#' returned
#' @return Seurat object with Jaccard index added to meta.data
#' @export
calc_overlap <- function(sobj_in, clonotype_col = "clonotype_id", cluster_col,
                         prefix = "jcrd_", return_seurat = TRUE) {

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
    Var1 = stringr::str_c(prefix, .data$Var1)
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
      Var2 = stringr::str_c(prefix, .data$Var2)
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
#' Gene usage is calculated as the percent of cells that express each V(D)J
#' gene present in gene_col. Cells that lack V(D)J data and have an NA present
#' in gene_col are excluded from this calculation.
#'
#' @param sobj_in Seurat object containing V(D)J data
#' @param gene_col meta.data column containing genes used for each clonotype
#' @param cluster_col meta.data column containing cell clusters to use when
#' calculating gene usage
#' @param chain Chain to use for calculating gene usage. Set to NULL to include
#' all chains.
#' @param chain_col meta.data column containing chains for each cell
#' @param sep Separator to use for expanding gene_col
#' @return data.frame containing gene usage summary
#' @export
calc_usage <- function(sobj_in, gene_col, cluster_col = NULL, chain = NULL,
                       chain_col = NULL, sep = ";") {

  # data.frame to calculate usage
  split_cols <- c(gene_col, chain_col)
  vdj_cols   <- c(".cell_id", cluster_col, split_cols)

  meta_df <- sobj_in@meta.data
  meta_df <- tibble::as_tibble(meta_df, rownames = ".cell_id")
  meta_df <- dplyr::filter(meta_df, !is.na(!!dplyr::sym(gene_col)))  # remove cells with no gene_col data
  meta_df <- dplyr::select(meta_df, dplyr::all_of(vdj_cols))

  res <- dplyr::mutate(meta_df, across(
    all_of(split_cols),
    ~ stringr::str_split(.x, sep)
  ))

  # Filter chains
  if (!is.null(chain)) {
    if (is.null(chain_col)) {
      stop("must specify chain_col")
    }

    res <- dplyr::mutate(res, !!dplyr::sym(gene_col) := purrr::map2(
      !!dplyr::sym(gene_col), !!dplyr::sym(chain_col), ~ {
        .x <- dplyr::if_else(
          any(.y %in% chain),
          list(.x[.y %in% chain]),
          list("None")
        )

        unlist(.x)
      }
    ))

    res <- dplyr::select(res, -all_of(chain_col))
  }

  res <- tidyr::unnest(res, cols = all_of(gene_col))
  res <- dplyr::distinct(res)

  # Count genes used
  res <- dplyr::group_by(res, !!dplyr::sym(gene_col))

  if (!is.null(cluster_col)) {
    res <- dplyr::group_by(res, !!dplyr::sym(cluster_col), .add = TRUE)
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
      values_from = freq,
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
      n_cells = clust_counts[!!dplyr::sym(cluster_col)],
      n_cells = as.integer(n_cells)
    )

    res <- dplyr::relocate(res, n_cells, .before = freq)
  }

  res <- dplyr::mutate(res, pct = (freq / n_cells) * 100)

  res
}


#' Summarize values for chains
#'
#' Summarize values present for each column provided to the data_cols argument.
#' For each cell, the function(s) provided to .fun will be applied to each
#' unique label in chain_col.
#'
#' @param sobj_in Seurat object containing V(D)J data
#' @param data_cols meta.data columns to summarize
#' @param .fun Function to use for summarizing data_cols
#' @param chain_col meta.data column(s) containing labels for each chain
#' expressed in the cell. These labels are used for grouping the summary
#' output. Set chain_col to NULL to group solely based on the cell ID.
#' @param include_cols Additional columns to include in the output data.frame
#' @param sep Separator to use for expanding data_cols and chain_col
#' @return data.frame containing summary results
#' @export
summarize_chains <- function(sobj_in, data_cols = c("umis", "reads"), .fun,
                             chain_col = "chains", include_cols = NULL, sep = ";") {

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
    ~ stringr::str_split(.x, sep)
  ))

  res <- tidyr::unnest(res, cols = dplyr::all_of(c(data_cols, chain_col)))

  # Summarize data_cols for each chain present for the cell
  res <- dplyr::mutate(res, dplyr::across(
    dplyr::all_of(data_cols),
    ~ convert_char(.x, as.numeric)
  ))

  grp_cols <- c(".cell_id", chain_col, include_cols)
  res      <- dplyr::group_by(res, !!!dplyr::syms(grp_cols))

  res <- dplyr::summarize(
    res,
    dplyr::across(dplyr::all_of(data_cols), .fun),
    .groups = "drop"
  )

  res
}


#' Attempt to convert character vector using provided function
#'
#' @param x Character vector to convert
#' @param fun Function to try
#' @return Value converted using fun
convert_char <- function(x, fun) {
  if (!is.character(x)) {
    return(x)
  }

  suppressWarnings(ifelse(!is.na(fun(x)) | is.na(x), fun(x), x))
}







