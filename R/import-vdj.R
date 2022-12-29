#' Import V(D)J data
#'
#' @param input Object containing single cell data, if set to NULL a data.frame
#' containing V(D)J results will be returned
#' @param vdj_dir Directory containing the output from cellranger vdj. A vector
#' or named vector can be given to load data from multiple runs. If a named
#' vector is given, the cell barcodes will be prefixed with the provided names.
#' This mimics the behavior of Seurat::Read10X().
#' @param prefix Prefix to add to new columns
#' @param data_cols Additional columns from filtered_contig_annotations.csv to
#' include in object.
#' @param filter_chains Only include chains with at least one productive and
#' full length contig.
#' @param filter_paired Only include clonotypes with paired chains. For TCR
#' data each clonotype must have at least one TRA and TRB chain, for BCR data
#' each clonotype must have at least one IGH chain and at least one IGK or IGL
#' chain.
#' @param define_clonotypes Define clonotype IDs based on V(D)J data. This is
#' useful if the V(D)J datasets being loaded do not have consistent clonotype
#' IDs, i.e., clonotype1 is not the same across samples. Possible values are:
#'
#' - 'cdr3aa', define clonotypes based on the CDR3 amino acid sequence
#' - 'cdr3nt', define clonotypes based on the CDR3 nucleotide sequence
#' - 'cdr3_gene', define clonotypes based on the combination of the CDR3
#' nucleotide sequence and the V(D)J genes.
#'
#' When defining clonotypes, only productive full length chains will be used.
#' Set to NULL (default) to use the clonotype IDs already present in the input
#' data.
#'
#' @param include_mutations Include information about the number of
#' insertions/deletions/mismatches for each chain. This requires the
#' concat_ref.bam file from cellranger vdj to be present the directory provided
#' to vdj_dir. If include_mutations is TRUE, filter_chains is also
#' automatically set TRUE since indel data is only available for productive
#' chains.
#' @param aggr_dir Path to cellranger aggr output. To include mutation
#' information for each chain, also provide paths to the original cellranger
#' vdj output directories using the vdj_dir argument.
#'
#' To correctly match cell barcodes to those in the object, gene expression
#' data for each sample must be loaded in the same order as the samples were
#' specified in the cellranger aggr config file. In addition, if loading
#' mutation data, sample paths provided to the vdj_dir argument must also be in
#' the same order as the samples were specified in the cellranger aggr config
#' file.
#'
#' @param quiet If `TRUE` progress updates will not be displayed
#' @param sep Separator to use for storing per cell V(D)J data
#' @return Single cell object or data.frame with added V(D)J data
#'
#' @examples
#' # Loading multiple datasets
#' # to ensure cell barcodes for the V(D)J data match those in the object
#' # load the datasets in the same order as the gene expression data
#' data_dir <- system.file("extdata/splen", package = "djvdj")
#'
#' vdj_dirs <- c(
#'   file.path(data_dir, "BL6_BCR"),
#'   file.path(data_dir, "MD4_BCR")
#' )
#'
#' res <- splen_so |>
#'   import_vdj(vdj_dir = vdj_dirs)
#'
#' head(slot(res, "meta.data"), 1)
#'
#' # Specifying cell prefixes using vector names
#' # cell barcode prefixes can also be specified by passing a named vector
#' vdj_dirs <- c(
#'   BL6 = file.path(data_dir, "BL6_BCR"),
#'   MD4 = file.path(data_dir, "MD4_BCR")
#' )
#'
#' res <- splen_so |>
#'   import_vdj(vdj_dir = vdj_dirs)
#'
#' head(slot(res, "meta.data"), 1)
#'
#' # Only include V(D)J data for paired chains
#' res <- splen_so |>
#'   import_vdj(
#'     vdj_dir = vdj_dirs,
#'     filter_paired = TRUE
#'   )
#'
#' head(slot(res, "meta.data"), 1)
#'
#' # Defining clonotypes
#' # this is useful if the original clonotype IDs are not consistent across
#' # datasets, i.e. clonotype1 is not the same for all samples
#' res <- splen_so |>
#'   import_vdj(
#'     vdj_dir = vdj_dirs,
#'     define_clonotypes = "cdr3_gene"
#'   )
#'
#' head(slot(res, "meta.data"), 1)
#'
#' # Include mutation information for each chain
#' # this information will be included if the file concat_ref.bam is present
#' # including mutation information will cause data import to be slower
#' res <- splen_so |>
#'   import_vdj(
#'     vdj_dir = vdj_dirs,
#'     include_mutations = TRUE
#'   )
#'
#' head(slot(res, "meta.data"), 1)
#'
#' @export
import_vdj <- function(input = NULL, vdj_dir = NULL, prefix = "",
                       data_cols = NULL, filter_chains = TRUE,
                       filter_paired = FALSE, define_clonotypes = NULL,
                       include_mutations = FALSE, aggr_dir = NULL,
                       quiet = FALSE, sep = ";") {

  # Set global variables based on prefix
  global$chain_col     <- paste0(prefix, "chains")
  global$clonotype_col <- paste0(prefix, "clonotype_id")
  global$sep           <- sep

  # Check input classes
  .check_args(
    environment(), data_cols = list(len_one = FALSE, allow_null = TRUE)
  )

  # Check input values
  # vdj_dir or aggr_dir must be provided
  load_aggr <- !is.null(aggr_dir)

  if (is.null(vdj_dir) && !load_aggr) {
    cli::cli_abort("`vdj_dir` or `aggr_dir` must be provided")
  }

  # Check that vdj_dir is also provided when loading mutation data for
  # aggr results
  if (load_aggr && is.null(vdj_dir) && include_mutations) {
    cli::cli_warn(
      "To include V(D)J mutation data when loading cellranger aggr results,
       paths to the original cellranger vdj output directories must be
       provided to the `vdj_dir` argument"
    )

    include_mutations <- FALSE
  }

  # When including indel data, only use productive full length chains
  if (!filter_chains && include_mutations) {
    filter_chains <- TRUE

    cli::cli_warn(
      "When `include_mutations` is `TRUE`, `filter_chains` is also
       automatically set `TRUE` since mutation data is only available for
       productive chains"
    )
  }

  # Sequence columns to include
  # lengths will be calculated for these columns
  # by default only include CDR3 sequences unless user specifies others
  # LIST AA COLUMN FIRST
  seq_cols  <- c("fwr1", "cdr1", "fwr2", "cdr2", "fwr3", "cdr3", "fwr4")
  seq_cols  <- purrr::map(seq_cols, ~ c(.x, paste0(.x, "_nt")))
  seq_cols  <- purrr::reduce(seq_cols, c)
  cdr3_cols <- grep("^cdr3", seq_cols, value = TRUE)
  seq_cols  <- seq_cols[seq_cols %in% c(cdr3_cols, data_cols)]
  data_cols <- data_cols[!data_cols %in% seq_cols]

  # V(D)J columns to include
  gene_cols  <- c("v_gene", "d_gene", "j_gene", "c_gene")
  count_cols <- c("reads", "umis")
  qc_cols    <- c("productive", "full_length")
  len_cols   <- paste0(seq_cols, "_length")

  # Columns containing per-cell info
  # cell_cols <- c("barcode", "clonotype_id")
  cell_cols <- c("barcode", "clonotype_id", "paired")

  # Optional aggr columns
  aggr_cols <- c("donor", "origin")

  # Columns containing per-chain info that needs to be collapsed for each cell
  # user provided columns are included here
  sep_cols <- c(
    gene_cols, "chains",
    seq_cols,  count_cols,
    qc_cols
  )

  data_cols <- data_cols[!data_cols %in% c(sep_cols, cell_cols, len_cols)]
  sep_cols  <- c(sep_cols, data_cols)

  # Set cell barcode prefixes
  # if input object is provided, must match barcodes
  .add_progress_step("Loading V(D)J data", quiet = quiet)

  if (!is.null(input)) {
    bcs <- .get_meta(input)[[global$cell_col]]

    prfx_df <- .extract_cell_prefix(bcs, strip_bcs = FALSE)
    prfx_df <- dplyr::distinct(prfx_df, .data$prfx, .data$sfx)

    prfxs <- prfx_df$prfx
    sfxs  <- prfx_df$sfx

    if (!is.null(names(vdj_dir))) {
      vdj_prfxs <- names(vdj_dir) <- paste0(names(vdj_dir), "_")

      if (any(duplicated(prfxs))) {
        cli::cli_abort(
          "To match the provided cell prefixes ({vdj_prfxs}) with those
           in the object ({prfxs}), the cell prefixes in the object
           cannot be duplicated"
        )
      }

      if (!all(names(vdj_dir) %in% prfxs)) {
        cli::cli_abort(
          "The provided cell prefixes ({vdj_prfxs}) do not match
           those in the input object ({prfxs})"
        )
      }

      sfxs  <- sfxs[match(names(vdj_dir), prfxs)]
      prfxs <- names(vdj_dir)
    }

  # If no prefixes, auto-generate, do not add prefix if only one sample
  # Read10X() will add the prefix, "1_", "2_", "3_", etc. for each sample
  } else if (!is.null(vdj_dir)) {
    prfxs <- names(vdj_dir)

    if (is.null(prfxs)) {
      prfxs <- ""

      if (length(vdj_dir) > 1) prfxs <- paste0(seq_along(vdj_dir), "_")
    }

    sfxs <- rep("-1", length(vdj_dir))
  }

  # Load V(D)J data and add cell prefixes
  if (!is.null(aggr_dir)) {
    cell_cols <- c(cell_cols, aggr_cols)

    contigs <- .load_aggr_data(aggr_dir, prfxs, sfxs)
    contigs <- list(contigs)

  } else {
    contigs <- .load_vdj_data(vdj_dir, prfxs, sfxs)
  }

  # vdj_cols should have all columns that should be included in output
  vdj_cols <- c(cell_cols, sep_cols)

  # For genes replace NAs
  # if a chain is missing a V(D)J segment, the gene name will be left empty
  # when read into R this results in an NA
  contigs <- purrr::map(contigs, ~ {
    dplyr::mutate(.x, across(all_of(gene_cols), tidyr::replace_na, "None"))
  })

  # Filter for productive full length chains
  if (filter_chains) {
    contigs <- purrr::map(contigs, dplyr::filter, !!!syms(qc_cols))
  }

  # Add indel info for each contig
  # if indel data is included, always filter for productive contigs since most
  # non-productive contigs are missing indel data
  if (include_mutations) {
    .add_progress_step("Calculating mutation frequencies", quiet = quiet)

    # Fix contig_ids in contigs
    contigs <- purrr::map(contigs, ~ {
      mutate(
        .x,
        contig_sfx = unlist(.str_extract_all(.data$contig_id, "_contig_[0-9]+$")),
        contig_id  = paste0(.data$barcode, .data$contig_sfx),
        contig_sfx = NULL
      )
    })

    # Load mutation data
    indels <- .load_muts(vdj_dir, prfxs, sfxs)

    if (!is.null(indels)) {
      if (!is.null(aggr_dir)) indels <- list(dplyr::bind_rows(indels))

      indel_cols <- names(indels[[1]])
      indel_cols <- indel_cols[indel_cols != "contig_id"]

      # Join indel data
      # SHOULD CHECK BARCODE OVERLAP HERE!!!
      # IF BARCODES DO NOT OVERLAP HERE, WILL RETURN ALL 0s
      indel_ctigs <- purrr::map2(
        contigs, indels, dplyr::left_join, by = "contig_id"
      )

      # Replace NAs with 0
      # contigs that did not have any mutations will have NAs
      indel_ctigs <- purrr::map(
        indel_ctigs,
        ~ mutate(.x, dplyr::across(all_of(indel_cols), tidyr::replace_na, 0))
      )

      contigs    <- indel_ctigs
      count_cols <- c(count_cols, indel_cols)
      sep_cols   <- c(sep_cols, indel_cols)
      vdj_cols   <- c(vdj_cols, indel_cols)
    }
  }

  # Classify input data as TCR or BCR
  .add_progress_step("Formatting V(D)J data", quiet = quiet)

  vdj_class <- purrr::map_chr(contigs, .classify_vdj)
  vdj_class <- unique(vdj_class)

  if (length(vdj_class) > 1) {
    cli::cli_abort(
      "Multiple data types detected ({vdj_class}), provided data must be
       either TCR or BCR. To add both TCR and BCR data to the same object,
       run {.fn import_vdj} separately for each and use the `prefix` argument to
       add distinct column names."
    )
  }

  # Identify paired chains
  contigs <- purrr::map(contigs, .identify_paired, vdj_class)

  # Calculate cell barcode overlap
  # use map to check each sample separately
  # bind contig data.frames
  overlap_stats <- purrr::imap(contigs, ~ .calc_overlap(input, .x, .y))

  if (all(purrr::map_chr(overlap_stats, ~ .x$Status) == "x")) {
    .print_import_summary(overlap_stats)

    cli::cli_abort(
      "Cell barcodes do not match those in the object,
       this will occur if you are loading the samples in the wrong order or are
       providing the wrong cell barcode prefixes. If loading results
       from cellranger aggr, check that gene expression data for each sample
       was loaded into the object in the same order as the samples were
       specified in the cellranger aggr config file."
    )
  }

  contigs <- dplyr::bind_rows(contigs)

  # Check for 'exact_subclonotype_id' columns, not included in all versions of
  # cellranger
  ex_sub_cols <- identical(vdj_class, "BCR") &&
    "exact_subclonotype_id" %in% colnames(contigs)

  if (ex_sub_cols) {
    cell_cols <- c(cell_cols, "exact_subclonotype_id")
    vdj_cols  <- c(vdj_cols, "exact_subclonotype_id")
  }

  # Calculate sequence lengths
  # report length 0 if there is no reported CDR3 sequence
  contigs <- dplyr::mutate(
    contigs,
    across(
      all_of(seq_cols), ~ ifelse(.x == "None", 0, nchar(.x)),
      .names = "{.col}_length"
    )
  )

  sep_cols <- c(sep_cols, len_cols)
  vdj_cols <- c(vdj_cols, len_cols)

  # Remove contigs that do not have an assigned clonotype_id
  n_remove <- contigs$clonotype_id
  n_remove <- n_remove[is.na(n_remove)]
  n_remove <- length(n_remove)

  if (n_remove > 0) {
    cli::cli_warn(
      "{n_remove} contig{?s} do not have an assigned clonotype_id,
       these contigs will be removed"
    )

    contigs <- dplyr::filter(contigs, !is.na(.data$clonotype_id))
  }

  # Select V(D)J columns to keep
  # check that all vdj_cols are in data
  # some columns could be duplicated if also provided to data_cols argument
  vdj_cols <- unique(vdj_cols)

  .check_obj_cols(contigs, vdj_cols, list_avail = TRUE)

  contigs <- dplyr::select(contigs, all_of(vdj_cols))

  # Check for NAs in data, additional NAs would indicate malformed input
  if (!all(stats::complete.cases(contigs))) {
    cli::cli_abort("Malformed input data, `NA`s are present, check input files")
  }

  # Check if sep is already present in sep_cols
  sep <- .check_sep(contigs, sep_cols, sep)

  # Sum contig reads, UMIs, and mutations for chains since some chains are
  # supported by multiple contigs
  # In the vloupe browser the UMI count is summed, but the summed read count
  # and summed mutations do not always match
  grp_cols <- vdj_cols[!vdj_cols %in% count_cols]
  contigs  <- dplyr::group_by(contigs, !!!syms(grp_cols))

  contigs <- dplyr::summarize(
    contigs, across(all_of(count_cols), sum), .groups = "drop"
  )

  # Filter paired chains
  if (filter_paired) contigs <- dplyr::filter(contigs, .data$paired)

  # Order chains and CDR3 sequences
  # when rows are collapsed, the cdr3 sequences must be in the same order for
  # every cell. This is required so the cdr3 columns can be used directly as
  # the clonotype ID
  contigs <- dplyr::arrange(
    contigs,
    .data$barcode, .data$chains, .data$cdr3_nt
  )

  # Extract isotypes from c_gene for IGH chain (for BCR data only)
  if (vdj_class %in% c("BCR", "Multi")) {
    contigs <- .extract_isotypes(contigs)

    cell_cols <- c(cell_cols, "isotype")
  }

  # Collapse chains into a single row for each cell
  # include columns containing per-cell info groups so they are included in the
  # summarized results
  sep_cols <- sep_cols[!sep_cols %in% cell_cols]
  contigs  <- dplyr::group_by(contigs, !!!syms(cell_cols))

  meta <- summarize(
    contigs,
    across(
      all_of(sep_cols),
      ~ paste0(as.character(.x), collapse = sep)
    ),
    n_chains = n(),
    .groups = "drop"
  )

  # Reorder columns
  meta <- dplyr::relocate(meta, "paired",          .after = "full_length")
  meta <- dplyr::relocate(meta, all_of(len_cols),  .after = last(seq_cols))
  meta <- dplyr::relocate(meta, "n_chains",        .after = "chains")
  meta <- dplyr::relocate(meta, all_of(gene_cols), .after = last(len_cols))

  if (vdj_class %in% c("BCR", "Multi")) {
    meta <- dplyr::relocate(meta, "isotype", .after = last(gene_cols))
  }

  # Check for duplicated cell barcodes
  if (any(duplicated(meta$barcode))) {
    cli::cli_abort(
      "Malformed input data, multiple clonotype_ids
       are associated with the same cell barcode"
    )
  }

  # Allow user to redefine clonotypes
  res <- tibble::column_to_rownames(meta, "barcode")

  if (!is.null(define_clonotypes)) {
    .add_progress_step("Defining clonotypes", quiet = quiet)

    clone_cols <- list(
      cdr3aa    = "cdr3",
      cdr3nt    = "cdr3_nt",
      cdr3_gene = c("cdr3_nt", gene_cols[gene_cols != "c_gene"])
    )

    if (!define_clonotypes %in% names(clone_cols)) {
      cli::cli_abort(
        "`define_clonotypes` must be {.or {names(clone_cols)}}"
      )
    }

    clone_cols <- clone_cols[[define_clonotypes]]

    filt_chains <- NULL

    if (filter_chains) filt_chains <- qc_cols

    res <- define_clonotypes(
      res, data_cols = clone_cols, filter_chains = filt_chains
    )
  }

  # Filter to only include cells with valid clonotype_id
  # cells with missing clonotype have a clonotype_id of 'None'
  res <- dplyr::filter(res, .data$clonotype_id != "None")

  if (nrow(res) == 0) {
    cli::cli_abort("No valid clonotypes present, check input data")
  }

  # Add prefix to V(D)J columns
  res <- dplyr::rename_with(res, ~ paste0(prefix, .x))

  # Add new meta.data to input object
  res <- .merge_meta(input, res)

  cli::cli_progress_done()

  if (!quiet) .print_import_summary(overlap_stats)

  res
}

#' Load V(D)J data
#'
#' @param vdj_dir Directory containing the output from cellranger vdj. A vector
#' or named vector can be given to load data from several runs. If a named
#' vector is given, the cell barcodes will be prefixed with the provided names.
#' This mimics the behavior of the Read10X function found in the Seurat
#' package.
#' @param contig_file cellranger vdj output file containing data for each
#' contig annotation
#' @param chk_none Value of 'None' will be replaced with FALSE for the
#' specified columns and converted to logical
#' @return List containing one data.frame for each path provided to vdj_dir
#' @importFrom readr read_csv cols
#' @noRd
.load_vdj_data <- function(vdj_dir, cell_prfxs, cell_sfxs,
                           contig_file = "filtered_contig_annotations.csv",
                           chk_none = c("productive", "full_length")) {

  col_spec <- readr::cols(
    v_gene = readr::col_character(),
    d_gene = readr::col_character(),
    j_gene = readr::col_character(),
    c_gene = readr::col_character()
  )

  # Check for file and return path
  res <- purrr::map_chr(vdj_dir, .get_vdj_path, file = contig_file)

  # Load data
  res <- purrr::map(res, ~ {
    readr::read_csv(
      .x,
      col_types      = col_spec,
      progress       = FALSE,
      show_col_types = FALSE
    )
  })

  # Replace 'None' in productive with FALSE
  res <- purrr::map(res, ~ {
    d <- dplyr::filter(.x, is_cell)

    d <- .replace_none(d, chk_none)

    d <- dplyr::rename(
      d,
      chains       = chain,
      clonotype_id = raw_clonotype_id
    )

    d
  })

  # Format cell barcode prefixes
  prfx_args <- list(
    df_in      = res,
    cell_prfxs = cell_prfxs,
    cell_sfxs  = cell_sfxs
  )

  res <- purrr::pmap(prfx_args, ~ {
    .format_cell_prefixes(..., bc_col = "barcode")
  })

  res
}

#' Load data from cellranger aggr
#'
#' @param aggr_dir Directory containing the output from cellranger aggr
#' @param contig_file cellranger aggr output file containing data for each
#' contig annotation
#' @param chk_none Value of 'None' will be replaced with FALSE for the
#' specified columns and converted to logical
#' @return data.frame
#' @noRd
.load_aggr_data <- function(aggr_dir, cell_prfxs, cell_sfxs,
                            contig_file = "filtered_contig_annotations.csv",
                            chk_none = c("productive", "full_length")) {

  col_spec <- readr::cols(
    v_gene = readr::col_character(),
    d_gene = readr::col_character(),
    j_gene = readr::col_character(),
    c_gene = readr::col_character()
  )

  # Check for file and return path
  res <- .get_vdj_path(aggr_dir, file = contig_file)

  # Load data
  res <- readr::read_csv(
    res,
    col_types      = col_spec,
    progress       = FALSE,
    show_col_types = FALSE
  )

  # Filter for contigs in cells
  res <- dplyr::filter(res, .data$is_cell)

  # Replace 'None' with FALSE for QC columns
  res <- .replace_none(res, chk_none)

  res <- dplyr::rename(res, chains = "chain", clonotype_id = "raw_clonotype_id")

  # Format cell barcode prefixes
  res <- .format_cell_prefixes(
    res,
    bc_col     = "barcode",
    cell_prfxs = cell_prfxs,
    cell_sfxs  = cell_sfxs
  )

  res
}

#' Format cell barcode prefixes
#'
#' @param df_in data.frame
#' @param bc_col Column containing cell barcodes
#' @param prfxs Named vector containing new cell prefixes
#' @return data.frame with formatted barcodes
#' @importFrom stringr str_remove
#' @noRd
.format_cell_prefixes <- function(df_in, bc_col = "barcode", cell_prfxs,
                                  cell_sfxs) {

  # Extract current cell prefixes
  bcs <- df_in[[bc_col]]

  prfx_df <- .extract_cell_prefix(bcs, strip_bcs = TRUE)

  # Match old and new prefixes
  new <- dplyr::distinct(prfx_df, .data$prfx, .data$sfx)

  if (nrow(new) != length(cell_prfxs)) {
    cli::cli_abort(
      "The number of provided cell prefixes does not match the number of
       unique prefixes present on barcodes"
    )
  }

  new$new_prfx <- cell_prfxs
  new$new_sfx  <- cell_sfxs

  prfx_df <- dplyr::left_join(prfx_df, new, by = c("prfx", "sfx"))

  # Format cell barcodes
  prfx_df <- dplyr::mutate(
    prfx_df,
    prfx = ifelse(is.na(.data$new_prfx), .data$prfx, .data$new_prfx),
    sfx  = ifelse(is.na(.data$new_sfx), .data$sfx, .data$new_sfx),
    bc   = paste0(.data$prfx, .data$bc, .data$sfx)
  )

  df_in[[bc_col]] <- prfx_df$bc

  df_in
}

.extract_cell_prefix <- function(bcs, strip_bcs, bc_len = 16) {
  bc_re  <- paste0("[ATGCN]{", bc_len, "}")
  sep_re <- "[^[:alnum:]]"

  p <- .extract_pattern(bcs, paste0("^.+", sep_re, "(?=", bc_re, ")"))
  s <- .extract_pattern(bcs, paste0("(?<=", bc_re, ")", sep_re, ".+$"))

  res <- tibble::tibble(
    bc   = bcs,
    prfx = p,
    sfx  = s
  )

  # Would be nice to implement base R version of str_remove that accepts a
  # vector of patterns
  if (strip_bcs) {
    res <- dplyr::mutate(
      res,
      bc = stringr::str_remove(.data$bc, paste0("^", .data$prfx)),
      bc = stringr::str_remove(.data$bc, paste0(.data$sfx, "$"))
    )
  }

  res
}

.extract_pattern <- function(x, pattern) {
  res <- .str_extract_all(x, pattern)
  res <- map_chr(res, ~ ifelse(purrr::is_empty(.x), "", .x))

  res
}

#' Replace 'None' with FALSE
#'
#' @param df_in data.frame
#' @param clmns Columns to replace 'None' and convert to logical
#' @return data.frame
#' @importFrom stringr str_replace
#' @noRd
.replace_none <- function(df_in, clmns) {

  clmns <- clmns[!map_lgl(df_in[clmns], is.logical)]

  if (purrr::is_empty(clmns)) return(df_in)

  res <- dplyr::mutate(
    df_in,
    dplyr::across(all_of(clmns), ~ {
      as.logical(stringr::str_replace(.x, "^None$", "FALSE"))
    })
  )

  res
}

#' Load mutation information for each contig
#'
#' @param vdj_dir Directory containing the output from cellranger vdj. A vector
#' or named vector can be given to load data from several runs. If a named
#' vector is given, the cell barcodes will be prefixed with the provided names.
#' This mimics the behavior of the Read10X function found in the Seurat
#' package.
#' @param bam_file bam file from cellranger vdj containing alignment data
#' comparing each contig with the germline reference
#' @return List containing one data.frame for each path provided to vdj_dir
#' @importFrom Rsamtools scanBam
#' @noRd
.load_muts <- function(vdj_dir, cell_prfxs, cell_sfxs,
                       bam_file  = "concat_ref.bam",
                       airr_file = "airr_rearrangement.tsv") {

  # Retrieve bam and airr file paths
  file_paths <- c(bam = bam_file, airr = airr_file)

  file_paths <- purrr::map(file_paths, ~ {
    fl <- .x

    purrr::map_chr(vdj_dir, .get_vdj_path, file = fl, warn = TRUE)
  })

  any_missing <- any(purrr::map_lgl(file_paths, ~ any(is.na(.x))))

  if (any_missing) {
    cli::cli_warn(
      "To add mutation data to object {bam_file} and {airr_file}
       must be present for all samples, check that these files are in the
       provided directory paths, mutation data not added to object"
    )

    return(NULL)
  }

  # Extract mutations from bam file
  mut_coords <- purrr::map(file_paths$bam, .extract_mut_coords)

  # Extract VDJ coords from AIRR
  vdj_coords <- purrr::map(file_paths$airr, .extract_vdj_coords)

  # Map mutations to VDJ segments
  res <- purrr::map2(mut_coords, vdj_coords, .map_muts)

  # Extract cell barcode from contig_id
  id_re <- "^.+(?=_contig_[0-9]+$)"

  res <- purrr::map(
    res,
    mutate,
    barcode = unlist(.str_extract_all(.data$contig_id, id_re))
  )

  # Format cell barcode prefixes
  prfx_args <- list(
    df_in      = res,
    cell_prfxs = cell_prfxs,
    cell_sfxs  = cell_sfxs
  )

  res <- purrr::pmap(prfx_args, .format_cell_prefixes, bc_col = "barcode")

  res <- purrr::map(
    res,
    mutate,
    contig_sfx = unlist(.str_extract_all(.data$contig_id, "_contig_[0-9]+$")),
    contig_id  = paste0(.data$barcode, .data$contig_sfx),
    contig_sfx = NULL,
    barcode    = NULL
  )

  res
}

.extract_mut_coords <- function(bam_file) {

  bam_info <- Rsamtools::scanBam(bam_file)[[1]]

  wdths <- as.data.frame(bam_info$seq@ranges)$width

  bam_info <- tibble::tibble(
    cigar     = bam_info$cigar,
    contig_id = bam_info$qname,
    len       = wdths
  )

  bam_info <- dplyr::filter(bam_info, grepl("_contig_[0-9]+$", .data$contig_id))

  # Get 0-based coordinates for mutations
  # set width of deletion coordinates as 0
  res <- dplyr::mutate(
    bam_info,
    n    = .str_extract_all(.data$cigar, "[0-9]+(?=[^0-9])"),
    type = .str_extract_all(.data$cigar, "(?<=[0-9])[^0-9]{1}")
  )

  res <- tidyr::unnest(res, all_of(c("n", "type")))
  res <- dplyr::group_by(res, .data$contig_id)

  res <- dplyr::mutate(
    res,
    n     = as.numeric(.data$n),
    idx   = ifelse(.data$type != "D", .data$n, 0),
    end   = cumsum(.data$idx),
    start = lag(.data$end, default = 0)
  )

  res <- dplyr::ungroup(res)
  res <- dplyr::filter(res, .data$type != "=")

  res <- dplyr::select(
    res,
    all_of(c("contig_id", "len", "start", "end", "type", "n"))
  )

  res
}

.extract_vdj_coords <- function(airr_file) {

  col_spec <- readr::cols(
    v_call  = readr::col_character(),
    v_cigar = readr::col_character(),
    d_call  = readr::col_character(),
    d_cigar = readr::col_character(),
    j_call  = readr::col_character(),
    j_cigar = readr::col_character(),
    c_call  = readr::col_character(),
    c_cigar = readr::col_character(),
    v_sequence_start = readr::col_double(),
    v_sequence_end   = readr::col_double(),
    d_sequence_start = readr::col_double(),
    d_sequence_end   = readr::col_double(),
    j_sequence_start = readr::col_double(),
    j_sequence_end   = readr::col_double(),
    c_sequence_start = readr::col_double(),
    c_sequence_end   = readr::col_double()
  )

  airr <- readr::read_tsv(
    airr_file,
    col_types      = col_spec,
    progress       = FALSE,
    show_col_types = FALSE
  )

  # Pull V(D)J gene coordinates from AIRR file
  # tidyr::extract is much faster than tidyr::separate
  coord_cols_re <- "^([vdjc])(?=_).*(?<=_)(start|end)$"

  res <- dplyr::select(
    airr,
    contig_id = "sequence_id",
    dplyr::matches(coord_cols_re, perl = TRUE)
  )

  if (ncol(res) == 1) {
    cli::cli_abort("V(D)J coordinates not found, check {.file airr_file}")
  }

  res <- tidyr::pivot_longer(res, -"contig_id")
  res <- dplyr::filter(res, !is.na(.data$value))
  res <- tidyr::extract(res, "name", c("seg", "pos"), coord_cols_re)
  res <- tidyr::pivot_wider(res, names_from = "pos")

  res <- dplyr::mutate(res,
    start = .data$start - 1, len = .data$end - .data$start
  )

  res <- dplyr::select(
    res,
    all_of(c("contig_id", "len", "start", "end", "seg"))
  )

  res
}

.map_muts <- function(mut_coords, vdj_coords) {

  mut_key <- c(I = "ins", D = "del", X = "mis")

  mut_coords <- dplyr::mutate(
    mut_coords,
    type = dplyr::recode(.data$type, !!!mut_key)
  )

  # If no vdj_coords, return mutation totals
  if (identical(vdj_coords, NA)) {
    res <- tidyr::pivot_wider(
      all_muts,
      names_from  = "type", values_from = "n", values_fill = 0
    )

    res <- dplyr::mutate(
      res,
      across(starts_with("all_"), ~ .x / .data$len, .names = "{.col}_freq")
    )

    return(res)
  }

  # Intersect mutations with VDJ gene coordinates for each contig
  # some annotations overlap each other! Example: AAACCTGAGAACTGTA-1_contig_1
  # left_join + mutate is much faster than valr::bed_intersect, probably due
  # to the extreme number of "chromosomes"
  vdj_muts <- dplyr::left_join(
    mut_coords, vdj_coords, by = "contig_id", suffix = c("", ".seg")
  )

  vdj_muts <- dplyr::filter(
    vdj_muts, .data$start < .data$end.seg & .data$end > .data$start.seg
  )

  vdj_muts <- dplyr::mutate(
    vdj_muts,
    len = .data$len.seg,

    new_start = ifelse(
      .data$start >= .data$start.seg, .data$start, .data$start.seg
    ),

    new_end = ifelse(
      .data$end <= .data$end.seg, .data$end, .data$end.seg
    ),

    new_end = ifelse(
      .data$type == mut_key[["D"]], .data$new_end + 1, .data$new_end
    ),

    n = ifelse(
      .data$type != mut_key[["D"]], .data$new_end - .data$new_start, .data$n
    )
  )

  # Identify junction indels
  jxn_muts <- filter(vdj_muts, .data$type %in% unname(mut_key[c("I", "D")]))

  jxn_muts <- mutate(
    jxn_muts,
    seg = case_when(
      .data$seg == "v" & .data$end.seg   == .data$new_end   ~ "vd",
      .data$seg == "d" & .data$start.seg == .data$new_start ~ "vd",
      .data$seg == "d" & .data$end.seg   == .data$new_end   ~ "dj",
      .data$seg == "j" & .data$start.seg == .data$new_start ~ "dj",
      TRUE ~ as.character(NA)
    )
  )

  jxn_muts <- dplyr::filter(jxn_muts, !is.na(.data$seg))
  jxn_muts <- dplyr::select(jxn_muts, -"len")

  vdj_muts <- bind_rows(vdj_muts, jxn_muts)

  # Summarize mutation counts
  vdj_muts <- dplyr::group_by(
    vdj_muts, .data$contig_id, .data$len, .data$type, .data$seg
  )

  vdj_muts <- dplyr::summarize(vdj_muts, n = sum(.data$n), .groups = "drop")

  # Summarize total mutations and total length per contig
  # for each mutation type, sum total for v, d, j, and c segments, exclude jxns
  all_muts <- dplyr::filter(vdj_muts, !.data$seg %in% c("vd", "dj"))
  all_muts <- dplyr::group_by(all_muts, .data$contig_id, .data$type)

  all_muts <- dplyr::summarize(
    all_muts,
    n       = sum(.data$n),
    len     = sum(.data$len),
    seg     = "all",
    .groups = "drop"
  )

  vdj_muts <- dplyr::bind_rows(vdj_muts, all_muts)
  res      <- tidyr::unite(
    vdj_muts, "type", all_of(c("seg", "type")), sep = "_"
  )

  # Set final output columns
  freq_cols <- mut_cols <- c("v", "d", "j", "c", "all")
  jxn_cols  <- c("vd", "dj")

  mut_cols <- purrr::map(mut_cols, paste0, "_", mut_key)
  mut_cols <- purrr::reduce(mut_cols, c)

  jxn_cols <- purrr::map(jxn_cols, paste0, "_", unname(mut_key[c("I", "D")]))
  jxn_cols <- purrr::reduce(jxn_cols, c)
  mut_cols <- c(mut_cols, jxn_cols)

  freq_cols <- purrr::map_chr(freq_cols, paste0, "_", mut_key[["X"]])

  # Calculate mismatch frequency
  freq <- dplyr::filter(res, .data$type %in% freq_cols)

  freq <- dplyr::mutate(
    freq,
    n    = round(.data$n / .data$len, 6),
    type = paste0(.data$type, "_freq"),
    len  = NULL
  )

  res <- dplyr::bind_rows(res, freq)
  res <- dplyr::select(res, -"len")

  res <- tidyr::pivot_wider(
    res,
    names_from  = "type",
    values_from = "n",
    values_fill = 0
  )

  # Add 0s for missing columns and set column order
  # these are segments with no mutations for any chain
  mut_cols <- c(mut_cols, paste0(freq_cols, "_freq"))

  missing_cols <- mut_cols[!mut_cols %in% names(res)]

  res[, missing_cols] <- 0

  res <- res[, c("contig_id", mut_cols)]

  res
}

#' Check for V(D)J data file in provided directory
#'
#' @param vdj_dir Directory containing the output from cellranger vdj
#' @param file Name of cellranger vdj output file
#' @param warn When the file is not found display a warning message instead of
#' an error
#' @return path to cellranger vdj output file
#' @noRd
.get_vdj_path <- function(vdj_dir, file, warn = FALSE) {

  path <- file.path(vdj_dir, file)

  if (!file.exists(path)) path <- paste0(path, ".gz")
  if (!file.exists(path)) path <- NA

  if (is.na(path)) {
    fn <- cli::cli_abort

    if (warn) fn <- cli::cli_warn

    if (!file.exists(vdj_dir)) {
      fn("{.file {vdj_dir}} does not exist")

    } else {
      fn("{file} not found in {.file {vdj_dir}}")
    }
  }

  path
}

#' Check for separator in data.frame
#'
#' @param df_in data.frame
#' @param sep_cols Names of columns to check for sep, if `NULL` all columns
#' will be checked
#' @param sep Separator to use for storing V(D)J data
#' @return Separator with white space stripped
#' @noRd
.check_sep <- function(df_in, sep_cols, sep) {
  if (is.null(sep_cols))  sep_cols <- colnames(df_in)
  if (is.null(sep))       return(sep)
  if (!is.character(sep)) cli::cli_abort("`sep` must be a character")

  # Strip whitespace from sep
  sep <- gsub("[[:space:]]", "", sep)

  has_sep <- grepl(sep, df_in[, sep_cols, drop = FALSE], fixed = TRUE)

  if (any(has_sep)) {
    cli::cli_abort(
      "The string '{sep}' is already present in the input data,
       select a different value for `sep`"
    )
  }

  sep
}

#' Determine whether TCR or BCR data were provided
#'
#' @param df_in data.frame containing V(D)J data formatted so that each row
#' represents a single contig
#' @param chain_col Column in input data containing chain identity
#' @return Character string indicating whether TCR or BCR data were provided
#' @noRd
.classify_vdj <- function(df_in, chain_col = "chains") {

  chains <- list(
    "TCR" = c("TRA", "TRB", "TRD", "TRG"),
    "BCR" = c("IGH", "IGK", "IGL")
  )

  n_chains <- purrr::imap(chains, ~ purrr::set_names(rep(.y, length(.x)), .x))
  n_chains <- purrr::flatten(n_chains)

  # Classify chains
  # remove values that do not match, such as chains with "None"
  n_chains <- n_chains[df_in[[chain_col]]]
  n_chains <- n_chains[!is.na(names(n_chains))]
  n_chains <- table(as.character(n_chains))

  # Error if no chains match
  if (is_empty(n_chains)) {
    chains <- unlist(chains, use.names = FALSE)

    cli::cli_abort(
      "None of the expected chains ({.or {chains}}) were found,
       unable to determine whether TCR or BCR data were provided"
    )
  }

  # Calculate fraction of BCR/TCR chains
  # set type if >50% match
  res <- n_chains / sum(n_chains)
  res <- names(res[res > 0.5])

  if (purrr::is_empty(res)) {
    res   <- "Multi"
    n_bcr <- n_chains[["BCR"]]
    n_tcr <- n_chains[["TCR"]]

    cli::cli_warn(
      "Equal number of BCR ({n_bcr}) and TCR ({n_tcr}) chains detected, unable
       to determine data type"
    )
  }

  res
}

#' Check cell barcode overlap with object
#'
#' @param input Single cell object
#' @param meta meta.data to check against object
#' @param nm Sample name to use for messages
#' @param pct_min Warn user if the percent overlap is less than pct_min
#' @return input data
#' @noRd
.calc_overlap <- function(input, meta, nm, pct_min = 25) {

  met_dat <- dplyr::distinct(meta, .data$barcode, .data$paired)

  met_cells   <- met_dat$barcode
  n_met_cells <- length(met_cells)
  n_met_pair  <- length(met_cells[met_dat$paired])

  if (is.null(input)) {
    n_obj_cells <- n_overlap <- pct_overlap <- NA

  } else {
    obj_meta    <- .get_meta(input)
    obj_cells   <- obj_meta[[global$cell_col]]
    n_obj_cells <- length(obj_cells)
    n_overlap   <- length(obj_cells[obj_cells %in% met_cells])
    pct_overlap <- round(n_overlap / n_met_cells, 0) * 100
  }

  status <- dplyr::case_when(
    n_overlap == 0        ~ "x",
    pct_overlap < pct_min ~ "!",
    TRUE                  ~ "v"
  )

  res <- list(
    "Status"   = status,
    "Sample"   = nm,
    "# cells"  = n_obj_cells,
    "# VDJ"    = n_met_cells,
    "# paired" = n_met_pair,
    "Overlap"  = n_overlap,
    "Percent"  = pct_overlap
  )

  res
}

.print_import_summary <- function(stats) {

  stats <- purrr::map(stats, ~ map(.x, ~ {
    if (is.na(.x)) .x <- "NA"
    .x
  }))

  stats <- map(stats, ~ {
    names(.x)[names(.x) == "Sample"] <- " "
    .x
  })

  # Calculate maximum char width for header and values in each column
  # exclude sample from header
  clmn_wdth <- dplyr::bind_rows(stats)
  clmn_wdth <- imap(clmn_wdth, ~ max(nchar(c(.x, .y))))

  nms <- names(clmn_wdth)
  nms <- nms[!nms %in% c("Status", "Percent")]

  # Format header
  header <- purrr::map2(nms, clmn_wdth[nms], .add_padding)
  header <- paste0(header, collapse = " | ")
  header <- paste0("\u00a0\u00a0", header, " |")

  cli::cli_rule()
  cli::cli_text(header)

  # Format rows
  purrr::walk(stats, ~ {
    rw  <- .x[nms]
    pct <- .x$Percent

    if (!identical(pct, "NA")) pct <- paste0(pct, "%")

    padded <- purrr::map2(rw, clmn_wdth[names(rw)], .add_padding)

    res <- paste0(padded, collapse = " | ")
    res <- paste0(res, " | ", cli::col_blue(pct))

    names(res) <- .x$Status

    cli::cli_bullets(res)
  })

  cli::cli_rule()
}

.add_padding <- function(x, n) {
  n_pad <- n - nchar(x)

  pad <- paste0(rep("\u00a0", n_pad), collapse = "")

  if (is.numeric(x) || identical(x, "NA")) {
    res <- paste0(pad, x)

  } else {
    res <- paste0(x, pad)
  }

  res
}


#' Identify clonotypes with paired chains
#'
#' @param df_in data.frame containing V(D)J data formatted so each row
#' represents a single contig
#' @return Input data.frame with paired column added
#' @noRd
.identify_paired <- function(df_in, vdj_class) {

  res <- dplyr::group_by(df_in, .data$barcode)

  if (identical(vdj_class, "TCR")) {
    res <- dplyr::mutate(
      res,
      paired = all(c("TRA", "TRB") %in% .data$chains)
    )

  } else if (identical(vdj_class, "BCR")) {
    res <- dplyr::mutate(
      res,
      paired = "IGH" %in% .data$chains & any(c("IGL", "IGK") %in% .data$chains)
    )

  } else {
    res <- dplyr::mutate(res, paired = FALSE)
  }

  res <- dplyr::ungroup(res)

  res
}

#' Add isotypes to V(D)J data
#'
#' @param df_in data.frame containing V(D)J data formatted so each row
#' represents a single contig
#' @param iso_col Column containing data to use for extracting isotypes
#' @param chain_col Column in input data containing chain identity
#' @return Input data.frame with isotype column added
#' @noRd
.extract_isotypes <- function(df_in, iso_col = "c_gene", chain_col = "chains") {

  # Pull data for isotypes
  isos <- df_in[[iso_col]]

  chains <- df_in[[chain_col]]

  idx <- chains == "IGH" & isos != "None"

  isos[idx] <- substr(isos[idx], 1, 4)

  isos[!idx] <- as.character(NA)

  # Identify cells with multiple isotypes
  iso_df <- df_in[, c("barcode", iso_col)]

  iso_df[iso_col] <- isos

  iso_df <- dplyr::distinct(iso_df, .data$barcode, .data$c_gene)
  iso_df <- stats::na.omit(iso_df)

  dups <- iso_df$barcode
  dups <- dups[duplicated(dups)]

  # Add isotypes to meta.data
  iso_df <- mutate(
    iso_df,
    isotype = ifelse(.data$barcode %in% dups, "Multi", !!sym(iso_col))
  )

  isos <- purrr::set_names(
    iso_df$isotype,
    iso_df$barcode
  )

  res <- mutate(
    df_in,
    isotype = unname(isos[.data$barcode]),
    isotype = tidyr::replace_na(.data$isotype, "None")
  )

  res
}

#' Add cli progress step
#'
#' @param msg Message for progress step
#' @param quiet If `TRUE` do nothing, if `FALSE` add progress step
#' @param envir Environment to set progress step
#' @param ... Additional arguments to pass to cli::cli_progress_step()
#' @noRd
.add_progress_step <- function(msg, quiet = FALSE, envir = parent.frame(),
                               ...) {

  if (!quiet) cli::cli_progress_step(msg, .envir = envir, ...)
}

#' Define clonotypes based on V(D)J data
#'
#' This will assign new clonotype IDs based on the combination of values
#' present in the provided columns
#'
#' @param input Single cell object or data.frame containing V(D)J data. If a
#' data.frame is provided, the cell barcodes should be stored as row names.
#' @param data_cols meta.data columns containing V(D)J data to use for defining
#' clonotypes
#' @param clonotype_col Name of column to use for storing clonotype IDs
#' @param filter_chains Column(s) to use for filtering chains prior to defining
#' clonotypes (e.g. productive, full_length). The column(s) must contain TRUE
#' or FALSE for each chain. If NULL, all chains are used when defining
#' clonotypes.
#' @param sep Separator used for storing per cell V(D)J data
#' @return Single cell object or data.frame with added clonotype IDs
#'
#' @examples
#' # Define clonotypes using the CDR3 nucleotide sequence
#' res <- define_clonotypes(
#'   vdj_so,
#'   data_cols = "cdr3_nt"
#' )
#'
#' head(slot(res, "meta.data"), 1)
#'
#' # Define clonotypes based on the combination of the CDR3 nucleotide sequence
#' # and the V and J genes
#' res <- define_clonotypes(
#'   vdj_sce,
#'   data_cols = c("cdr3_nt", "v_gene", "j_gene")
#' )
#'
#' head(slot(res, "colData"), 1)
#'
#' # Modify the name of the column used to store clonotype IDs
#' res <- define_clonotypes(
#'   vdj_so,
#'   data_cols = "cdr3_nt",
#'   clonotype_col = "NEW_clonotype_id"
#' )
#'
#' head(slot(res, "meta.data"), 1)
#'
#' # When defining clonotypes only use chains that are productive
#' res <- define_clonotypes(
#'   vdj_sce,
#'   data_cols = "cdr3_nt",
#'   filter_chains = "productive"
#' )
#'
#' head(slot(res, "colData"), 1)
#'
#' @export
define_clonotypes <- function(input, data_cols, clonotype_col = "clonotype_id",
                              filter_chains = c("productive", "full_length"),
                              sep = global$sep) {

  # Check that columns are present in object
  .check_obj_cols(input, data_cols, filter_chains)

  # Check input classes
  .check_args(
    environment(), filter_chains = list(Class = "character", len_one = FALSE)
  )

  # Get meta.data
  # overwrite exising clonotype_col if it has the same name
  meta <- .get_meta(input)

  meta <- dplyr::select(meta, -any_of(clonotype_col))

  all_cols <- c(data_cols, filter_chains)

  # Remove cells with NAs for any data_cols
  vdj <- dplyr::filter(
    meta,
    dplyr::if_all(dplyr::all_of(all_cols), ~ !is.na(.x))
  )

  vdj <- dplyr::select(vdj, all_of(c(global$cell_col, all_cols)))

  # Only use values in data_cols that are TRUE for all filter_chains columns
  # first identify contigs TRUE for all filter_chains columns
  # subset each data_cols column based on .clone_idx
  if (!is.null(filter_chains)) {
    clmns <- syms(filter_chains)

    vdj <- tibble::column_to_rownames(vdj, global$cell_col)

    vdj <- mutate_vdj(
      vdj,
      .clone_idx = list(purrr::reduce(list(!!!clmns), ~ .x & .y)),

      dplyr::across(
        dplyr::all_of(data_cols),
        ~ paste0(.x[.data$.clone_idx], collapse = ""),
        .names = ".clone_{.col}"
      ),

      clonotype_col = NULL,
      data_cols = all_cols
    )

    data_cols <- paste0(".clone_", data_cols)

    vdj <- .get_meta(vdj)
  }

  # Add new clonotype IDs
  vdj <- dplyr::mutate(
    vdj,
    .new_clone = paste(!!!syms(data_cols), sep = ""),
    .new_id    = rank(.data$.new_clone, ties.method = "min"),

    !!sym(clonotype_col) := ifelse(
      .data$.new_clone == "", "None", paste0("clonotype", .data$.new_id)
    )
  )

  # Add new clonotype IDs to meta.data
  vdj  <- dplyr::select(vdj, all_of(c(global$cell_col, clonotype_col)))
  meta <- dplyr::left_join(meta, vdj, by = global$cell_col)

  res <- .add_meta(input, meta)

  res
}
