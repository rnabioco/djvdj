#' Import V(D)J data
#'
#' @param input Object containing single cell data, if set to NULL a data.frame
#' containing V(D)J results will be returned
#' @param vdj_dir Directory containing the output from cellranger vdj. A vector
#' or named vector can be given to load data from multiple runs. If a named
#' vector is given, the cell barcodes will be prefixed with the provided names.
#' This mimics the behavior of Seurat::Read10X(). Cell barcode prefixes can
#' also be provided using the cell_prefix argument.
#' @param prefix Prefix to add to new columns
#' @param cell_prefix Prefix to add to cell barcodes, this is helpful when
#' loading data from multiple runs into a single object. If NULL, cell barcode
#' prefixes are automatically generated in a similar manner as
#' Seurat::Read10X().
#'
#' For the V(D)J data to be successfully added to the object, the cell prefixes
#' must match the prefixes that are already present in the object. If the cell
#' prefixes are incorrect, import_vdj will be unable to assign the V(D)J data
#' to the correct cells.
#'
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
#' - cdr3aa, define clonotypes based on the CDR3 amino acid sequence
#' - cdr3nt, define clonotypes based on the CDR3 nucleotide sequence
#' - cdr3_gene, define clonotypes based on the combination of the CDR3
#' nucleotide sequence and the V(D)J genes.
#'
#' When defining clonotypes, only productive full length chains will be used.
#' Set to NULL (default) to use the clonotype IDs already present in the input
#' data.
#'
#' @param include_indels Include the number of insertions/deletions for each
#' chain. This requires the concat_ref.bam file from cellranger vdj to be
#' present the directory provided to vdj_dir. If include_indels is TRUE,
#' filter_chains is also automatically set TRUE since indel data is only
#' available for productive chains.
#' @param sep Separator to use for storing per cell V(D)J data
#' @return Single cell object or data.frame with added V(D)J data
#'
#' @examples
#' # Loading multiple datasets
#' vdj_dir <- c(
#'   system.file("extdata/bcr_1", package = "djvdj"),
#'   system.file("extdata/bcr_2", package = "djvdj")
#' )
#'
#' vdj_so <- import_vdj(tiny_so, vdj_dir)
#'
#' head(vdj_so@meta.data, 1)
#'
#' # Specifying cell prefixes
#' # if cell prefixes are not specified when loading multiple datasets,
#' # prefixes will be automatically generated in a similar manner as
#' # Seurat::Read10X
#' vdj_so <- import_vdj(
#'   tiny_so,
#'   vdj_dir = vdj_dir,
#'   cell_prefix = c("1", "2")
#' )
#'
#' head(vdj_so@meta.data, 1)
#'
#' # Specifying cell prefixes using vector names
#' # if a named vector is passed, the names will be used as the cell prefixes
#' vdj_dir <- c(
#'   "1" = system.file("extdata/bcr_1", package = "djvdj"),
#'   "2" = system.file("extdata/bcr_2", package = "djvdj")
#' )
#'
#' vdj_so <- import_vdj(tiny_so, vdj_dir)
#'
#' head(vdj_so@meta.data, 1)
#'
#' # Only include V(D)J data for productive full length chains
#' vdj_so <- import_vdj(
#'   tiny_so,
#'   vdj_dir = vdj_dir,
#'   filter_chains = TRUE
#' )
#'
#' head(vdj_so@meta.data, 1)
#'
#' # Only include V(D)J data for cells with paired chains
#' vdj_so <- import_vdj(
#'   tiny_so,
#'   vdj_dir = vdj_dir,
#'   filter_paired = TRUE
#' )
#'
#' head(vdj_so@meta.data, 1)
#'
#' # Defining clonotypes
#' # this is useful if the original clonotype IDs are not consistent across
#' # datasets, i.e. clonotype1 is not the same for all samples
#' vdj_so <- import_vdj(
#'   tiny_so,
#'   vdj_dir = vdj_dir,
#'   define_clonotypes = "cdr3_gene"
#' )
#'
#' head(vdj_so@meta.data, 1)
#'
#' # Omit indel information for each chain
#' # this information will be included if the file concat_ref.bam is present
#' # to speed up data import, omit indel information
#' vdj_so <- import_vdj(
#'   tiny_so,
#'   vdj_dir = vdj_dir,
#'   include_indels = FALSE
#' )
#'
#' head(vdj_so@meta.data, 1)
#'
#' # Loading both BCR and TCR data
#' # load each dataset separately and assign unique column names using the
#' # prefix argument
#' bcr_dir <- system.file("extdata/bcr_1", package = "djvdj")
#' tcr_dir <- system.file("extdata/tcr_1", package = "djvdj")
#'
#' vdj_so <- import_vdj(
#'   tiny_so,
#'   vdj_dir     = bcr_dir,
#'   prefix      = "bcr_",
#'   cell_prefix = "1"
#' )
#'
#' vdj_so <- import_vdj(
#'   vdj_so,
#'   vdj_dir     = tcr_dir,
#'   prefix      = "tcr_",
#'   cell_prefix = "1"
#' )
#'
#' head(vdj_so@meta.data, 1)
#'
#' # Using import_vdj outside of Seurat
#' # SingleCellExperiment objects are also compatible, or if an input object is
#' # omitted, a data.frame containing the V(D)J data will be returned
#' vdj_sce <- import_vdj(tiny_sce, vdj_dir)
#'
#' head(vdj_sce@colData, 1)
#'
#' vdj_df <- import_vdj(vdj_dir = vdj_dir)
#'
#' head(vdj_df, 1)
#'
#' @export
import_vdj <- function(input = NULL, vdj_dir, prefix = "", cell_prefix = NULL, filter_chains = TRUE,
                       filter_paired = FALSE, define_clonotypes = NULL, include_indels = TRUE, sep = ";") {

  # V(D)J columns to include
  cdr3_cols  <- c("cdr3", "cdr3_nt")
  count_cols <- c("reads", "umis")
  gene_cols  <- c("v_gene", "d_gene", "j_gene", "c_gene")
  indel_cols <- c("n_insertion", "n_deletion", "n_mismatch")
  qc_cols    <- c("productive", "full_length")

  sep_cols <- c(
    gene_cols, "chains",
    cdr3_cols, count_cols,
    qc_cols
  )

  vdj_cols <- c(
    "barcode", "clonotype_id",
    sep_cols
  )

  # Format/check cell prefixes
  vdj_dir <- .format_cell_prefixes(vdj_dir, cell_prefix)

  # Load V(D)J data and add cell prefixes
  contigs <- .load_vdj_data(vdj_dir)

  # Add indel info for each contig
  # if indel data is included, always filter for productive contigs since most
  # non-productive contigs are missing indel data
  if (include_indels) {
    indels <- .load_vdj_indels(vdj_dir)

    if (!is.null(indels)) {
      indel_ctigs <- purrr::map(contigs, dplyr::filter, !!!syms(qc_cols))

      indel_ctigs <- purrr::map2(
        indel_ctigs, indels,
        dplyr::left_join,
        by = "contig_id"
      )

      # Check for NAs after merging with contigs
      # do not include indels if any chains are missing data
      missing_indels <- purrr::map_lgl(
        indel_ctigs,
        ~ any(!stats::complete.cases(.x[, indel_cols]))
      )

      if (any(missing_indels)) {
        warning(
          "Some chains are missing indel data, check input files. ",
          "Indel results will not be included in the object."
        )

      } else {
        contigs  <- indel_ctigs
        sep_cols <- c(sep_cols, indel_cols)
        vdj_cols <- c(vdj_cols, indel_cols)

        if (!filter_chains) {
          warning(
            "When include_indels is TRUE, filter_chains is also automatically ",
            "set TRUE since indel data is only available for productive chains."
          )
        }

        # Set filter_chains FALSE since indel_ctigs has already been filtered
        filter_chains <- FALSE
      }
    }
  }

  # Filter for productive contigs
  if (filter_chains) {
    contigs <- purrr::map(
      contigs,
      dplyr::filter,
      .data$productive, .data$full_length
    )
  }

  # Classify input data as TCR or BCR
  vdj_class <- purrr::map_chr(contigs, .classify_vdj)
  vdj_class <- unique(vdj_class)

  if (length(vdj_class) > 1) {
    vdj_class <- paste0(vdj_class, collapse = ", ")

    stop(
      "Multiple data types detected (", vdj_class, "), provided data must be ",
      "either TCR or BCR. To add both TCR and BCR data to the same object, ",
      "run import_vdj separately for each and use the 'prefix' argument to ",
      "add different column names."
    )
  }

  # Check cell barcode overlap
  # use map to check each sample separately
  # give warning for low overlap
  # bind contig data.frames
  contigs <- purrr::imap_dfr(contigs, ~ .check_overlap(input, .x, .y))

  # Remove contigs that do not have an assigned clonotype_id
  n_remove <- contigs$clonotype_id
  n_remove <- n_remove[is.na(n_remove)]
  n_remove <- length(n_remove)

  if (n_remove > 0) {
    warning(n_remove, " contigs do not have an assigned clonotype_id, these contigs will be removed.")

    contigs <- dplyr::filter(contigs, !is.na(.data$clonotype_id))
  }

  # Check for NAs in data, additional NAs would indicate malformed input.
  if (!all(stats::complete.cases(contigs))) {
    stop("Malformed input data, NAs are present, check input files.")
  }

  # Select V(D)J columns to keep
  contigs <- dplyr::select(contigs, all_of(vdj_cols))

  # Check if sep is already present in sep_cols
  sep <- .check_sep(contigs, sep_cols, sep)

  # Sum contig reads and UMIs for chains
  # some chains are supported by multiple contigs
  grp_cols <- vdj_cols[!vdj_cols %in% count_cols]
  contigs  <- dplyr::group_by(contigs, !!!syms(grp_cols))

  contigs <- dplyr::summarize(
    contigs,
    across(all_of(count_cols), sum),
    .groups = "drop"
  )

  # Determine which clonotypes are paired
  contigs <- .identify_paired(contigs)

  if (filter_paired) {
    contigs <- dplyr::filter(contigs, .data$paired)
  }

  # Calculate CDR3 length
  # report length 0 if there is no reported CDR3 sequence
  contigs <- dplyr::mutate(
    contigs,
    across(
      all_of(cdr3_cols),
      ~ ifelse(.x == "None", 0, nchar(.x)),
      .names = "{.col}_length"
    )
  )

  sep_cols <- c(sep_cols, paste0(cdr3_cols, "_length"))

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
    contigs <- dplyr::group_by(contigs, .data$isotype)
  }

  # Collapse chains into a single row for each cell
  # include isotype, clonotype_id, and paired as groups so they are included in
  # the summarized results
  contigs <- dplyr::group_by(
    contigs,
    .data$barcode, .data$clonotype_id, .data$paired,
    .add = TRUE
  )

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
  meta <- dplyr::relocate(meta, .data$paired, .after = .data$full_length)

  if (vdj_class %in% c("BCR", "Multi")) {
    meta <- dplyr::relocate(meta, .data$isotype, .after = .data$chains)
  }

  # Check for duplicated cell barcodes
  if (any(duplicated(meta$barcode))) {
    stop("Malformed input data, multiple clonotype_ids are associated with the same cell barcode.")
  }

  # Allow user to redefine clonotypes
  res <- tibble::column_to_rownames(meta, "barcode")

  if (!is.null(define_clonotypes)) {
    clone_cols <- list(
      cdr3aa    = "cdr3",
      cdr3nt    = "cdr3_nt",
      cdr3_gene = c("cdr3_nt", gene_cols[gene_cols != "c_gene"])
    )

    if (!define_clonotypes %in% names(clone_cols)) {
      stop("define_clonotypes must be one of ", paste0(names(clone_cols), collapse = ", "), ".")
    }

    clone_cols <- clone_cols[[define_clonotypes]]

    res <- define_clonotypes(
      res,
      vdj_cols      = clone_cols,
      filter_chains = qc_cols
    )
  }

  # Filter to only include cells with valid clonotype_id
  # cells with missing clonotype have a clonotype_id of 'None'
  res <- dplyr::filter(res, .data$clonotype_id != "None")

  if (nrow(res) == 0) {
    warning("No valid clonotypes present, check input data.")
  }

  # Add prefix to V(D)J columns
  res <- dplyr::rename_with(res, ~ paste0(prefix, .x))

  # Add new meta.data to input object
  res <- .merge_meta(input, res)

  res
}

#' Format and add cell prefixes
#'
#' @param vdj_dir Directory containing the output from cellranger vdj. A vector
#' or named vector can be given to load data from several runs. If a named
#' vector is given, the cell barcodes will be prefixed with the provided names.
#' This mimics the behavior of the Read10X function found in the Seurat
#' package. Cell barcode prefixes can also be provided using the cell_prefix
#' argument.
#' @param cell_prefix Prefix to add to cell barcodes, this is helpful when
#' loading data from multiple runs into a single object. If set to NULL, cell
#' barcode prefixes will be automatically generated in a similar way as the
#' Read10X function found in the Seurat package.
#' @param sep Separator to use when appending prefixes to cell barcodes, set to
#' NULL to not add a separator
#' @return Paths provided to vdj_dir with cell prefixes added as names
#' @noRd
.format_cell_prefixes <- function(vdj_dir, cell_prefix, sep = "_") {

  res <- vdj_dir

  # If vdj_dir is not named, check cell_prefix
  if (is.null(names(res))) {

    # If no prefixes, auto-generate, do not add prefix if only one sample
    # Read10X() will add the prefix, "1_", "2_", "3_", etc. for each sample
    if (is.null(cell_prefix)) {
      cell_prefix <- ""

      if (length(res) > 1) {
        cell_prefix <- paste0(seq_along(res), "_")
      }
    }

    # Check there is a cell prefix provided for each path
    if (length(res) != length(cell_prefix)) {
      stop("cell_prefix must be the same length as vdj_dir (", length(res),").")
    }

    names(res) <- cell_prefix
  }

  # Check for NAs in cell prefixes
  if (any(is.na(names(res)))) {
    stop("Cell prefixes cannot include NAs.")
  }

  # Check for duplicated cell prefixes
  if (any(duplicated(names(res)))) {
    dups <- duplicated(names(res))
    dups <- names(res)[dups]
    dups <- paste0(dups, collapse = ", ")

    warning("Some cell barcode prefixes are duplicated: ", dups)
  }

  # Add separator if one is not included in cell prefixes
  if (!is.null(sep)) {
    sep_regex <- paste0(sep, "$")

    nms <- names(res) != "" & !grepl(sep_regex, names(res))

    names(res)[nms] <- paste0(names(res)[nms], sep)
  }

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
#' @return List containing one data.frame for each path provided to vdj_dir
#' @importFrom readr read_csv cols
#' @noRd
.load_vdj_data <- function(vdj_dir, contig_file = "filtered_contig_annotations.csv") {

  # Check for file and return path
  res <- purrr::map_chr(vdj_dir, .get_vdj_path, file = contig_file)

  # Load data
  res <- purrr::map(
    res,
    readr::read_csv,
    col_types = readr::cols(),
    progress  = FALSE
  )

  # Add cell prefixes and replace 'None' in productive with FALSE
  res <- purrr::imap(res, ~ {
    d <- dplyr::filter(.x, is_cell)
    d <- dplyr::rowwise(d)

    d <- dplyr::mutate(
      d,
      dplyr::across(.cols = c(full_length, productive), ~ {
        ifelse(
          !is.logical(.x),
          as.logical(gsub("None", "FALSE", .x)),
          .x
        )
      }),

      barcode = paste0(.y, .data$barcode)
    )

    d <- dplyr::ungroup(d)

    d <- dplyr::rename(
      d,
      chains       = .data$chain,
      clonotype_id = .data$raw_clonotype_id
    )

    d
  })

  res
}

#' Load insertion/deletion information for each contig
#'
#' @param vdj_dir Directory containing the output from cellranger vdj. A vector
#' or named vector can be given to load data from several runs. If a named
#' vector is given, the cell barcodes will be prefixed with the provided names.
#' This mimics the behavior of the Read10X function found in the Seurat
#' package.
#' @param bam_file bam file from cellranger vdj containing alignment data comparing
#' each contig with the germline reference
#' @return List containing one data.frame for each path provided to vdj_dir
#' @importFrom Rsamtools scanBam
#' @noRd
.load_vdj_indels <- function(vdj_dir, bam_file = "concat_ref.bam") {

  .extract_indels <- function(bam_lst) {

    res <- tibble::tibble(
      cigar     = bam_lst[[1]]$cigar,
      contig_id = bam_lst[[1]]$qname
    )

    res <- dplyr::filter(res, grepl("_contig_[0-9]+$", .data$contig_id))

    # Add indel columns,
    res <- dplyr::mutate(
      res,
      n_insertion = .extract_pat(.data$cigar, "[0-9]+(?=I)"),
      n_deletion  = .extract_pat(.data$cigar, "[0-9]+(?=D)"),
      n_mismatch  = .extract_pat(.data$cigar, "[0-9]+(?=X)"),
    )

    res <- dplyr::select(res, -.data$cigar)

    res
  }

  .extract_pat <- function(string, pat) {
    res <- purrr::map_dbl(string, ~ {
      strg  <- .x
      mtch  <- gregexpr(pat, strg, perl = TRUE)[[1]]
      strts <- as.integer(mtch)
      lens  <- attr(mtch, "match.length", exact = TRUE)

      n <- purrr::map2_int(strts, lens, ~ {
        x <- substr(strg, .x, .x + .y - 1)
        x <- as.integer(x)
        x
      })

      # Multiple indel sites will result in vector with length > 1
      # sum bp from all sites and replace NAs with 0
      n <- tidyr::replace_na(sum(n), 0)

      n
    })

    res
  }

  # Check for file and return path
  # if bam is missing for any sample, return NULL
  # do not want extra NAs in the V(D)J data columns
  res <- purrr::map_chr(vdj_dir, .get_vdj_path, file = bam_file, warn = TRUE)

  if (any(is.na(res))) {
    return(NULL)
  }

  # Load data
  res <- purrr::map(res, Rsamtools::scanBam)

  # Extract indel data
  res <- purrr::map(res, .extract_indels)

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

  # Check vdj_dir for file
  # try adding "outs" to path if file not found
  path <- case_when(
    file.exists(file.path(vdj_dir, file))         ~ file.path(vdj_dir, file),
    file.exists(file.path(vdj_dir, "outs", file)) ~ file.path(vdj_dir, "outs", file)
  )

  if (is.na(path)) {
    fun <- stop

    if (warn) {
      fun <- warning
    }

    fun(file, " not found in ", vdj_dir, ".")
  }

  path
}

#' Check for separator in data.frame
#'
#' @param df_in data.frame
#' @param sep_cols Names of columns to check for sep
#' @param sep Separator to use for storing V(D)J data
#' @return Separator with white space stripped
#' @noRd
.check_sep <- function(df_in, sep_cols, sep) {
  if (is.null(sep_cols)) {
    sep_cols <- colnames(df_in)
  }

  if (is.null(sep)) {
    return(sep)
  }

  if (!is.character(sep)) {
    stop("'sep' must be a character.")
  }

  sep <- gsub("[[:space:]]", "", sep)

  has_sep <- grepl(sep, df_in[, sep_cols, drop = FALSE], fixed = TRUE)

  if (any(has_sep)) {
    stop("The string '", sep, "' is already present in the input data, select a different value for 'sep'.")
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
  .count_chain_class <- function(chains, .df = df_in, .clmn = chain_col) {
    dat <- .df[[.clmn]]
    res <- dat %in% chains
    res <- length(dat[res])

    res
  }

  chains <- list(
    "TCR" = c("TRA", "TRB", "TRD", "TRG"),
    "BCR" = c("IGH", "IGK", "IGL")
  )

  # Count occurrences of BCR/TCR chains
  n_chains <- purrr::map_int(chains, .count_chain_class)

  if (all(n_chains == 0)) {
    chains <- unlist(chains, use.names = FALSE)
    chains <- paste0(chains, collapse = ", ")

    stop("None of the expected chains (", chains, ") were found in the '", chain_col, "' column, unable to determine whether TCR or BCR data were provided.")
  }

  # Calculate fraction of BCR/TCR chains
  # set type if >75% match
  res <- n_chains / sum(n_chains)
  res <- names(res[res > 0.5])

  if (purrr::is_empty(res)) {
    res <- "Multi"

    warning("Equal number of BCR (", n_chains[["BCR"]], ") and TCR (", n_chains[["TCR"]], ") chains detected, unable to determine data type.")
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
.check_overlap <- function(input, meta, nm, pct_min = 25) {

  if (is.null(input)) {
    return(meta)
  }

  obj_meta  <- .get_meta(input)
  obj_cells <- obj_meta[[CELL_COL]]
  met_cells <- unique(meta$barcode)

  n_overlap   <- length(obj_cells[obj_cells %in% met_cells])
  pct_overlap <- round(n_overlap / length(met_cells), 2) * 100

  if (nm != "") {
    nm <- paste0(nm, " ")
  }

  if (identical(n_overlap, 0L)) {
    stop(nm, "cell barcodes do not match those in the object, are you using the correct cell barcode prefixes?")
  }

  if (pct_overlap < pct_min) {
    warning("Only ", pct_overlap, "% (", n_overlap, ") of ", nm, "cell barcodes overlap with the provided object")
  }

  meta
}

#' Identify clonotypes with paired chains
#'
#' @param df_in data.frame containing V(D)J data formatted so each row
#' represents a single contig
#' @return Input data.frame with paired column added
#' @noRd
.identify_paired <- function(df_in) {

  res <- dplyr::group_by(df_in, .data$barcode)

  res <- dplyr::mutate(
    res,
    paired = (all(c("TRA", "TRB") %in% .data$chains)) | ("IGH" %in% .data$chains & any(c("IGL", "IGK") %in% .data$chains))
  )

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

  iso_df <- dplyr::distinct(iso_df, barcode, c_gene)
  iso_df <- stats::na.omit(iso_df)

  dups <- iso_df$barcode
  dups <- dups[duplicated(dups)]

  # Add isotypes to meta.data
  iso_df <- mutate(
    iso_df,
    isotype = ifelse(barcode %in% dups, "Multi", !!sym(iso_col))
  )

  isos <- purrr::set_names(
    iso_df$isotype,
    iso_df$barcode
  )

  res <- mutate(
    df_in,
    isotype = unname(isos[barcode]),
    isotype = tidyr::replace_na(isotype, "None")
  )

  res
}

#' Define clonotypes based on V(D)J data
#'
#' This will assign new clonotype IDs based on the combination of values
#' present in the provided columns
#'
#' @param input Single cell object or data.frame containing V(D)J data. If a
#' data.frame is provided, the cell barcodes should be stored as row names.
#' @param vdj_cols meta.data columns containing V(D)J data to use for defining
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
#'   vdj_cols = "cdr3_nt"
#' )
#'
#' head(res@meta.data, 1)
#'
#' # Define clonotypes based on the combination of the CDR3 nucleotide sequence
#' # and the V and J genes
#' res <- define_clonotypes(
#'   vdj_sce,
#'   vdj_cols = c("cdr3_nt", "v_gene", "j_gene")
#' )
#'
#' head(res@colData, 1)
#'
#' # Modify the name of the column used to store clonotype IDs
#' res <- define_clonotypes(
#'   vdj_so,
#'   vdj_cols = "cdr3_nt",
#'   clonotype_col = "NEW_clonotype_id"
#' )
#'
#' head(res@meta.data, 1)
#'
#' # When defining clonotypes only use chains that are productive
#' res <- define_clonotypes(
#'   vdj_sce,
#'   vdj_cols = "cdr3_nt",
#'   filter_chains = "productive"
#' )
#'
#' head(res@colData, 1)
#'
#' @export
define_clonotypes <- function(input, vdj_cols, clonotype_col = "clonotype_id",
                              filter_chains = c("productive", "full_length"), sep = ";") {

  # Get meta.data
  meta <- .get_meta(input)

  if (!all(vdj_cols %in% colnames(meta))) {
    stop("Not all vdj_cols (", paste0(vdj_cols, collapse = ", "), ") are present in meta.data.")
  }

  # Only use values in vdj_cols that are TRUE for all filter_chains columns
  # first identify contigs TRUE for all filter_chains columns
  # subset each vdj_cols column based on .clone_idx
  if (!is.null(filter_chains)) {

    clmns <- syms(filter_chains)

    meta <- mutate_vdj(
      input,

      .clone_idx = list(
        purrr::reduce(list(!!!clmns), ~ .x & .y)
      ),

      dplyr::across(
        dplyr::all_of(vdj_cols),
        ~ paste0(.x[.data$.clone_idx], collapse = ""),
        .names = ".clone_{.col}"
      ),

      clonotype_col = NULL,

      vdj_cols = c(vdj_cols, filter_chains)
    )

    vdj_cols <- paste0(".clone_", vdj_cols)

    meta <- .get_meta(meta)
  }

  # Add new clonotype IDs
  meta <- dplyr::mutate(
    meta,
    .new_clone = paste(!!!syms(vdj_cols), sep = ""),
    .new_id    = rank(.data$.new_clone, ties.method = "min"),

    !!sym(clonotype_col) := ifelse(
      .data$.new_clone == "",
      "None",
      paste0("clonotype", .data$.new_id)
    )
  )

  # Remove temporary columns
  meta <- dplyr::select(
    meta,
    -dplyr::matches("^.new_(clone|id)"),
    -dplyr::starts_with(".clone_")
  )

  # Format results
  res <- .add_meta(input, meta)

  res
}


