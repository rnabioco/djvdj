#' Import V(D)J data
#'
#' @param input Object containing single cell data, if set to NULL a data.frame
#' containing V(D)J results will be returned
#' @param vdj_dir Directory containing the output from cellranger vdj. A vector
#' or named vector can be given to load data from several runs. If a named
#' vector is given, the cell barcodes will be prefixed with the provided names.
#' This mimics the behavior of the Read10X function found in the Seurat
#' package. Cell barcode prefixes can also be provided using the cell_prefix
#' argument.
#' @param prefix Prefix to add to new columns
#' @param cell_prefix Prefix to add to cell barcodes, this is helpful when
#' loading data from multiple runs into a single object. If set to NULL, cell
#' barcode prefixes will be automatically generated in a similar way as the
#' Read10X function found in the Seurat package.
#' @param filter_contigs Only include chains with at least one productive
#' contig
#' @param sep Separator to use for storing per cell V(D)J data
#' @return Single cell object or data.frame with added V(D)J data
#' @export
import_vdj <- function(input = NULL, vdj_dir, prefix = "", cell_prefix = NULL, filter_contigs = TRUE,
                       sep = ";") {

  # VDJ columns
  count_cols <- c("reads", "umis")
  cdr3_cols  <- c("cdr3", "cdr3_nt")

  sep_cols <- c(
    "v_gene",   "d_gene",
    "j_gene",   "c_gene",
    "chains",   cdr3_cols,
    count_cols, "productive",
    "full_length"
  )

  vdj_cols <- c(
    "barcode", "clonotype_id",
    sep_cols
  )

  # Format/check cell prefixes
  vdj_dir <- .format_cell_prefixes(vdj_dir, cell_prefix)

  # Load V(D)J data and add cell prefixes
  contigs <- .load_vdj_data(vdj_dir)

  # Classify input data as TCR or BCR
  vdj_class <- purrr::map_chr(contigs, .classify_vdj)

  if (length(unique(vdj_class)) > 1) {
    stop("Provided data must be either TCR or BCR. To add both TCR and BCR data to the same object, run import_vdj separately for each and use the 'prefix' argument to add different column names.")
  }

  # Check cell barcode overlap
  # give warning for low overlap
  contigs <- purrr::imap_dfr(contigs, ~ .check_overlap(input, .x, .y))

  # Filter for productive contigs
  if (filter_contigs) {
    contigs <- dplyr::filter(contigs, .data$productive, .data$full_length)
  }

  # Remove contigs that do not have an assigned clonotype_id
  n_remove <- contigs$raw_clonotype_id
  n_remove <- n_remove[is.na(n_remove)]
  n_remove <- length(n_remove)

  if (n_remove > 0) {
    warning(n_remove, " contigs do not have an assigned clonotype_id, these contigs will be removed")

    contigs <- dplyr::filter(contigs, !is.na(.data$raw_clonotype_id))
  }

  # Select V(D)J columns to keep
  contigs <- dplyr::select(contigs, all_of(vdj_cols))

  # Check if sep is already in sep_cols
  sep <- gsub("[[:space:]]", "", sep)

  if (any(grepl(sep, contigs[, sep_cols]))) {
    stop("The string '", sep, "' is already present in the V(D)J data, select a different value for sep")
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

  # Calculate CDR3 length
  contigs <- dplyr::mutate(
    contigs,
    across(all_of(cdr3_cols), nchar, .names = "{.col}_length")
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
  if (vdj_class == "BCR") {
    contigs <- .extract_isotypes(contigs)
    contigs <- dplyr::group_by(contigs, .data$isotype)
  }

  # Collapse chains into a single row for each cell
  # include isotype and clonotype_id as groups so they are included in the
  # summarized results
  contigs <- dplyr::group_by(
    contigs,
    .data$barcode, .data$clonotype_id,
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

  # Check for duplicated cell barcodes
  if (any(duplicated(meta$barcode))) {
    stop("Malformed inport data, multiple clonotype_ids are associated with the same cell barcode")
  }

  # Format meta.data
  # currently the number of cell barcodes that overlap with the object is
  # checked with .merge_meta, ideally each sample would be checked separately
  res <- tibble::column_to_rownames(meta, "barcode")
  res <- dplyr::rename_with(res, ~ paste0(prefix, .x))

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
#' @param sep Separator to use when appending prefixes to cell barcodes. Set to
#' NULL to not add a separator.
#' @return Paths provided to vdj_dir with cell prefixes added as names.
.format_cell_prefixes <- function(vdj_dir, cell_prefix, sep = "_") {

  res <- vdj_dir

  # If vdj_dir is not named, check cell_prefix
  if (is.null(names(res))) {

    # If no cell prefixes are provided, auto generate
    if (is.null(cell_prefix)) {
      cell_prefix    <- paste0(seq_along(res), "_")
      cell_prefix[1] <- ""
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
#' contig annotation.
#' @return List containing one data.frame for each path provided to vdj_dir.
.load_vdj_data <- function(vdj_dir, contig_file = "filtered_contig_annotations.csv") {

  # Check vdj_dir for contig_file
  # try adding "outs" to path if contig_file is not found
  res <- purrr::map_chr(vdj_dir, ~ {
    path <- case_when(
      file.exists(file.path(.x, contig_file))         ~ file.path(.x, contig_file),
      file.exists(file.path(.x, "outs", contig_file)) ~ file.path(.x, "outs", contig_file)
    )

    if (is.na(path)) {
      stop(res, " not found in ", vdj_dir, ".")
    }

    path
  })

  # Load data
  res <- purrr::map(res, readr::read_csv, col_types = readr::cols())

  # Add cell prefixes
  res <- purrr::imap(res, ~ {
    .x <- dplyr::mutate(
      .x,
      barcode      = paste0(.y, .data$barcode),
      clonotype_id = paste0(.y, raw_clonotype_id)
    )

    dplyr::rename(.x, chains = .data$chain)
  })

  res
}

#' Determine whether TCR or BCR data were provided
#'
#' @param input data.frame containing V(D)J data formatted so that each row
#' represents a single contig.
#' @param chain_col Column in input data containing chain identity.
#' @return Character string indicating whether TCR or BCR data were provided.
.classify_vdj <- function(input, chain_col = "chains") {
  get_chain_class <- function(chains) {
    any(chains %in% input[[chain_col]])
  }

  chains <- list(
    "TCR" = c("TRA", "TRB", "TRD", "TRG"),
    "BCR" = c("IGH", "IGK", "IGL")
  )

  res <- purrr::map_lgl(chains, get_chain_class)
  res <- names(res[res])

  if (length(res) == 0) {
    chains <- unlist(chains, use.names = FALSE)
    chains <- paste0(chains, collapse = ", ")

    stop("None of the expected chains (", chains, ") were found in the '", chain_col, "' column, unable to determine whether TCR or BCR data were provided.")

  } else if (length(res) > 1) {
    stop("Malformed input data, both TCR and BCR chains are present.")
  }

  res
}

#' Check cell barcode overlap with object
#'
#' @param input Object containing single cell data.
#' @param meta meta.data to check against object.
#' @param nm Sample name to use for messages.
#' @param pct_min Warn user if the percent overlap is less than pct_min.
#' @return input data.
.check_overlap <- function(input, meta, nm, pct_min = 25) {

  if (is.null(input)) {
    return(meta)
  }

  obj_meta  <- .get_meta(input)
  obj_cells <- obj_meta$.cell_id
  met_cells <- unique(meta$barcode)

  n_overlap   <- length(obj_cells[obj_cells %in% met_cells])
  pct_overlap <- round(n_overlap / length(met_cells), 2) * 100

  if (nm != "") {
    nm <- paste0(nm, " ")
  }

  if (n_overlap == 0) {
    stop(nm, "cell barcodes do not match those in the object, are you using the correct cell barcode prefixes?")
  }

  if (pct_overlap < pct_min) {
    warning("Only ", pct_overlap, "% (", n_overlap, ") of ", nm, "cell barcodes overlap with the provided object")
  }

  meta
}

#' Add isotypes to V(D)J data
#'
#' @param input data.frame containing V(D)J data formatted so that each row
#' represents a single contig.
#' @param iso_pat Regular expression to use for extracting isotypes from
#' iso_col.
#' @param iso_col Column containing data to use for extracting isotypes.
#' @return input data.frame with isotype column added.
.extract_isotypes <- function(input, iso_pat = "^IGH[ADEGM]", iso_col = "c_gene") {

  res <- dplyr::group_by(input, .data$barcode)

  res <- dplyr::mutate(
    res,
    isotype = list(
      substr(
        !!sym(iso_col), 1,
        attr(regexpr(iso_pat, !!sym(iso_col)), "match.length", exact = TRUE)
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

  res <- dplyr::ungroup(res)

  res
}


