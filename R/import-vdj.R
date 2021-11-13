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
#' @param filter_contigs Only include chains with at least one productive and
#' full length contig. It is highly recommended that filter_contigs is set to
#' TRUE when using the define_clonotypes argument.
#' @param filter_paired Only include clonotypes with paired chains. For TCR
#' data each clonotype must have at least one TRA and TRB chain, for BCR data
#' each clonotype must have at least one IGH chain and at least one IGK or IGL
#' chain.
#' @param define_clonotypes Define clonotype IDs based on V(D)J data. This is
#' useful if the V(D)J datasets being loaded do not have consistent clonotype
#' IDs, e.g., clonotype1 is not the same across samples.
#' 'cdr3aa' will use the CDR3 amino acid sequence, 'cdr3nt' will use the CDR3
#' nucleotide sequence, and 'cdr3_gene' will use the combination of the CDR3
#' nucleotide sequence and the V(D)J genes. Set to NULL (default) to use the
#' clonotype IDs already present in the input data.
#' @param sep Separator to use for storing per cell V(D)J data
#' @return Single cell object or data.frame with added V(D)J data
#' @export
import_vdj <- function(input = NULL, vdj_dir, prefix = "", cell_prefix = NULL, filter_contigs = TRUE,
                       filter_paired = FALSE, define_clonotypes = NULL, sep = ";") {

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
  vdj_class <- unique(vdj_class)

  if (length(vdj_class) > 1) {
    vdj_class <- paste0(vdj_class, collapse = ", ")

    stop("Multiple data types detected (", vdj_class, "), provided data must be either TCR or BCR. To add both TCR and BCR data to the same object, run import_vdj separately for each and use the 'prefix' argument to add different column names.")
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
    warning(n_remove, " contigs do not have an assigned clonotype_id, these contigs will be removed.")

    contigs <- dplyr::filter(contigs, !is.na(.data$raw_clonotype_id))
  }

  # Select V(D)J columns to keep
  contigs <- dplyr::select(contigs, all_of(vdj_cols))

  # Check if sep is already in sep_cols
  sep <- gsub("[[:space:]]", "", sep)

  if (any(grepl(sep, contigs[, sep_cols]))) {
    stop("The string '", sep, "' is already present in the input data, select a different value for 'sep'.")
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

  # Determine which clonotypes are paired
  contigs <- .identify_paired(contigs)

  if (filter_paired) {
    contigs <- dplyr::filter(contigs, .data$paired)
  }

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

  meta <- dplyr::relocate(meta, .data$paired, .after = .data$full_length)
  meta <- dplyr::relocate(meta, .data$isotype, .after = .data$chains)

  # Check for duplicated cell barcodes
  if (any(duplicated(meta$barcode))) {
    stop("Malformed inport data, multiple clonotype_ids are associated with the same cell barcode.")
  }

  # Allow user to redefine clonotypes
  res <- tibble::column_to_rownames(meta, "barcode")

  if (!is.null(define_clonotypes)) {

    if (!filter_contigs) {
      warning("It is highly recommended that filter_contigs be set to TRUE when using the define_clonotypes argument.")
    }

    clone_cols <- list(
      cdr3aa    = "cdr3",
      cdr3nt    = "cdr3_nt",
      cdr3_gene = c("cdr3_nt", "v_gene", "d_gene", "j_gene")
    )

    if (!define_clonotypes %in% names(clone_cols)) {
      stop("define_clonotypes must be one of ", paste0(names(clone_cols), collapse = ", "), ".")
    }

    clone_cols <- clone_cols[[define_clonotypes]]

    res <- define_clonotypes(res, vdj_cols = clone_cols)
  }

  # Filter to only include cells with valid clonotype_id
  # cells with missing clonotype have a clonotype_id of 'None'
  res <- dplyr::filter(res, grepl("^clonotype", .data$clonotype_id))

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
#' contig annotation
#' @return List containing one data.frame for each path provided to vdj_dir
.load_vdj_data <- function(vdj_dir, contig_file = "filtered_contig_annotations.csv") {

  # Check vdj_dir for contig_file
  # try adding "outs" to path if contig_file is not found
  res <- purrr::map_chr(vdj_dir, ~ {
    path <- case_when(
      file.exists(file.path(.x, contig_file))         ~ file.path(.x, contig_file),
      file.exists(file.path(.x, "outs", contig_file)) ~ file.path(.x, "outs", contig_file)
    )

    if (is.na(path)) {
      stop(contig_file, " not found in ", .x, ".")
    }

    path
  })

  # Load data
  res <- purrr::map(res, readr::read_csv, col_types = readr::cols())

  # Add cell prefixes and replace 'None' in productive with FALSE
  res <- purrr::imap(res, ~ {
    .x <- dplyr::rowwise(.x)
    .x <- dplyr::mutate(
      .x,
      barcode      = paste0(.y, .data$barcode),
      clonotype_id = paste0(.y, raw_clonotype_id),
      dplyr::across(c(full_length, productive), ~ {
        ifelse(
          is.character(.x),
          as.logical(gsub("None", "False", .x)),
          .x
        )
      })
    )

    .x <- dplyr::ungroup(.x)

    dplyr::rename(.x, chains = .data$chain)
  })

  res
}

#' Determine whether TCR or BCR data were provided
#'
#' @param df_in data.frame containing V(D)J data formatted so that each row
#' represents a single contig
#' @param chain_col Column in input data containing chain identity
#' @return Character string indicating whether TCR or BCR data were provided
.classify_vdj <- function(df_in, chain_col = "chains") {
  .count_chain_class <- function(chains) {
    dat <- df_in[[chain_col]]
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

  if (length(res) == 0) {
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

#' Identify clonotypes with paired chains
#'
#' @param df_in data.frame containing V(D)J data formatted so each row
#' represents a single contig
#' @return Input data.frame with paired column added
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
#' @param iso_pat Regular expression to use for extracting isotypes from
#' iso_col
#' @param iso_col Column containing data to use for extracting isotypes
#' @return Input data.frame with isotype column added
.extract_isotypes <- function(df_in, iso_pat = "^IGH[ADEGM]", iso_col = "c_gene") {

  res <- dplyr::group_by(df_in, .data$barcode)

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

#' Define clonotypes based on V(D)J data
#'
#' @param input Single cell object or data.frame containing V(D)J data. If a
#' data.frame is provided, the cell barcodes should be stored as row names.
#' @param vdj_cols meta.data columns containing V(D)J data to use for defining
#' clonotypes
#' @param clonotype_col Name of column to use for storing clonotype IDs
#' @param filter_contigs To only use productive and full length chains when
#' defining clonotypes, specify column names. If set to NULL chains are not
#' filtered.
#' @param sep Separator used for storing per cell V(D)J data
#' @return Single cell object or data.frame with added clonotype IDs
#' @export
define_clonotypes <- function(input, vdj_cols, clonotype_col = "clonotype_id",
                              filter_contigs = c("productive", "full_length"), sep = ";") {

  # Get meta.data
  meta <- .get_meta(input)

  if (!all(vdj_cols %in% colnames(meta))) {
    stop("Not all vdj_cols (", paste0(vdj_cols, collapse = ", "), ") are present in meta.data.")
  }

  # Only use values in vdj_cols that are TRUE for all filter_contigs columns
  # first identify contigs TRUE for all filter_contigs columns
  # subset each vdj_cols column based on .clone_idx
  if (!is.null(filter_contigs)) {
    meta <- mutate_vdj(
      input,
      .clone_idx = list(purrr::reduce(list(!!!syms(filter_contigs)), ~ .x & .y)),
      dplyr::across(
        dplyr::all_of(vdj_cols),
        ~ paste0(.x[.data$.clone_idx], collapse = ""),
        .names = ".clone_{.col}"
      ),
      clonotype_col = NULL,
      vdj_cols      = c(vdj_cols, "productive", "full_length")
    )

    vdj_cols <- paste0(".clone_", vdj_cols)

    meta <- .get_meta(meta)
  }

  # Add new clonotype IDs
  meta <- meta %>%
    mutate(
      .new_clone            = paste(!!!syms(vdj_cols), sep = ""),
      .new_id               = rank(.new_clone, ties.method = "min"),
      !!sym(clonotype_col) := ifelse(.new_clone == "", "None", paste0("clonotype", .data$.new_id))
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


