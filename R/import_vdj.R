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

  # SHOULD CHECK BARCODE OVERLAP BEFORE BINDING ROWS

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
  contigs <- dplyr::rowwise(contigs)

  contigs <- dplyr::mutate(
    contigs,
    across(
      all_of(cdr3_cols),
      ~ length(strsplit(.x, "")[[1]]),
      .names = "{.col}_length"
    )
  )

  contigs <- dplyr::ungroup(contigs)

  sep_cols <- c(sep_cols, paste0(cdr3_cols, "_length"))

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
  res <- tibble::column_to_rownames(meta, "barcode")
  res <- dplyr::rename_with(res, ~ paste0(prefix, .x))

  if (!is.null(input)) {
    res <- .merge_meta(input, res)
  }

  res
}
