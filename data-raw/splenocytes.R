library(djvdj)
library(Seurat)
library(here)
library(tidyverse)
library(Matrix)
library(Rsamtools)
library(clustifyr)
library(clustifyrdata)   # USE EXPERIMENTHUB INSTEAD


# Functions ----

# Demultiplex samples
demux_so <- function(dat_dir) {

  # Create Seurat object
  dat_dir <- file.path(
    dat_dir,
    "JH179_GEX-JH181_ADT_JH181_HTO/outs/filtered_feature_bc_matrix"
  )

  mat <- Read10X(dat_dir)

  so <- mat$`Gene Expression` %>%
    CreateSeuratObject()

  # Add HTO data and demultiplex
  hto_mat <- mat$`Antibody Capture`
  hto_mat <- hto_mat[rownames(hto_mat) != "_HEL", ]

  so[["HTO"]] <- CreateAssayObject(hto_mat)

  so <- so %>%
    NormalizeData(assay = "HTO", normalization.method = "CLR") %>%
    HTODemux(assay = "HTO") %>%
    subset(hash.ID != "Doublet" & hash.ID != "Negative")

  # Preprocess object
  so <- so %>%
    preprocess_so()

  so
}

# Preprocess object
preprocess_so <- function(so_in) {
  res <- so_in %>%
    NormalizeData() %>%
    FindVariableFeatures(nfeatures = 1000) %>%
    ScaleData() %>%
    RunPCA() %>%
    RunUMAP(dims = 1:30) %>%
    FindNeighbors(dims = 1:30) %>%
    FindClusters(resolution = 1)

  # Annotate B/T cells
  res <- res %>%
    clustify(
      cluster_col = "seurat_clusters",
      ref_mat = ref_immgen
    )

  res <- res %>%
    mutate_meta(
      mutate,
      cell_type = str_remove(type, " \\(.+$"),
      cell_type = recode(cell_type, Tgd = "T cells"),
      cell_type = ifelse(
        cell_type %in% c("B cells", "T cells"),
        cell_type,
        "other"
      )
    )

  # Only include variable genes in object
  res <- res %>%
    subset(features = VariableFeatures(res))

  res
}

# Downsample object
downsample_so <- function(so_in, n, b_frac = 0.8) {

  if (is.null(n) || is.null(b_frac)) return(so_in)

  # Sample B and T cells
  n_b <- round(n * b_frac, 1)
  n_t <- n - n_b

  set.seed(42)
  b_bcs <- so %>%
    subset(cell_type == "B cells") %>%
    colnames() %>%
    sample(n_b)

  set.seed(42)
  t_bcs <- so %>%
    subset(cell_type == "T cells") %>%
    colnames() %>%
    sample(n_t)

  bcs <- c(b_bcs, t_bcs)

  # Subset object
  res <- so_in[, colnames(so_in) %in% bcs]

  res
}

# Write matrices
write_mtx <- function(mtx, bcs, old_rna = NULL, old_vdj = NULL, new_rna = NULL,
                      new_vdj = NULL) {

  # Write RNA data
  if (!is.null(old_rna) && !is.null(new_rna)) {
    dir.create(new_rna, mode = "0750", recursive = TRUE)

    feats <- here(old_rna, "features.tsv.gz") %>%
      read_tsv(col_names = FALSE) %>%
      filter(
        X3 == "Gene Expression",
        X2 %in% rownames(mtx)
      ) %>%
      group_by(X2) %>%
      dplyr::slice(1) %>%
      ungroup()

    feats %>%
      write_tsv(here(new_rna, "features.tsv.gz"), col_names = FALSE)

    bcs %>%
      as.data.frame() %>%
      write_tsv(here(new_rna, "barcodes.tsv.gz"), col_names = FALSE)

    mtx[feats$X2, bcs] %>%
      writeMM(here(new_rna, "matrix.mtx.gz"))
  }

  # Write VDJ data
  if (!is.null(old_vdj) && !is.null(new_vdj)) {
    dir.create(new_vdj, mode = "0750", recursive = TRUE)

    # Format and write BCR contig data
    vdj_dat <- here(old_vdj, "filtered_contig_annotations.csv") %>%
      read_csv() %>%
      filter(barcode %in% bcs)

    vdj_dat %>%
      write_csv(here(new_vdj, "filtered_contig_annotations.csv.gz"))

    here(old_vdj, "airr_rearrangement.tsv") %>%
      read_tsv() %>%
      filter(cell_id %in% bcs) %>%
      write_tsv(here(new_vdj, "airr_rearrangement.tsv.gz"))

    # Set bam filtering rules
    filt_fn <- function(x) {

      clns <- x$rname %>%
        str_remove("_concat_ref_[0-9]+$")

      cells <- x$qname %>%
        str_remove("_contig_[0-9]+$")

      clns %in% clones & (cells %in% bcs | grepl("^clonotype", cells))
    }

    filt_rules <- FilterRules(list(rname = filt_fn))

    clones <- vdj_dat$raw_clonotype_id

    # Filter bam file
    bam <- here(new_vdj, "concat_ref.bam")

    indexBam(here(old_vdj, "concat_ref.bam"))

    here(old_vdj, "concat_ref.bam") %>%
      filterBam(bam, filter = filt_rules)

    # Remove bam indices
    here(old_vdj, "concat_ref.bam.bai") %>%
      file.remove()

    str_c(bam, ".bai") %>%
      file.remove()
  }
}

write_mtxs <- function(so_in, dat_dir, write_dir) {

  # Split RNA matrices
  bcs <- so_in@meta.data %>%
    mutate(
      hash.ID = recode(
        as.character(hash.ID),
        "CD45-and-MHC-class-I-C0301" = "BL6",
        "CD45-and-MHC-class-I-C0302" = "MD4"
      )
    ) %>%
    split(.$hash.ID) %>%
    map(rownames)

  # Write new matrices
  old_rna_dir <- file.path(
    dat_dir,
    "JH179_GEX-JH181_ADT_JH181_HTO/outs/filtered_feature_bc_matrix"
  )

  old_vdj_dirs <- c(
    BCR = file.path(dat_dir, "JH180_BCR/outs/"),
    TCR = file.path(dat_dir, "JH180_TCR/outs/")
  )

  rna_dirs <- c(
    BL6 = here(write_dir, "BL6_GEX/filtered_feature_bc_matrix"),
    MD4 = here(write_dir, "MD4_GEX/filtered_feature_bc_matrix")
  )

  bcr_dirs <- c(
    BL6 = here(write_dir, "BL6_BCR"),
    MD4 = here(write_dir, "MD4_BCR")
  )

  tcr_dirs <- c(
    BL6 = here(write_dir, "BL6_TCR"),
    MD4 = here(write_dir, "MD4_TCR")
  )

  # Write RNA and BCR data
  names(rna_dirs) %>%
    walk(~ {
      so_in@assays$RNA@counts %>%
        write_mtx(
          bcs     = bcs[[.x]],
          old_rna = old_rna_dir,
          old_vdj = old_vdj_dirs[["BCR"]],
          new_rna = rna_dirs[[.x]],
          new_vdj = bcr_dirs[[.x]]
        )

      # Write TCR data
      so_in@assays$RNA@counts %>%
        write_mtx(
          bcs     = bcs[[.x]],
          old_vdj = old_vdj_dirs[["TCR"]],
          new_vdj = tcr_dirs[[.x]]
        )
    })
}

# Write Seurat object
write_so <- function(dat_dir, write_dir, file, n_grps = 5) {
  rna_dirs <- c(
    BL6 = here(dat_dir, "BL6_GEX/filtered_feature_bc_matrix"),
    MD4 = here(dat_dir, "MD4_GEX/filtered_feature_bc_matrix")
  )

  # Create Seurat objects
  so <- rna_dirs %>%
    Read10X() %>%
    CreateSeuratObject(min.cells = 1) %>%
    preprocess_so()

  # Add mock groups
  so <- so %>%
    mutate_meta(~ {
      .x %>%
        group_by(orig.ident) %>%
        mutate(sample = str_c(orig.ident, "-", ntile(orig.ident, n_grps))) %>%
        ungroup()
    })

  splen_meta <- so@meta.data

  # Save Seurat meta.data
  usethis::use_data(splen_meta, compress = "xz")
}


# Write files ----

# Demultiplex samples
dat_dir <- "~/Projects/Smith_AVIDseq/results/2020-01-18/"

so <- demux_so(dat_dir)

# Write files
params <- list(
  file   = "splen",
  n      = 500,
  frac   = 0.8,
  n_grps = 3
)

params %>%
  pwalk(~ {
    params <- list(...)

    mtx_dir <- here("inst/extdata", params$file)

    so <- so %>%
      downsample_so(n = params$n, b_frac = params$frac)

    write_mtxs(so, dat_dir, mtx_dir)

    write_so(
      dat_dir   = mtx_dir,
      write_dir = "data",
      file      = params$file,
      n_grps    = params$n_grps
    )
  })
