
library(djvdj)
library(Seurat)
library(here)
library(tidyverse)
library(Matrix)
library(Rsamtools)



# ---- WRITE RAW MATRICES

# Create Seurat object and demultiplet samples
# Save a separate matrix for each sample

dat_dir <- c("~/Projects/Smith_AVIDseq/results/2020-01-18/")

# Function to write matrices
write_mtx <- function(mtx, bcs, old_rna = NULL, old_vdj = NULL, new_rna = NULL,
                      new_vdj = NULL) {

  # Write RNA data
  if (!is.null(old_rna) && !is.null(new_rna)) {
    dir.create(new_rna, mode = "0750", recursive = TRUE)

    here(old_rna, "features.tsv.gz") %>%
      read_tsv(col_names = FALSE) %>%
      filter(X3 == "Gene Expression") %>%
      write_tsv(here(new_rna, "features.tsv.gz"), col_names = FALSE)

    bcs %>%
      as.data.frame() %>%
      write_tsv(here(new_rna, "barcodes.tsv.gz"), col_names = FALSE)

    mtx[, bcs] %>%
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
      write_csv(here(new_vdj, "filtered_contig_annotations.csv"))

    here(old_vdj, "airr_rearrangement.tsv") %>%
      read_tsv() %>%
      filter(cell_id %in% bcs) %>%
      write_tsv(here(new_vdj, "airr_rearrangement.tsv"))

    # Set bam filtering rules
    filt_fn <- function(x) {
      clns <- x$rname %>%
        str_remove("_concat_ref_[0-9]+$")

      clns %in% clones
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

# Create Seurat object
old_rna_dir <- file.path(dat_dir, "JH179_GEX-JH181_ADT_JH181_HTO/outs/filtered_feature_bc_matrix")
mat <- Read10X(old_rna_dir)

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

# Split RNA matrices
write_dir <- here("data/splenocytes")

bcs <- so@meta.data %>%
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
    so@assays$RNA@counts %>%
      write_mtx(
        bcs     = bcs[[.x]],
        old_rna = old_rna_dir,
        old_vdj = old_vdj_dirs[["BCR"]],
        new_rna = rna_dirs[[.x]],
        new_vdj = bcr_dirs[[.x]]
      )

    # Write TCR data
    so@assays$RNA@counts %>%
      write_mtx(
        bcs     = bcs[[.x]],
        old_vdj = old_vdj_dirs[["TCR"]],
        new_vdj = tcr_dirs[[.x]]
      )
  })



# ---- SAVE SEURAT OBJECTS

# Create Seurat objects
splen_so <- rna_dirs %>%
  Read10X() %>%
  CreateSeuratObject() %>%
  NormalizeData() %>%
  FindVariableFeatures(nfeatures = 1000) %>%
  ScaleData() %>%
  RunPCA() %>%
  RunUMAP(dims = 1:30) %>%
  FindNeighbors(dims = 1:30) %>%
  FindClusters(resolution = 1)

# Annotate B/T cells
splen_so <- splen_so %>%
  mutate_meta(
    mutate,
    cell_type = case_when(
      seurat_clusters %in% c(3, 15, 12, 0, 19) ~ "B cells",
      seurat_clusters %in% c(7, 20, 14, 4, 21, 17, 1, 2, 11, 10, 8, 5) ~ "T cells",
      TRUE ~ "other"
    )
  )

# Add mock groups
splen_so <- splen_so %>%
  mutate_meta(~ {
    .x %>%
      group_by(orig.ident) %>%
      mutate(sample = str_c(orig.ident, "-", ntile(orig.ident, 5))) %>%
      ungroup()
  })

# Save object
save(splen_so, file = here("data/splenocytes/splen_so.rda"))
