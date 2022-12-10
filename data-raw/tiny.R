library(djvdj)
library(tidyverse)
library(usethis)
library(Seurat)
library(SingleCellExperiment)
library(Rsamtools)

# Parameters
dat_dir <- "~/Projects/Smith_AVIDseq"
out_dir <- "inst/extdata/"

mat_path <- c(
  str_c(dat_dir, "/results/2020-01-18/JH179_GEX-JH181_ADT_JH181_HTO/outs/filtered_feature_bc_matrix"),
  str_c(dat_dir, "/results/2020-01-18/JH179_GEX-JH181_ADT_JH181_HTO/outs/filtered_feature_bc_matrix")
)

bcr_path <- c(
  str_c(dat_dir, "/results/2020-01-18/JH180_BCR/outs/"),
  str_c(dat_dir, "/results/2020-01-18/JH180_BCR/outs/")
)

tcr_path <- str_c(dat_dir, "/results/2020-01-18/JH180_TCR/outs/")





### PROCESS GEX ----

# Create Suerat object
mat <- Read10X(mat_path)

tiny_so <- mat$`Gene Expression` %>%
  CreateSeuratObject() %>%
  mutate_meta(
    mutate,
    orig.ident = ifelse(str_detect(.cell_id, "^2_"), "avid_2", "avid_1")
  )

# Identify cells for tiny objects
# want to include some cells with non-productive chains
bcr_contigs <- bcr_path %>%
  map(str_c, "filtered_contig_annotations.csv") %>%
  map(read_csv)

bad_chain_cells <- bcr_contigs[[1]] %>%
  filter(!productive) %>%
  head(10) %>%
  pull(barcode) %>%
  unique() %>%
  str_c("1_", .)

all_cells <- colnames(tiny_so)

set.seed(100)

tiny_cells <- all_cells %>%
  sample(200) %>%
  c(bad_chain_cells)

tiny_cells <- all_cells[all_cells %in% tiny_cells]

# Subset object
# keep cell order the same
set.seed(100)

tiny_so <- tiny_so %>%
  subset(
    features = sample(rownames(tiny_so), 200),
    cells    = tiny_cells
  )

# Normalize object
tiny_so <- tiny_so %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData()

# Run PCA, cluster, UMAP
tiny_so <- tiny_so %>%
  RunPCA(assay = "RNA") %>%
  FindNeighbors(
    reduction = "pca",
    dims      = 1:40
  ) %>%
  FindClusters(
    resolution = 0.3,
    verbose    = FALSE
  ) %>%
  RunUMAP(
    assay = "RNA",
    dims  = 1:40
  )

tiny_so <- tiny_so %>%
  AddMetaData(FetchData(., c("UMAP_1", "UMAP_2")))





### FORMAT BCR DATA ----

# Format and write BCR contig data
bcr_contigs <- bcr_contigs %>%
  imap(~ filter(.x, str_c(.y, "_", barcode) %in% tiny_cells))

names(bcr_contigs) <- str_c("bcr_", 1:2) %>%
  str_c(out_dir, ., "/outs/filtered_contig_annotations.csv")

clones <- bcr_contigs %>%
  map(~ unique(.x$raw_clonotype_id)) %>%
  purrr::reduce(c)

bcr_contigs %>%
  iwalk(write_csv)

# Format and write BCR bams
bcr_bam <- bcr_path %>%
  str_c("/concat_ref.bam")

names(bcr_bam) <- str_c("bcr_", 1:2) %>%
  str_c(out_dir, ., "/outs/concat_ref.bam")

filt_fn <- function(x) {
  clns <- x$rname %>%
    str_remove("_concat_ref_[0-9]+$")

  clns %in% clones
}

filt_rules <- FilterRules(list(rname = filt_fn))

bcr_bam %>%
  walk(indexBam)

bcr_bam %>%
  iwalk(filterBam, filter = filt_rules)

names(bcr_bam) %>%
  str_c(".bai") %>%
  file.remove()





### FORMAT BAD DATA ----

# Format bad BCR contig data
# this is to test:
# * low overlap warning
# * missing clonotype_id
# * non-overlapping indel data
bad_path <- str_c(out_dir, "bad_bcr_1/outs/")

bad_contigs <- bcr_contigs[[1]] %>%
  head(1) %>%
  bind_rows(bcr_contigs[[2]]) %>%
  mutate(raw_clonotype_id = ifelse(rownames(.) == 1, NA, raw_clonotype_id))

bad_contigs <- set_names(
  list(bad_contigs),
  str_c(bad_path, "filtered_contig_annotations.csv")
)

bad_contigs %>%
  iwalk(write_csv)

# Use TCR concat_ref.bam
bad_bam <- set_names(
  str_c(tcr_path, "concat_ref.bam"),
  str_c(bad_path, "concat_ref.bam")
)

bad_bam %>%
  walk(indexBam)

bad_bam %>%
  iwalk(filterBam, filter = filt_rules)

names(bad_bam) %>%
  str_c(".bai") %>%
  file.remove()





### FORMAT TCR DATA ----

# Format and write TCR contig data
# add some NAs to raw_clonotype_id for testing
tcr_contigs <- tcr_path %>%
  str_c("filtered_contig_annotations.csv") %>%
  read_csv()

tcr_contigs <- tcr_contigs %>%
  filter(str_c("1_", barcode) %in% tiny_cells)

tcr_contigs %>%
  write_csv(str_c(out_dir, "tcr_1/outs/filtered_contig_annotations.csv"))

# Format and write TCR bam
clones <- unique(tcr_contigs$raw_clonotype_id)

tcr_bam <- set_names(
  str_c(tcr_path, "concat_ref.bam"),
  str_c(out_dir, "tcr_1/outs/concat_ref.bam")
)

tcr_bam %>%
  walk(indexBam)

tcr_bam %>%
  iwalk(filterBam, filter = filt_rules)

names(tcr_bam) %>%
  str_c(".bai") %>%
  file.remove()





### SAVE OBJECTS ----

# Add V(D)J data to Seurat object
vdj_so <- tiny_so %>%
  import_vdj(
    vdj_dir       = bcr_path,
    filter_chains = TRUE
  )

# Create tiny SingleCellExperiment object
tiny_sce <- SingleCellExperiment(
  list(counts = GetAssayData(tiny_so, "counts")),
  colData = tiny_so@meta.data
)

vdj_sce <- tiny_sce %>%
  import_vdj(
    vdj_dir       = bcr_path,
    filter_chains = TRUE
  )

# Save objects
use_data(tiny_so,  overwrite = TRUE)
use_data(vdj_so,   overwrite = TRUE)
use_data(tiny_sce, overwrite = TRUE)
use_data(vdj_sce,  overwrite = TRUE)
