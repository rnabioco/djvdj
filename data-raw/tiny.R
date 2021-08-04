library(tidyverse)
library(Seurat)
library(djvdj)

# Parameters
dat_dir <- "~/Projects/Smith_AVIDseq"

mat_path <- c(
  str_c(dat_dir, "/results/JH179_GEX-JH181_ADT_JH181_HTO/outs/filtered_feature_bc_matrix"),
  str_c(dat_dir, "/2020-07-17/JH191_GEX/outs/filtered_feature_bc_matrix")
)

vdj_path <- c(
  str_c(dat_dir, "/results/JH180_BCR/outs"),
  str_c(dat_dir, "/2020-07-17/BCR/outs")
)

names(vdj_path) <- c("", "2_")

# Create Suerat object
mat <- Read10X(mat_path)

tiny_so <- mat$`Gene Expression` %>%
  CreateSeuratObject() %>%
  mutate_meta(
    mutate,
    orig.ident = ifelse(str_detect(.cell_id, "^2_"), "avid_2", "avid_1")
  )

# Subset object
# keep cell order the same
all_cells <- colnames(tiny_so)

tiny_cells <- all_cells %>%
  sample(200)

tiny_cells <- all_cells[all_cells %in% tiny_cells]

tiny_so <- tiny_so %>%
  subset(
    features = sample(rownames(tiny_so), 200),
    cells    = tiny_cells
  )

# Normalize object
tiny_so <- tiny_so %>%
  PercentageFeatureSet(
    pattern  = "^mt-",
    col.name = "Percent_mito"
  ) %>%
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

# Format contig data
# add some NAs to raw_clonotype_id
contigs <- vdj_path %>%
  map(str_c, "/filtered_contig_annotations.csv") %>%
  map(read_csv)

contigs <- contigs %>%
  imap(~ filter(.x, str_c(.y, barcode) %in% tiny_cells))

contigs[[2]] <- contigs[[2]] %>%
  mutate(
    raw_clonotype_id = ifelse(contig_id == "ACGAGCCTCGGATGTT-1_contig_3", NA, raw_clonotype_id)
  )

names(contigs) <- str_c("bcr_", 1:2) %>%
  str_c("inst/extdata/", ., "/outs/filtered_contig_annotations.csv")

contigs %>%
  iwalk(write_csv)

# Add V(D)J data to object
tiny_vdj <- tiny_so %>%
  import_vdj(
    vdj_dir        = vdj_path,
    filter_contigs = TRUE
  )

# Save objects
usethis::use_data(
  tiny_so,
  overwrite = TRUE
)

usethis::use_data(
  tiny_vdj,
  overwrite = TRUE
)
