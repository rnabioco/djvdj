library(tidyverse)
library(Seurat)
library(clustifyr)
library(clustifyrdata)
library(here)


# Load matrices
data_dir <- "~/Projects/Pelanda_scRNAseq/results"

samples <- c(
  "1_ABA_015_Expression",
  "2_ABCtr_015_Expression",
  "1_ABA_019_GEX",
  "2_ABCtr_019_GEX",
  "3_ABA_07_GEX",
  "4_ABCtr_07_GEX"
)

paths <- file.path(
  data_dir, samples,
  "outs/filtered_feature_bc_matrix"
)

names(paths) <- samples %>%
  str_remove("_GEX$|_Expression$")

mat <- Read10X(paths)

# Create Seurat object
so_bcr <- mat %>%
  CreateSeuratObject(
    min.cells    = 3,
    min.features = 200,
    names.delim  = "_",
    names.field  = 1:3
  )

# Filter Seurat object
so_bcr <- so_bcr %>%
  PercentageFeatureSet(
    pattern  = "^MT-",
    col.name = "percent.mt"
  )

so_bcr <- so_bcr %>%
  subset(
    percent.mt < 20 &
    nCount_RNA > 500 &
    nCount_RNA < 50000 &
    nFeature_RNA > 500
  )

# Normalize counts
so_bcr <- so_bcr %>%
  NormalizeData() %>%

  FindVariableFeatures(
    selection.method = "vst",
    nfeatures        = 3000,
    verbose          = FALSE
  ) %>%

  ScaleData(verbose = FALSE)

# Dimensional reduction
seed_val <- 20201021

so_bcr <- so_bcr %>%
  RunPCA(
    npcs     = 40,
    seed.use = seed_val,
    verbose  = FALSE
  ) %>%

  RunUMAP(
    reduction   = "pca",
    dims        = 1:30,
    n.neighbors = 20L,
    min.dist    = 0.5,
    seed.use    = seed_val,
    verbose     = FALSE
  )

# Clustering
so_bcr <- so_bcr %>%
  FindNeighbors(
    reduction = "pca",
    dims      = 1:30,
    k.param   = 20L,
    verbose   = FALSE
  ) %>%

  FindClusters(
    resolution  = 0.3,
    random.seed = seed_val,
    verbose     = FALSE
  )

# Assign cell types
so_bcr <- so_bcr %>%
  clustify(
    cluster_col = "seurat_clusters",
    ref_mat     = ref_hema_microarray
  )

# Save object
save(so_bcr, file = here("data/so_bcr.rda"))
