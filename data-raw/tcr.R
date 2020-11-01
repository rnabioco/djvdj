library(tidyverse)
library(Seurat)
library(here)


# Load matrices
data_dir <- "~/Projects/Rincon_scVDJseq/results"

samples <- c(
  "WT_DN3", "WT_DN4",
  "KI_DN3", "KI_DN4"
) %>%
  str_c("_GE")

paths <- file.path(
  data_dir, samples,
  "outs/filtered_feature_bc_matrix"
)

names(paths) <- samples

mat <- Read10X(paths)

# Create Seurat object
so_tcr <- mat %>%
  CreateSeuratObject(
    min.cells    = 3,
    min.features = 200,
    names.delim  = "_",
    names.field  = 1:2
  )

so_tcr$genotype <- str_split(so_tcr$orig.ident, "_", simplify = T)[, 1]
so_tcr$stage    <- str_split(so_tcr$orig.ident, "_", simplify = T)[, 2]

# Filter Seurat object
so_tcr <- so_tcr %>%
  PercentageFeatureSet(
    pattern  = "^mt-",
    col.name = "percent.mt"
  )

so_tcr <- so_tcr %>%
  subset(
    percent.mt < 7.5 &
    nCount_RNA > 500 &
    nCount_RNA < 50000 &
    nFeature_RNA > 500
  )

# Normalize counts
so_tcr <- so_tcr %>%
  NormalizeData() %>%

  FindVariableFeatures(
    selection.method = "vst",
    nfeatures        = 3000,
    verbose          = FALSE
  ) %>%

  ScaleData(verbose = FALSE)

# Dimensional reduction
seed_val <- 20200901

so_tcr <- so_tcr %>%
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
  ) %>%

  AddMetaData(
    metadata = FetchData(., c("UMAP_1", "UMAP_2")),
    col.name = c("UMAP_1", "UMAP_2")
  )

# Clustering
so_tcr <- so_tcr %>%
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

# Save object
save(so_tcr, file = here("data/so_tcr.rda"))
