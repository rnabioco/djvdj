library(tidyverse)
library(Seurat)
library(clustifyr)
library(clustifyrdata)
library(djvdj)

# Parameters
mat_path <- "~/Projects/Smith_AVIDseq/2020-07-17/JH191_GEX/outs/filtered_feature_bc_matrix"
vdj_path <- "~/Projects/Smith_AVIDseq/2020-07-17/BCR/outs"
gene_min <- 100
gene_max <- 5000
mito_max <- 15

hto_ids <- list(
  BL6 = "HTO1",
  MD4 = "HTO2"
)

hto_lims <- list(
  BL6 = 0.3,
  MD4 = 2
)

# Create Seurat object using gene expression data
mat_list <- Read10X(mat_path)
rna_mat <- mat_list[["Gene Expression"]]

so <- rna_mat %>%
  CreateSeuratObject(
    project = "AVID-seq",
    min.cells = 5
  )

# Create ADT matrix
adt_mat <- so %>%
  colnames() %>%
  mat_list[["Antibody Capture"]][, .]

# Create HTO matrix and set names
hto_mat <- adt_mat[rownames(adt_mat) %in% hto_ids, , drop = F]   # Selected HTO antibodies
adt_mat <- adt_mat[!rownames(adt_mat) %in% hto_ids, , drop = F]  # Removed HTO data from ADT matrix

hto_names <- set_names(names(hto_ids), hto_ids)
rownames(hto_mat) <- hto_names[rownames(hto_mat)]

# Add ADT and HTO assays
so[["ADT"]] <- CreateAssayObject(adt_mat)
so[["HTO"]] <- CreateAssayObject(hto_mat)

# Normalize object
so <- so %>%
  NormalizeData(assay = "RNA") %>%
  NormalizeData(
    assay = "ADT",
    normalization.method = "CLR"
  ) %>%
  NormalizeData(
    assay = "HTO",
    normalization.method = "CLR"
  ) %>%
  FindVariableFeatures(assay = "RNA") %>%
  ScaleData(assay = "RNA")

# Divide into BL6 and MD4 based on HTO signal
so@meta.data <- so %>%
  AddMetaData(metadata = FetchData(., c("hto_BL6", "hto_MD4"))) %>%
  .@meta.data %>%
  as_tibble(rownames = "cell_id") %>%
  mutate(
    mouse = case_when(
      hto_MD4 > hto_lims$MD4 & hto_BL6 < hto_lims$BL6 ~ "MD4",
      hto_BL6 > hto_lims$BL6 & hto_MD4 < hto_lims$MD4 ~ "BL6",
      hto_MD4 > hto_lims$MD4 & hto_BL6 > hto_lims$BL6 ~ "dbl_pos",
      TRUE                                            ~ "dbl_neg"
    )
  ) %>%
  column_to_rownames("cell_id")

# Filter object
so <- so %>%
  PercentageFeatureSet(
    pattern  = "^mt-",
    col.name = "Percent_mito"
  ) %>%
  subset(
    nFeature_RNA > gene_min &
    nFeature_RNA < gene_max &
    Percent_mito < mito_max &
    mouse %in% names(hto_lims)
  )

# Run PCA, cluster, run UMAP
so <- so %>%
  RunPCA(assay = "RNA") %>%
  FindNeighbors(
    assay     = "RNA",
    reduction = "pca",
    dims      = 1:40
  ) %>%
  FindClusters(
    resolution = 0.3,
    verbose = F
  ) %>%
  RunUMAP(
    assay = "RNA",
    dims = 1:40
  )

# Clustify cells
so <- so %>%
  clustify(
    ref_mat = ref_tabula_muris_drop,
    cluster_col = "seurat_clusters"
  )

so@meta.data <- so@meta.data %>%
  as_tibble(rownames = "cell_id") %>%
  mutate(
    type = if_else(
      str_detect(type, "B cell|T cell"),
      str_extract(type, "B cell|T cell"),
      type
    ),
    type_mouse = if_else(str_detect(type, "^B cell$|^T cell$"), type, "other"),
    type_mouse = str_c(type_mouse, " ", mouse)
  ) %>%
  mutate(type = str_remove(type, "-[a-zA-Z]+$")) %>%
  column_to_rownames("cell_id")

# Test data
tiny_cells <- colnames(so) %>%
  sample(100)

tiny_so <- so %>%
  subset(
    features = VariableFeatures(so)[1:100],
    cells = tiny_cells
  )

tiny_vdj <- tiny_so %>%
  import_vdj(
    vdj_dir = vdj_path,
    filter_contigs = TRUE
  )

# Test contigs
contigs <- read_csv(file.path(vdj_path, "filtered_contig_annotations.csv"))

contigs %>%
  filter(barcode %in% tiny_cells) %>%
  write_csv("inst/outs/filtered_contig_annotations.csv")

# Save objects
usethis::use_data(
  tiny_so,
  overwrite = TRUE
)

usethis::use_data(
  tiny_vdj,
  overwrite = TRUE
)

