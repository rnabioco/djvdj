library(dplyr)
library(tidyr)
library(purrr)
library(tibble)
library(stringr)
library(SingleCellExperiment)
library(djvdj)

### TESTING LARGE DATA ###
#
# # generating test data
# load("data/avid/so_avid.rda")
#
# test_vdj <- so_avid %>%
#   import_vdj("data/avid/bcr/")
#
# test_vdj <- test_vdj@meta.data
#
# walk(1:5, ~ {
#   test_vdj <<- bind_rows(test_vdj, test_vdj)
# })
#
# test_vdj <- test_vdj %>%
#   mutate(clonotype_id = str_c(clonotype_id, row_number(clonotype_id) %% 5))
#
# # testing functions
# test_vdj %>%
#   calc_abundance(cluster_col = "seurat_clusters")
#
# test_vdj %>%
#   plot_abundance(
#     cluster_col = "seurat_clusters",
#     type = "line",
#     n_clonotypes = 2
#   )
#
# test_vdj %>%
#   plot_similarity(
#     cluster_col = "seurat_clusters"
#   )
#
# test_vdj %>%
#   plot_diversity(
#     cluster_col = "seurat_clusters",
#   )
#
# clmns <- c(
#   "v_gene", "d_gene", "chains",
#   "umis", "reads", "cdr3_length",
#   "cdr3_nt_length", "productive",
#   "full_length"
# )
#
# tictoc::tic()
# x <- fetch_vdj(test_vdj, clmns)
# tictoc::toc()
#
# tictoc::tic()
# y <- summarize_vdj(test_vdj, clmns)
# tictoc::toc()
#
# tictoc::tic()
# x <- fetch_vdj(test_vdj, clonotype_col = NULL)
# tictoc::toc()
#
# tictoc::tic()
# y <- fetch_vdj(test_vdj, clonotype_col = "clonotype_id")
# tictoc::toc()

