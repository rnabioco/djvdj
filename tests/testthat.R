library(testthat)
library(dplyr)
library(tidyr)
library(purrr)
library(tibble)
library(stringr)
library(SingleCellExperiment)
library(djvdj)

vdj_dirs <- c(
  BL6 = system.file("extdata/splen/BL6_BCR", package = "djvdj"),
  MD4 = system.file("extdata/splen/MD4_BCR", package = "djvdj")
)

test_so <- splen_so |>
  import_vdj(vdj_dirs, define_clonotypes = "cdr3_gene")

test_check("djvdj")





### TESTING LARGE DATA ###

# # Generate test data
# n_reps <- 4
#
# vdj_dirs <- c(
#   BL6 = system.file("extdata/splen/BL6_BCR", package = "djvdj"),
#   MD4 = system.file("extdata/splen/MD4_BCR", package = "djvdj")
# )
#
# test_vdj <- splen_so |>
#   import_vdj(vdj_dirs, define_clonotypes = "cdr3_gene")
#
# test_vdj <- test_vdj@meta.data
#
# walk(seq_len(n_reps), ~ {
#   test_vdj <<- bind_rows(test_vdj, test_vdj)
# })
#
# test_vdj <- test_vdj %>%
#   mutate(clonotype_id = str_c(clonotype_id, row_number(clonotype_id) %% n_reps))

# # testing functions
# usage_args <- list(
#   input       = test_vdj,
#   # data_cols   = c("v_gene", "j_gene"),
#   data_cols   = c("v_gene"),
#   # cluster_col = c("sample"),
#   chain       = "IGK"
# )
#
# # profvis::profvis({
#   tic()
#   x <- purrr::lift_dl(calc_gene_usage)(usage_args)
#   # x <- purrr::lift_dl(calc_gene_usage, return_df = TRUE)(usage_args)
#   toc()
# # })
#
# # profvis::profvis({
#   tic()
#   y <- purrr::lift_dl(old_calc_gene_usage)(usage_args)
#   toc()
# # })
#
# identical(
#   y %>% arrange(v_gene, n_cells, freq, pct),
#   x %>% select(-any_of("shared")) %>% arrange(v_gene, n_cells, freq, pct)
# )
#
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

