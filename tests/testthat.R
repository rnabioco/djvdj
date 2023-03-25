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





# ### TESTING LARGE DATA ###
# 
# # Generate test data
# n_reps <- 10
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
