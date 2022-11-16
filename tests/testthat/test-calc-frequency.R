# Test data
df_1 <- vdj_so@meta.data

df_2 <- vdj_so@meta.data |>
  as_tibble(rownames = ".cell_id")

# Check all calc_frequency arguments
arg_lst <- list(
  input       = list(vdj_so, vdj_sce, df_1, df_1),
  data_col    = "cdr3",
  cluster_col = list(NULL, "seurat_clusters"),
  prefix      = c("", "X"),
  return_df   = c(TRUE, FALSE)
)

test_all_args(
  arg_lst = arg_lst,
  .fn     = calc_frequency,
  desc    = "calc_frequency args",
  chk     = expect_silent
)

test_that("Check Seurat output", {
  res <- vdj_so |>
    calc_frequency(
      data_col = "cdr3",
      return_df     = FALSE
    )

  expect_s4_class(res, "Seurat")
  expect_identical(colnames(res), colnames(vdj_so))
})

# Check data.frame output
test_that("calc_frequency df out", {
  res <- vdj_so |>
    calc_frequency(
      data_col = "cdr3",
      return_df     = TRUE
    )

  expect_s3_class(res, "data.frame")
})

# Check data.frame input
test_that("calc_frequency df in", {
  res <- vdj_so@meta.data |>
    calc_frequency(data_col = "cdr3")

  expect_s3_class(res, "data.frame")
})

