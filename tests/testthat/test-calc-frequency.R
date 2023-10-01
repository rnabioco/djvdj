# Test data
df_1 <- vdj_sce@colData

df_2 <- vdj_sce@colData |>
  as.data.frame() |>
  as_tibble(rownames = ".cell_id")

# Check all calc_frequency arguments
arg_lst <- list(
  input       = list(vdj_sce, df_1, df_1),
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

test_that("Check SCE output", {
  res <- vdj_sce |>
    calc_frequency(
      data_col  = "cdr3",
      return_df = FALSE
    )

  expect_s4_class(res, "SingleCellExperiment")
  expect_identical(colnames(res), colnames(vdj_sce))
})

# Check data.frame output
test_that("calc_frequency df out", {
  res <- vdj_sce |>
    calc_frequency(
      data_col  = "cdr3",
      return_df = TRUE
    )

  expect_s3_class(res, "data.frame")
})

# Check data.frame input
test_that("calc_frequency df in", {
  res <- vdj_sce@colData |>
    calc_frequency(data_col = "cdr3")

  expect_s3_class(res, "data.frame")
})

