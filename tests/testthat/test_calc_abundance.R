
# Check all calc_abundance arguments
arg_lst <- list(
  input         = list(tiny_vdj),
  clonotype_col = "cdr3",
  cluster_col   = list(NULL, "seurat_clusters"),
  prefix        = c("", "X"),
  return_df     = c(TRUE, FALSE)
)

test_all_args(
  arg_lst = arg_lst,
  .fn     = calc_abundance,
  desc    = "calc_abundance args",
  chk     = expect_silent
)

test_that("Check Seurat output", {
  res <- tiny_vdj %>%
    calc_abundance(
      clonotype_col = "cdr3",
      return_df     = FALSE
    )

  expect_s4_class(res, "Seurat")
  expect_identical(colnames(res), colnames(tiny_vdj))
})

# Check data.frame output
test_that("calc_abundance df out", {
  res <- tiny_vdj %>%
    calc_abundance(
      clonotype_col = "cdr3",
      return_df     = TRUE
    )

  expect_s3_class(res, "data.frame")
})

# Check data.frame input
test_that("calc_abundance df in", {
  res <- tiny_vdj@meta.data %>%
    calc_abundance(clonotype_col = "cdr3")

  expect_s3_class(res, "data.frame")
})

# Check empty clonotype_col
test_that("calc_abundance NULL clonotype_col", {
  .fn<- function() {
    tiny_vdj %>%
      calc_abundance(return_df = FALSE)
  }

  expect_error(fn())
})
