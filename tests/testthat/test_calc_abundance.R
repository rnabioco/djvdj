
test_that("Check Seurat output", {
  res <- tiny_vdj %>%
    calc_abundance(
      clonotype_col = "clonotype_id",
      return_seurat = TRUE
    )

  expect_s4_class(res, "Seurat")                       # class
  expect_identical(colnames(res), colnames(tiny_vdj))  # cells in object
})

test_that("Check tibble output", {
  res <- tiny_vdj %>%
    calc_abundance(
      return_seurat = FALSE
    )

  expect_s3_class(res, "tbl")
})

test_that("Check empty clonotype_col", {
  fn <- function() {
    tiny_vdj %>%
      calc_abundance(return_seurat = TRUE)
  }

  expect_error(fn())
})
