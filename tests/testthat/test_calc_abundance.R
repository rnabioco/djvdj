
# Check all calc_abundance arguments
lst <- list(
  sobj_in       = list(tiny_vdj),
  clonotype_col = "cdr3",
  cluster_col   = list(NULL, "seurat_clusters"),
  prefix        = c("", "X"),
  return_seurat = c(TRUE, FALSE)
)

test_all_args(
  lst    = lst,
  .fn    = calc_abundance,
  ttl    = "calc_abundance args",
  chk_fn = expect_silent
)

test_that("Check Seurat output", {
  res <- tiny_vdj %>%
    calc_abundance(
      clonotype_col = "cdr3",
      return_seurat = TRUE
    )

  expect_s4_class(res, "Seurat")
  expect_identical(colnames(res), colnames(tiny_vdj))
})

# Check tibble output
test_that("calc_abundance tbl out", {
  res <- tiny_vdj %>%
    calc_abundance(
      clonotype_col = "cdr3",
      return_seurat = FALSE
    )

  expect_s3_class(res, "tbl")
})

# Check empty clonotype_col
test_that("calc_abundance NULL clonotype_col", {
  .fn<- function() {
    tiny_vdj %>%
      calc_abundance(return_seurat = TRUE)
  }

  expect_error(fn())
})
