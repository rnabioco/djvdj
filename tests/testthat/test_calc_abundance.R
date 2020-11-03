
test_that("Check all calc_abundance arguments", {
  args_df <- list(
    cluster_col   = list(NULL, "seurat_clusters"),
    prefix        = c("", "X"),
    return_seurat = c(TRUE, FALSE)
  ) %>%
    expand.grid()

  res <- pmap(
    args_df,
    check_args,
    .fn           = calc_abundance,
    sobj_in       = tiny_vdj,
    clonotype_col = "cdr3"
  )
})

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
      clonotype_col = "clonotype_id",
      return_seurat = FALSE
    )

  expect_s3_class(res, "tbl")
})

test_that("Check empty clonotype_col", {
  .fn<- function() {
    tiny_vdj %>%
      calc_abundance(return_seurat = TRUE)
  }

  expect_error(fn())
})
