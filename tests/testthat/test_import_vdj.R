
test_that("Single import without names", {
  res <- tiny_so %>%
    import_vdj(vdj_dir = "")

  expect_s4_class(res, "Seurat")                       # class
  expect_identical(colnames(res), colnames(tiny_so))   # cells in object
})
