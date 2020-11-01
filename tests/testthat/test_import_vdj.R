
test_that("Single import with full path", {
  res <- tiny_so %>%
    import_vdj(vdj_dir = system.file("extdata/outs", package = "djvdj"))

  expect_s4_class(res, "Seurat")                       # class
  expect_identical(colnames(res), colnames(tiny_so))   # cells in object
})

test_that("Single import without full path", {
  res <- tiny_so %>%
    import_vdj(vdj_dir = system.file("extdata", package = "djvdj"))

  expect_s4_class(res, "Seurat")                       # class
  expect_identical(colnames(res), colnames(tiny_so))   # cells in object
})
