
test_that("Single import with full path", {
  res <- tiny_so %>%
    import_vdj(
      vdj_dir = system.file("extdata/outs", package = "djvdj")
    )

  expect_s4_class(res, "Seurat")
  expect_identical(colnames(res), colnames(tiny_so))
})

test_that("Single import without full path", {
  res <- tiny_so %>%
    import_vdj(
      vdj_dir = system.file("extdata", package = "djvdj")
    )

  expect_s4_class(res, "Seurat")
  expect_identical(colnames(res), colnames(tiny_so))
})

test_that("Bad separator", {
  fn <- function() {
    res <- tiny_so %>%
      import_vdj(
        vdj_dir = system.file("extdata", package = "djvdj"),
        sep = "A"
      )
  }

  expect_error(fn())
})
