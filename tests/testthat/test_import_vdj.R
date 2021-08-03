
# Contig files for test
ctigs <- c(
  system.file("extdata/bcr_1", package = "djvdj"),
  system.file("extdata/bcr_2", package = "djvdj")
)

test_that("Single import with full path", {
  res <- tiny_so %>%
    import_vdj(vdj_dir = file.path(ctigs[1], "outs"))

  expect_s4_class(res, "Seurat")
  expect_identical(colnames(res), colnames(tiny_so))
})

test_that("Single import without full path", {
  res <- tiny_so %>%
    import_vdj(vdj_dir = ctigs[1])

  expect_s4_class(res, "Seurat")
  expect_identical(colnames(res), colnames(tiny_so))
})

test_that("Multiple named import", {
  dat <- ctigs

  names(dat) <- c("", "2_")

  res <- tiny_so %>%
    import_vdj(vdj_dir = dat)

  expect_s4_class(res, "Seurat")
  expect_identical(colnames(res), colnames(tiny_so))
})

test_that("Filtered contigs", {
  res <- tiny_so %>%
    import_vdj(
      vdj_dir        = ctigs,
      filter_contigs = TRUE
    )

  expect_false(any(grepl("FALSE", res$productive)))
  expect_false(any(grepl("FALSE", res$full_length)))
  expect_s4_class(res, "Seurat")
  expect_identical(colnames(res), colnames(tiny_so))
})

test_that("Unfiltered contigs", {
  suppressWarnings({
    res <- tiny_so %>%
      import_vdj(
        vdj_dir        = ctigs,
        filter_contigs = FALSE
      )
  })

  expect_true(any(grepl("FALSE", res$productive)))
  expect_true(any(grepl("FALSE", res$full_length)))
  expect_s4_class(res, "Seurat")
  expect_identical(colnames(res), colnames(tiny_so))
})

test_that("Multiple unnamed import", {
  res <- tiny_so %>%
    import_vdj(ctigs)

  expect_s4_class(res, "Seurat")
  expect_identical(colnames(res), colnames(tiny_so))
})

test_that("tibble output", {
  res <- import_vdj(vdj_dir = ctigs)

  expect_s3_class(res, "tbl")
})

test_that("Bad cell barcode prefixes", {
  dat <- setNames(ctigs, c("A", "B"))

  fn <- function() {
    res <- tiny_so %>%
      import_vdj(vdj_dir = dat, filter_contigs = FALSE)
  }

  expect_error(fn())
})

test_that("Low overlap warning", {
  dat <- c("2_" = ctigs[2])

  fn <- function() {
    res <- tiny_so %>%
      import_vdj(vdj_dir = dat)
  }

  expect_warning(fn())
})

test_that("Missing clonotype_id warning", {
  fn <- function() {
    res <- tiny_so %>%
      import_vdj(
        vdj_dir        = ctigs,
        filter_contigs = FALSE
      )
  }

  expect_warning(fn())
})

test_that("Bad separator", {
  fn <- function() {
    res <- tiny_so %>%
      import_vdj(
        vdj_dir = ctigs,
        sep     = "A"
      )
  }

  expect_error(fn())
})
