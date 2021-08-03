
# Test inputs
ctigs <- c(
  system.file("extdata/bcr_1", package = "djvdj"),
  system.file("extdata/bcr_2", package = "djvdj")
)

ctigs_2 <- ctigs_3 <- ctigs

names(ctigs_2) <- c("", "2_")
names(ctigs_3) <- c("", "2")

# Check arguments for different path inputs
arg_lst <- list(
  single = list(
    sobj_in     = list(tiny_so),
    vdj_dir     = list(ctigs[1], ctigs_2[1], ctigs_3[1], paste0(ctigs[1], "/outs")),
    cell_prefix = list(NULL, ""),
    prefix      = c("", "PREFIX")
  ),
  multi = list(
    sobj_in     = list(tiny_so),
    vdj_dir     = list(ctigs, ctigs_2, ctigs_3, paste0(ctigs, "/outs")),
    cell_prefix = list(NULL, c("", "2_"), c("", "2")),
    prefix      = c("", "PREFIX")
  )
)

arg_lst %>%
  iwalk(~ {
    test_all_args(
      arg_lst = .x,
      .fn     = import_vdj,
      desc    = paste("import_vdj", .y, "path class"),
      chk     = expr(expect_s4_class(.res, "Seurat"))
    )

    test_all_args(
      arg_lst = .x,
      .fn     = import_vdj,
      desc    = paste("import_vdj", .y, "path cells"),
      chk     = expr(expect_identical(colnames(.res), colnames(tiny_so)))
    )
  })

# Check bad path
test_that("import_vdj bad path", {
  fn <- function() {
    res <- tiny_so %>%
      import_vdj(vdj_dir = "BAD_PATH")
  }

  expect_error(fn())
})

# Check filtered contigs
test_that("import_vdj filtered contigs", {
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

# Check unfiltered contigs
test_that("import_vdj unfiltered contigs", {
  res <- tiny_so %>%
    import_vdj(
      vdj_dir        = ctigs[1],
      filter_contigs = FALSE
    )

  expect_true(any(grepl("FALSE", res$productive)))
  expect_true(any(grepl("FALSE", res$full_length)))
  expect_s4_class(res, "Seurat")
  expect_identical(colnames(res), colnames(tiny_so))
})

# Check tibble output
test_that("import_vdj tibble output", {
  res <- import_vdj(vdj_dir = ctigs)

  expect_s3_class(res, "tbl")
})

# Check bad barcode prefixes
test_that("import_vdj bad prefixes", {
  dat <- setNames(ctigs, c("A", "B"))

  fn <- function() {
    res <- tiny_so %>%
      import_vdj(vdj_dir = dat)
  }

  expect_error(fn())
})

# Check low overlap warning
test_that("import_vdj low overlap", {
  dat <- c("2_" = ctigs[2])

  fn <- function() {
    res <- tiny_so %>%
      import_vdj(vdj_dir = dat)
  }

  expect_warning(fn())
})

# Check missing clonotype_id warning
test_that("import_vdj missing clonotype_id", {
  fn <- function() {
    res <- tiny_so %>%
      import_vdj(
        vdj_dir        = ctigs,
        filter_contigs = FALSE
      )
  }

  expect_warning(fn())
})

# Check bad separator
test_that("import_vdj bad sep", {
  fn <- function() {
    res <- tiny_so %>%
      import_vdj(
        vdj_dir = ctigs,
        sep     = "A"
      )
  }

  expect_error(fn())
})

# Check duplicated cell barcode prefixes
test_that("import_vdj duplicate cell prefix", {
  prfxs <- rep("", 2)
  dat   <- set_names(ctigs, prfxs)

  fn <- function() {
    res <- tiny_so %>%
      import_vdj(vdj_dir = dat)
  }

  expect_warning(fn())

  fn <- function() {
    res <- tiny_so %>%
      import_vdj(
        vdj_dir     = ctigs,
        cell_prefix = prfxs
      )
  }

  expect_warning(fn())
})

# Check barcode prefix length
test_that("import_vdj cell prefix length", {
  fn <- function() {
    res <- tiny_so %>%
      import_vdj(
        vdj_dir     = ctigs,
        cell_prefix = "A"
      )
  }

  expect_error(fn())
})

# Check cell barcode prefix NAs
test_that("import_vdj cell prefix NAs", {
  prfxs <- c(NA, "2")
  dat   <- set_names(ctigs, prfxs)

  fn <- function() {
    res <- tiny_so %>%
      import_vdj(vdj_dir = dat)
  }

  expect_error(fn())

  fn <- function() {
    res <- tiny_so %>%
      import_vdj(
        vdj_dir     = ctigs,
        cell_prefix = prfxs
      )
  }

  expect_error(fn())
})








