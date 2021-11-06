# Test data
ctigs <- c(
  system.file("extdata/bcr_1", package = "djvdj"),
  system.file("extdata/bcr_2", package = "djvdj")
)

ctigs_2 <- ctigs_3 <- ctigs

names(ctigs_2) <- c("", "2_")
names(ctigs_3) <- c("", "2")

vdj_cols <- c(
  "v_gene", "d_gene",
  "j_gene", "c_gene",
  "chains", "clonotype_id",
  "cdr3",   "cdr3_nt",
  "reads",  "productive",
  "umis",   "full_length"
)

df_1 <- vdj_so@meta.data

df_2 <- vdj_so@meta.data %>%
  as_tibble(rownames = ".cell_id")

# Check arguments for different path inputs with Seurat
arg_lst <- list(
  single = list(
    input         = list(tiny_so),
    vdj_dir       = list(ctigs[1], ctigs_2[1], ctigs_3[1], paste0(ctigs[1], "/outs")),
    cell_prefix   = list(NULL, ""),
    filter_paired = c(TRUE, FALSE),
    prefix        = c("", "PREFIX")
  ),
  multi = list(
    input         = list(tiny_so),
    vdj_dir       = list(ctigs, ctigs_2, ctigs_3, paste0(ctigs, "/outs")),
    cell_prefix   = list(NULL, c("", "2_"), c("", "2")),
    filter_paired = c(TRUE, FALSE),
    prefix        = c("", "PREFIX")
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

# Check arguments for different path inputs with SCE
pluck(arg_lst, 1, 1) <- list(tiny_sce)
pluck(arg_lst, 2, 1) <- list(tiny_sce)

arg_lst %>%
  iwalk(~ {
    test_all_args(
      arg_lst = .x,
      .fn     = import_vdj,
      desc    = paste("import_vdj", .y, "path class"),
      chk     = expr(expect_s4_class(.res, "SingleCellExperiment"))
    )

    test_all_args(
      arg_lst = .x,
      .fn     = import_vdj,
      desc    = paste("import_vdj", .y, "path cells"),
      chk     = expr(expect_identical(colnames(.res), colnames(tiny_so)))
    )
  })

# Check arguments for different path inputs with data.frame
pluck(arg_lst, 1, 1) <- list(df_1)
pluck(arg_lst, 2, 1) <- list(df_1)

arg_lst %>%
  iwalk(~ {
    test_all_args(
      arg_lst = .x,
      .fn     = import_vdj,
      desc    = paste("import_vdj", .y, "path class"),
      chk     = expr(expect_s3_class(.res, "data.frame"))
    )

    test_all_args(
      arg_lst = .x,
      .fn     = import_vdj,
      desc    = paste("import_vdj", .y, "path cells"),
      chk     = expr(expect_identical(rownames(.res), rownames(vdj_so@meta.data)))
    )
  })

# Check arguments for different path inputs with tibble
pluck(arg_lst, 1, 1) <- list(df_2)
pluck(arg_lst, 2, 1) <- list(df_2)

arg_lst %>%
  iwalk(~ {
    test_all_args(
      arg_lst = .x,
      .fn     = import_vdj,
      desc    = paste("import_vdj", .y, "path class"),
      chk     = expr(expect_s3_class(.res, "data.frame"))
    )

    test_all_args(
      arg_lst = .x,
      .fn     = import_vdj,
      desc    = paste("import_vdj", .y, "path cells"),
      chk     = expr(expect_identical(rownames(.res), rownames(vdj_so@meta.data)))
    )
  })

# Check bad path
test_that("import_vdj bad path", {
  fn <- function() {
    res <- tiny_so %>%
      import_vdj(vdj_dir = "BAD_PATH")
  }

  expect_error(fn(), "filtered_contig_annotations.csv not found in BAD_PATH")
})

# Check column prefix
test_that("import_vdj column prefix", {
  prfx     <- "TEST_"
  new_cols <- paste0(prfx, vdj_cols)

  res <- tiny_so %>%
    import_vdj(
      vdj_dir = ctigs,
      prefix  = prfx
    )

  expect_false(any(vdj_cols %in% colnames(res@meta.data)))
  expect_true(all(new_cols %in% colnames(res@meta.data)))
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

# Check only paired chains
test_that("import_vdj only paired chains", {
  res <- tiny_so %>%
    import_vdj(
      vdj_dir       = ctigs,
      filter_paired = TRUE
    )

  expect_true(all(res$paired, na.rm = TRUE))
  expect_s4_class(res, "Seurat")
  expect_identical(colnames(res), colnames(tiny_so))
})

# Check include unpaired chains
test_that("import_vdj include unpaired chains", {
  res <- tiny_so %>%
    import_vdj(
      vdj_dir       = ctigs,
      filter_paired = FALSE
    )

  expect_true(any(!res$paired, na.rm = TRUE))
  expect_s4_class(res, "Seurat")
  expect_identical(colnames(res), colnames(tiny_so))
})

# Check define_clonotypes options
test_that("import_vdj define_clonotypes options", {
  opts <- c(cdr3 = "cdr3aa", cdr3_nt = "cdr3nt")

  opts %>%
    iwalk(~ {
      res <- tiny_so %>%
        import_vdj(
          vdj_dir = ctigs,
          define_clonotypes = .x
        )

      expect_identical(n_distinct(res$clonotype_id), n_distinct(res[[.y]]))
      expect_s4_class(res, "Seurat")
      expect_identical(colnames(res), colnames(tiny_so))
    })

  res <- tiny_so %>%
    import_vdj(
      vdj_dir = ctigs,
      define_clonotypes = "cdr3_gene"
    )

  x <- n_distinct(res$clonotype_id)
  y <- nrow(distinct(res@meta.data, cdr3_nt, v_gene, d_gene, j_gene))

  expect_identical(x, y)
  expect_s4_class(res, "Seurat")
  expect_identical(colnames(res), colnames(tiny_so))
})

# Check define_clonotypes bad define_clonotypes
test_that("import_vdj bad define_clonotypes", {
  expect_error(
    tiny_so %>% import_vdj(vdj_dir = ctigs, define_clonotypes = "BAD"),
    "define_clonotypes must be one of"
  )
})

# Check data.frame output
test_that("import_vdj df out", {
  res <- import_vdj(vdj_dir = ctigs)

  expect_s3_class(res, "data.frame")
})

# Check bad barcode prefixes
test_that("import_vdj bad prefixes", {
  dat <- setNames(ctigs, c("A", "B"))

  fn <- function() {
    res <- tiny_so %>%
      import_vdj(vdj_dir = dat)
  }

  expect_error(fn(), "barcodes do not match those in the object, are you using the correct cell barcode prefixes?")
})

# Check low overlap warning
# test_that("import_vdj low overlap", {
#   dat <- c("2_" = ctigs[2])
#
#   fn <- function() {
#     res <- tiny_so %>%
#       import_vdj(vdj_dir = ctigs)
#   }
#
#   expect_warning(fn())
# })

# Check missing clonotype_id warning
test_that("import_vdj missing clonotype_id", {
  fn <- function() {
    res <- tiny_so %>%
      import_vdj(
        vdj_dir        = ctigs,
        filter_contigs = FALSE
      )
  }

  expect_warning(fn(), "contigs do not have an assigned clonotype_id, these contigs will be removed")
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

  expect_error(fn(), "is already present in the V(D)J data, select a different value for sep", fixed = TRUE)
})

# Check duplicated cell barcode prefixes
# test_that("import_vdj duplicate cell prefix", {
#   prfxs <- rep("", 2)
#   dat   <- set_names(ctigs, prfxs)
#
#   fn <- function() {
#     res <- tiny_so %>%
#       import_vdj(vdj_dir = dat)
#   }
#
#   expect_warning(fn())
#
#   fn <- function() {
#     res <- tiny_so %>%
#       import_vdj(
#         vdj_dir     = ctigs,
#         cell_prefix = prfxs
#       )
#   }
#
#   expect_warning(fn())
# })

# Check barcode prefix length
test_that("import_vdj cell prefix length", {
  fn <- function() {
    res <- tiny_so %>%
      import_vdj(
        vdj_dir     = ctigs,
        cell_prefix = "A"
      )
  }

  expect_error(fn(), "cell_prefix must be the same length as vdj_dir")
})

# Check cell barcode prefix NAs
test_that("import_vdj cell prefix NAs", {
  prfxs <- c(NA, "2")
  dat   <- set_names(ctigs, prfxs)

  fn <- function() {
    res <- tiny_so %>%
      import_vdj(vdj_dir = dat)
  }

  expect_error(fn(), "Cell prefixes cannot include NAs.")

  fn <- function() {
    res <- tiny_so %>%
      import_vdj(
        vdj_dir     = ctigs,
        cell_prefix = prfxs
      )
  }

  expect_error(fn(), "Cell prefixes cannot include NAs.")
})
