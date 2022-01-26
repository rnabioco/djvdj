# Test data
ctigs <- c(
  system.file("extdata/bcr_1", package = "djvdj"),
  system.file("extdata/bcr_2", package = "djvdj")
)

ctigs_2 <- ctigs_3 <- ctigs

names(ctigs_2) <- c("1", "2")

tcr_ctigs <- system.file("extdata/tcr_1", package = "djvdj")
bad_ctigs <- system.file("extdata/bad_bcr_1", package = "djvdj")

vdj_cols <- c(
  "v_gene",      "d_gene",
  "j_gene",      "c_gene",
  "chains",      "clonotype_id",
  "cdr3",        "cdr3_nt",
  "reads",       "umis",
  "productive",  "full_length",
  "n_insertion", "n_deletion",
  "n_mismatch"
)

df_1 <- tiny_so@meta.data

df_2 <- tiny_so@meta.data %>%
  as_tibble(rownames = ".cell_id")

# Check arguments for different path inputs with Seurat
arg_lst <- list(
  single = list(
    input          = list(tiny_so),
    vdj_dir        = list(ctigs_2[1], c("1_" = paste0(ctigs[1], "/outs"))),
    cell_prefix    = list(NULL, "1"),
    filter_paired  = c(TRUE, FALSE),
    include_indels = FALSE
  ),
  multi = list(
    input          = list(tiny_so),
    vdj_dir        = list(ctigs, ctigs_2, paste0(ctigs, "/outs")),
    cell_prefix    = list(NULL, c("1", "2")),
    filter_paired  = c(TRUE, FALSE),
    include_indels = FALSE
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

# Check CDR3 length calculation
test_that("import_vdj CDR3 lengths", {
  res <- tiny_so@meta.data %>%
    import_vdj(ctigs)

  dat    <- fetch_vdj(res)
  lens_1 <- as.double(nchar(dat$cdr3))
  lens_2 <- as.double(nchar(dat$cdr3_nt))

  expect_true(length(unique(dat$cdr3_length)) > 1)
  expect_true(length(unique(dat$cdr3_nt_length)) > 1)
  expect_identical(lens_1, dat$cdr3_length)
  expect_identical(lens_2, dat$cdr3_nt_length)
})

# Check number of chains represented in each column
test_that("import_vdj number of chains represented by each column", {
  res <- tiny_so@meta.data %>%
    import_vdj(ctigs)

  dat <- res %>%
    select(all_of(vdj_cols)) %>%
    select(-clonotype_id) %>%
    mutate(across(everything(), ~ lengths(gregexpr(";", .x))))

  dat %>%
    purrr::reduce(expect_identical)
})

# Check include_indels
test_that("import_vdj include_indels", {
  res <- tiny_so %>%
    import_vdj(
      vdj_dir = ctigs,
      prefix  = "PREFIX_"
    )

  dat <- res@meta.data
  dat <- dat[, !grepl("PREFIX_", colnames(dat))]

  expect_identical(dat, tiny_so@meta.data)
  expect_s4_class(res, "Seurat")
  expect_identical(colnames(res), colnames(tiny_so))

  res <- tiny_so %>%
    import_vdj(
      vdj_dir     = ctigs[1],
      cell_prefix = "1",
      prefix      = "PREFIX_"
    )

  dat <- res@meta.data
  dat <- dat[, !grepl("PREFIX_", colnames(dat))]

  expect_identical(dat, tiny_so@meta.data)
  expect_s4_class(res, "Seurat")
  expect_identical(colnames(res), colnames(tiny_so))
})

# Check arguments for SCE input
test_that("import_vdj SingleCellExperiment", {
  res <- tiny_sce %>%
    import_vdj(
      vdj_dir = ctigs,
      prefix  = "PREFIX_"
    )

  dat <- res@colData
  dat <- dat[, !grepl("PREFIX_", colnames(dat))]

  expect_identical(dat, tiny_sce@colData)
  expect_s4_class(res, "SingleCellExperiment")
  expect_identical(colnames(res), colnames(tiny_sce))
})

# Check arguments for data.frame input
test_that("import_vdj data.frame", {
  res <- df_1 %>%
    import_vdj(
      vdj_dir = ctigs,
      prefix  = "PREFIX_"
    )

  dat <- res[, !grepl("PREFIX_", colnames(res))]

  expect_identical(dat, df_1)
  expect_s3_class(res, "data.frame")
  expect_identical(rownames(res), rownames(df_1))

  res <- df_2 %>%
    import_vdj(
      vdj_dir = ctigs,
      prefix  = "PREFIX_"
    )

  dat <- res[, !grepl("PREFIX_", colnames(res))]

  expect_identical(dat, df_1)
  expect_identical(c(".cell_id", colnames(dat)), colnames(df_2))
  expect_s3_class(res, "data.frame")
  expect_identical(rownames(res), df_2$.cell_id)
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
test_that("import_vdj filter_chains", {
  res <- tiny_so %>%
    import_vdj(ctigs)

  expect_false(any(grepl("FALSE", res$productive)))
  expect_false(any(grepl("FALSE", res$full_length)))
  expect_s4_class(res, "Seurat")
  expect_identical(colnames(res), colnames(tiny_so))

  res <- tiny_so %>%
    import_vdj(
      vdj_dir        = ctigs[1],
      cell_prefix    = "1",
      include_indels = FALSE,
      filter_chains  = FALSE
    )

  expect_true(any(grepl("FALSE", res$productive)))
  expect_true(any(grepl("FALSE", res$full_length)))
  expect_s4_class(res, "Seurat")
  expect_identical(colnames(res), colnames(tiny_so))

  fn <- function() {
    tiny_so %>%
      import_vdj(
        vdj_dir       = ctigs,
        filter_chains = FALSE
      )
  }

  expect_warning(fn(), "When include_indels is TRUE, filter_chains is also automatically set TRUE")
})

# Check filter_paired
test_that("import_vdj filter_paired", {
  res <- tiny_so %>%
    import_vdj(
      vdj_dir       = ctigs,
      filter_paired = TRUE
    )

  expect_true(all(res$n_chains > 1, na.rm = TRUE))
  expect_true(all(res$paired, na.rm = TRUE))
  expect_s4_class(res, "Seurat")
  expect_identical(colnames(res), colnames(tiny_so))

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
test_that("import_vdj low overlap", {
  dat <- c("1" = bad_ctigs)

  fn <- function() {
    res <- tiny_so %>%
      import_vdj(
        vdj_dir = dat,
        include_indels = TRUE
      )
  }

  expect_warning(fn(), "Only.+cell barcodes overlap")
})

# Check missing clonotype_id warning
test_that("import_vdj missing clonotype_id", {
  fn <- function() {
    res <- tiny_so %>%
      import_vdj(
        vdj_dir        = bad_ctigs,
        cell_prefix    = "1",
        filter_chains  = FALSE,
        include_indels = TRUE
      )
  }

  expect_warning(fn(), "contigs do not have an assigned clonotype_id, these contigs will be removed")
  expect_warning(fn(), "Some chains are missing indel data")
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

  expect_error(fn(), "is already present in the input data", fixed = TRUE)
})

# Check duplicated cell barcode prefixes
test_that("import_vdj duplicate cell prefix", {
  prfxs    <- rep("1", 2)
  ctigs_in <- c(ctigs[1], ctigs[1])
  dat      <- purrr::set_names(ctigs_in, prfxs)

  fn <- function() {
    res <- tiny_so %>%
      import_vdj(vdj_dir = dat)
  }

  expect_warning(fn(), "Some cell barcode prefixes are duplicated")

  fn <- function() {
    res <- tiny_so %>%
      import_vdj(
        vdj_dir     = ctigs_in,
        cell_prefix = prfxs
      )
  }

  expect_warning(fn(), "Some cell barcode prefixes are duplicated")
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

  expect_error(fn(), "cell_prefix must be the same length as vdj_dir")
})

# Check cell barcode prefix NAs
test_that("import_vdj cell prefix NAs", {
  prfxs <- c(NA, "2")
  dat   <- purrr::set_names(ctigs, prfxs)

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

# Check BCR and TCR
test_that("import_vdj BCR and TCR", {
  fn <- function() {
    res <- tiny_so %>%
      import_vdj(
        vdj_dir        = c(tcr_ctigs, ctigs),
        include_indels = FALSE
      )
  }

  expect_error(fn(), "Multiple data types detected")
})

# Check .classify_vdj
test_that(".classify_vdj", {
  dat <- vdj_so %>%
    fetch_vdj(vdj_cols = "chains")

  expect_error(
    dat %>%
      mutate(chains = paste0("A", chains)) %>%
      .classify_vdj(),
    "None of the expected chains.+were found"
  )

  expect_error(
    dat %>%
      mutate(chains = paste0(chains, "A")) %>%
      .classify_vdj(),
    "None of the expected chains.+were found"
  )

  expect_warning(
    dat %>%
      filter(chains %in% c("IGH", "IGK")) %>%
      group_by(chains) %>%
      dplyr::slice(1:10) %>%
      ungroup() %>%
      mutate(chains = gsub("IGH", "TRA", chains)) %>%
      .classify_vdj(),
    "Equal number of BCR.+and TCR.+chains detected"
  )
})




