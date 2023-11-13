# Test data
# Contig paths
ctigs <- c(
  file.path(data_dir, "BL6_BCR"),
  file.path(data_dir, "MD4_BCR")
)

ctigs_2 <- ctigs_3 <- ctigs

names(ctigs_2) <- c("BL6", "MD4")

tcr_ctigs <- c(
  file.path(data_dir, "BL6_TCR"),
  file.path(data_dir, "MD4_TCR")
)

vdj_cols <- c(
  "v_gene",      "d_gene",
  "j_gene",      "c_gene",
  "chains",      "clonotype_id",
  "cdr3",        "cdr3_nt",
  "reads",       "umis",
  "productive",  "full_length"
)

df_1 <- splen_meta

df_2 <- splen_meta |>
  as_tibble(rownames = ".cell_id")

# Check arguments for different path inputs with Seurat
arg_lst <- list(
  input          = list(so),
  vdj_dir        = list(ctigs, ctigs_2),
  filter_paired  = c(TRUE, FALSE),
  include_mutations = FALSE
)

arg_lst |>
  test_all_args(
    .fn     = import_vdj,
    desc    = paste("import_vdj multi path class"),
    chk     = expr(expect_s4_class(.res, "Seurat"))
  )

arg_lst |>
  test_all_args(
    .fn     = import_vdj,
    desc    = paste("import_vdj multi path cells"),
    chk     = expr(expect_identical(colnames(.res), colnames(so)))
  )

# Check CDR3 length calculation
test_that("import_vdj CDR3 lengths", {
  res <- so@meta.data |>
    import_vdj(ctigs, include_mutations = FALSE)

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
  res <- so@meta.data |>
    import_vdj(ctigs, include_mutations = FALSE)

  dat <- res |>
    select(all_of(vdj_cols)) |>
    select(-clonotype_id) |>
    mutate(across(everything(), ~ lengths(gregexpr(";", .x))))

  dat |>
    purrr::reduce(expect_identical)
})

# Check include_mutations
# this checks for bug that duplicated read and umi counts with include_mutations
test_that("import_vdj include_mutations", {
  res_1 <- import_vdj(
    vdj_dir = ctigs,
    include_mutations = TRUE
  )

  res_1 <- res_1 |>
    dplyr::select(-matches("^(all|[vdjc]+)_(ins|del|mis)"))

  res_2 <- import_vdj(
    vdj_dir = ctigs,
    include_mutations = FALSE
  )

  expect_identical(res_1, res_2)
})

# Check include_constant
test_that("import_vdj include_mutations include_constant", {

  check_constant <- function(so, include_constant) {
    res <- so |>
      import_vdj(
        vdj_dir = ctigs,
        prefix  = "PREFIX_",
        include_mutations = TRUE,
        include_constant = include_constant
      )

    dat <- res@meta.data
    types <- c("ins", "del", "mis")

    all_types <- lapply(types, function(type) {
      dat_all <- dat[, grepl(paste0("PREFIX_all_", type, "$"), colnames(dat))]

      suppressWarnings(dat_all <- as.double(dat_all))

      dat_all[is.na(dat_all)] <- 0

      if(include_constant){
        dat_ind <- dat[ , grepl(paste0("PREFIX_[v|d|j|c]_", type, "$"), colnames(dat))]

      } else {
        dat_ind <- dat[ , grepl(paste0("PREFIX_[v|d|j]_", type, "$"), colnames(dat))]
      }

      suppressWarnings(dat_ind <- apply(dat_ind, 2, as.double))
      dat_ind[is.na(dat_ind)] <- 0

      return(data.frame(sum_val = rowSums(dat_ind), all_val = dat_all))
    })

    all_types <- do.call(rbind, all_types)
    return(all_types)
  }

  constant_included <- check_constant(so, include_constant = TRUE)
  constant_excluded <- check_constant(so, include_constant = FALSE)

  expect_identical(constant_included$sum_val, constant_included$all_val)
  expect_identical(constant_excluded$sum_val, constant_excluded$all_val)
})

# Check arguments for data.frame input
test_that("import_vdj data.frame", {
  res <- df_1 |>
    import_vdj(
      vdj_dir = ctigs,
      prefix  = "PREFIX_",
      include_mutations = FALSE
    )

  dat <- res[, !grepl("PREFIX_", colnames(res))]

  expect_identical(dat, df_1)
  expect_s3_class(res, "data.frame")
  expect_identical(rownames(res), rownames(df_1))

  res <- df_2 |>
    import_vdj(
      vdj_dir = ctigs,
      prefix  = "PREFIX_",
      include_mutations = FALSE
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
    res <- so |>
      import_vdj(vdj_dir = "BAD_PATH")
  }

  expect_error(fn(), "does not exist")
})

# Check column prefix
test_that("import_vdj column prefix", {
  prfx     <- "TEST_"
  new_cols <- paste0(prfx, vdj_cols)

  res <- so |>
    import_vdj(
      vdj_dir = ctigs,
      prefix  = prfx,
      include_mutations = FALSE
    )

  expect_false(any(vdj_cols %in% colnames(res@meta.data)))
  expect_true(all(new_cols %in% colnames(res@meta.data)))
})

# Check filtered contigs
test_that("import_vdj filter_chains", {
  res <- so |>
    import_vdj(ctigs, include_mutations = FALSE)

  expect_false(any(grepl("FALSE", res$productive)))
  expect_false(any(grepl("FALSE", res$full_length)))
  expect_s4_class(res, "Seurat")
  expect_identical(colnames(res), colnames(so))

  fn <- function() {
    so |>
      import_vdj(
        vdj_dir       = ctigs,
        filter_chains = FALSE,
        include_mutations = TRUE
      )
  }

  expect_warning(fn(), "When `include_mutations` is `TRUE`")
})

# Check filter_paired
test_that("import_vdj filter_paired", {
  res <- so |>
    import_vdj(
      vdj_dir       = ctigs,
      filter_paired = TRUE,
      include_mutations = FALSE
    )

  expect_true(all(res$n_chains > 1, na.rm = TRUE))
  expect_true(all(res$paired, na.rm = TRUE))
  expect_s4_class(res, "Seurat")
  expect_identical(colnames(res), colnames(so))

  res <- so |>
    import_vdj(
      vdj_dir       = ctigs,
      filter_paired = FALSE,
      include_mutations = FALSE
    )

  expect_true(any(!res$paired, na.rm = TRUE))
  expect_s4_class(res, "Seurat")
  expect_identical(colnames(res), colnames(so))
})

# Check define_clonotypes options
test_that("import_vdj define_clonotypes options", {
  opts <- c(cdr3 = "cdr3aa", cdr3_nt = "cdr3nt")

  opts |>
    iwalk(~ {
      res <- so |>
        import_vdj(
          vdj_dir = ctigs,
          define_clonotypes = .x,
          include_mutations = FALSE
        )

      expect_identical(n_distinct(res$clonotype_id), n_distinct(res[[.y]]))
      expect_s4_class(res, "Seurat")
      expect_identical(colnames(res), colnames(so))
    })

  res <- so |>
    import_vdj(
      vdj_dir = ctigs,
      define_clonotypes = "cdr3_gene",
      include_mutations = FALSE
    )

  x <- n_distinct(res$clonotype_id)
  y <- nrow(distinct(res@meta.data, cdr3_nt, v_gene, d_gene, j_gene))

  expect_identical(x, y)
  expect_s4_class(res, "Seurat")
  expect_identical(colnames(res), colnames(so))
})

# Check define_clonotypes bad define_clonotypes
test_that("import_vdj bad define_clonotypes", {
  expect_error(
    so |>
      import_vdj(
        vdj_dir = ctigs,
        define_clonotypes = "BAD",
        include_mutations = FALSE
      ),
    "`define_clonotypes` must be"
  )
})

# Check data.frame output
test_that("import_vdj df out", {
  res <- import_vdj(vdj_dir = ctigs, include_mutations = FALSE)

  expect_s3_class(res, "data.frame")
})

# Check bad barcode prefixes
test_that("import_vdj bad prefixes", {
  dat <- setNames(ctigs, c("A", "B"))

  fn <- function() {
    res <- so |>
      import_vdj(vdj_dir = dat)
  }

  expect_error(fn(), "do not match those")
})

# Check bad separator
test_that("import_vdj bad sep", {
  fn <- function() {
    res <- so |>
      import_vdj(
        vdj_dir = ctigs,
        sep     = "A",
        include_mutations = FALSE
      )
  }

  expect_error(fn(), "is already present in the input data", fixed = TRUE)
})

# Check BCR and TCR
test_that("import_vdj BCR and TCR", {
  fn <- function() {
    res <- so |>
      import_vdj(
        vdj_dir        = c(tcr_ctigs[1], ctigs[1]),
        include_mutations = FALSE
      )
  }

  expect_error(fn(), "Multiple data types detected")
})

# # SHOULD ALSO CHECK WHEN TOO MANY VDJ DIRS ARE PROVIDED
# test_that("import_vdj too may vdj_dirs provided", {
#   fn <- function() {
#     res <- so |>
#       import_vdj(
#         vdj_dir        = c(tcr_ctigs[1], ctigs[1], ctigs[2]),
#         include_mutations = FALSE
#       )
#   }
#
#   expect_error(fn(), "Multiple data types detected")
# })

# Check .classify_vdj
test_that(".classify_vdj", {
  dat <- vdj_sce |>
    fetch_vdj(data_cols = "chains")

  expect_error(
    dat |>
      mutate(chains = paste0("A", chains)) |>
      .classify_vdj(),
    "None of the expected chains"
  )

  expect_error(
    dat |>
      mutate(chains = paste0(chains, "A")) |>
      .classify_vdj(),
    "None of the expected chains"
  )

  expect_warning(
    dat |>
      filter(chains %in% c("IGH", "IGK")) |>
      group_by(chains) |>
      dplyr::slice(1:10) |>
      ungroup() |>
      mutate(chains = gsub("IGH", "TRA", chains)) |>
      .classify_vdj(),
    "Equal number of BCR.+and TCR.+chains detected"
  )
})

