# Check NAs in result
check_nas <- function(df_in) {

  vdj_cols <- df_in |>
    as.data.frame() |>
    .get_vdj_cols(
      clone_col = "clonotype_id",
      cols_in   = NULL,
      sep       = ";"
    )

  df_in <- df_in |>
    as_tibble() |>
    select(all_of(vdj_cols$vdj)) |>
    mutate(across(all_of(vdj_cols$vdj), is.na))

  c_nas <- df_in$clonotype_id

  expect_true(all(map_lgl(df_in, ~ identical(.x, c_nas))))
}

# Check with no filtering
test_that("filter_vdj no filtering", {
  res <- vdj_sce |>
    filter_vdj(any(chains %in% c("IGH", "IGK", "IGL")))

  expect_identical(res, vdj_sce)

  res <- vdj_sce |>
    filter_vdj(chains %in% c("IGH", "IGK", "IGL"))

  expect_identical(res, vdj_sce)

  res <- vdj_sce |>
    filter_vdj(nCount_RNA < Inf)

  expect_identical(res, vdj_sce)
})

# Check per-chain filtering
test_that("filter_vdj per-chain", {
  res <- vdj_sce |>
    filter_vdj(umis > 10) |>
    fetch_vdj("umis")

  clmns <- colnames(res)
  clmns <- clmns[clmns != ".cell_id"]

  expect_false(any(res$umis <= 10, na.rm = TRUE))
  expect_identical(clmns, colnames(vdj_sce@colData))

  res <- vdj_sce |>
    filter_vdj(chains == "IGH")

  check_nas(res@colData)

  expect_false(any(grepl("IGK", res$chains)))
  expect_false(any(grepl("IGL", res$chains)))
  expect_identical(colnames(res@colData), colnames(vdj_sce@colData))
})

# Check NULL sep
test_that("filter_vdj NULL sep", {
  res <- vdj_sce |>
    filter_vdj(chains == "IGH;IGK", sep = NULL)

  check_nas(res@colData)

  expect_true(all(res$chains == "IGH;IGK", na.rm = TRUE))
  expect_identical(colnames(res@colData), colnames(vdj_sce@colData))
})

# Check per-cell filtering
test_that("filter_vdj per-cell", {
  res <- vdj_sce@colData |>
    filter_vdj(nCount_RNA > 2000)

  check_nas(res)
  expect_true(nrow(filter(res, nCount_RNA <= 2000 & !is.na(clonotype_id))) == 0)
  expect_identical(colnames(res), colnames(vdj_sce@colData))
})

# Check for .KEEP
test_that("filter_vdj .KEEP check", {
  res <- vdj_sce |>
    filter_vdj(orig.ident == "avid_1")

  expect_false(".KEEP" %in% colnames(res@colData))

  fn <- function() {
    vdj_sce |>
      filter_vdj(
        orig.ident == "avid_1",
        sep = "BED_SEP"
      )
  }

  res <- suppressWarnings(fn())

  expect_warning(fn())
  expect_false(".KEEP" %in% colnames(res@colData))
})

# Check all cells filtered
test_that("filter_vdj all cells filtered", {
  res <- vdj_sce |>
    filter_vdj(chains == "BAD")

  vdj_cols <- res@colData |>
    as.data.frame() |>
    .get_vdj_cols(
      clone_col = "clonotype_id",
      cols_in   = NULL,
      sep       = ";"
    )

  res <- res@colData |>
    as.data.frame() |>
    filter(if_all(all_of(vdj_cols$vdj), is.na))

  expect_identical(nrow(res), ncol(vdj_sce))
})

# Check bad lengths
test_that("filter_vdj bad length", {
  expect_error(
    vdj_sce |>
      filter_vdj(c("IGH", "IGK") %in% chains),
    "Filtering condition must return `TRUE`/`FALSE` for each chain"
  )
})





