
# Check Seurat return without filtering
test_that("filter_vdj Seurat out no filt", {
  res <- filter_vdj(
    sobj_in = tiny_vdj,
    filt    = any(chains %in% c("IGH", "IGK", "IGL")),
  )

  expect_s4_class(res, "Seurat")
  expect_identical(colnames(res), colnames(tiny_vdj))                          # cells in object
  expect_identical(res@meta.data, tiny_vdj@meta.data)                          # meta.data

  new_na <- colnames(res)[is.na(res$clonotype_id)]
  old_na <- colnames(tiny_vdj)[is.na(tiny_vdj$clonotype_id)]
  expect_identical(old_na, new_na)                                             # cells w/o VDJ data
})

# Check data.frame return without filtering
test_that("filter_vdj df out no filt", {
  res <- filter_vdj(
    sobj_in = tiny_vdj,
    filt    = any(chains %in% c("IGH", "IGK", "IGL"))
  )

  res <- res@meta.data

  expect_identical(rownames(res), colnames(tiny_vdj))                          # cells in object
  expect_identical(sort(colnames(res)), sort(colnames(tiny_vdj@meta.data)))    # meta.data columns are not in the same order

  new_na <- rownames(res)[is.na(res$clonotype_id)]
  old_na <- colnames(tiny_vdj)[is.na(tiny_vdj$clonotype_id)]
  expect_identical(old_na, new_na)                                             # cells w/o VDJ data
})

# Check Seurat return for all cells
test_that("filter_vdj Seurat return VDJ cells", {
  res <- filter_vdj(
    sobj_in       = tiny_vdj,
    filt          = any(chains %in% c("IGH", "IGK", "IGL")),
    clonotype_col = NULL,
    filter_cells  = TRUE
  )

  old_na <- colnames(tiny_vdj)[!is.na(tiny_vdj$clonotype_id)]

  expect_s4_class(res, "Seurat")
  expect_identical(colnames(res), old_na)                                      # cells w/ VDJ
})

# Check Seurat return for only VDJ cells
test_that("filter_vdj Seurat return all cells", {
  res <- tiny_vdj %>%
    filter_vdj(
      filt          = any(chains %in% c("IGH", "IGK", "IGL")),
      clonotype_col = "cdr3_nt",
      filter_cells  = TRUE
    )

  expect_s4_class(res, "Seurat")
  expect_identical(colnames(res), colnames(tiny_vdj))
})

# Check NULL clonotype_col and FALSE filter_cells
test_that("filter_vdj NULL clonotype_col FALSE filter_cells", {
  fn <- function() {
    res <- tiny_vdj %>%
      filter_vdj(
        filt          = any(chains %in% c("IGH", "IGK", "IGL")),
        clonotype_col = NULL,
        filter_cells  = FALSE
      )
  }

  expect_error(fn())
})

# Check for .KEEP
test_that("filter_vdj .KEEP check", {
  res <- tiny_vdj %>%
    filter_vdj(filt = orig.ident == "avid_1")

  expect_false(".KEEP" %in% colnames(res@meta.data))

  fn <- function() {
    tiny_vdj %>%
      filter_vdj(
        filt = orig.ident == "avid_1",
        sep  = "BED_SEP"
      )
  }

  res <- suppressWarnings(fn())

  expect_warning(fn())
  expect_false(".KEEP" %in% colnames(res@meta.data))
})





