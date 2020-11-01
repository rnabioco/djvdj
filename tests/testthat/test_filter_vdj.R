
test_that("Seurat return without filtering", {
  res <- filter_vdj(
    sobj_in = tiny_vdj,
    filt    = any(chains %in% c("IGH", "IGK", "IGL")),
  )

  expect_s4_class(res, "Seurat")                       # class
  expect_identical(colnames(res), colnames(tiny_vdj))  # cells in object
  expect_identical(res@meta.data, tiny_vdj@meta.data)  # meta.data

  new_na <- colnames(res)[is.na(res$clonotype_id)]
  old_na <- colnames(tiny_vdj)[is.na(tiny_vdj$clonotype_id)]
  expect_identical(old_na, new_na)                     # cells w/o VDJ data
})

test_that("data.frame return without filtering", {
  res <- filter_vdj(
    sobj_in       = tiny_vdj,
    filt          = any(chains %in% c("IGH", "IGK", "IGL"))
  )

  res <- res@meta.data

  expect_identical(rownames(res), colnames(tiny_vdj))                        # cells in object
  expect_identical(sort(colnames(res)), sort(colnames(tiny_vdj@meta.data)))  # meta.data columns are not in the same order

  new_na <- rownames(res)[is.na(res$clonotype_id)]
  old_na <- colnames(tiny_vdj)[is.na(tiny_vdj$clonotype_id)]
  expect_identical(old_na, new_na)                                           # cells w/o VDJ data
})

test_that("Seurat return for only cells with VDJ data", {
  res <- filter_vdj(
    sobj_in       = tiny_vdj,
    filt          = any(chains %in% c("IGH", "IGK", "IGL")),
    clonotype_col = NULL,
    filter_cells  = TRUE
  )

  old_na <- colnames(tiny_vdj)[!is.na(tiny_vdj$clonotype_id)]

  expect_s4_class(res, "Seurat")           # class
  expect_identical(colnames(res), old_na)  # cells w/ VDJ data
})


