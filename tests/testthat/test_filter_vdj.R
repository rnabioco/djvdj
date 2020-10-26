
test_that("Check Seurat return without filtering", {
  res <- filter_vdj(
    sobj_in = test_vdj,
    filt    = any(chains %in% c("IGH", "IGK", "IGL")),
  )

  expect_true("Seurat" %in% class(res))                # class
  expect_identical(colnames(res), colnames(test_vdj))  # cells in object
  expect_identical(res@meta.data, test_vdj@meta.data)  # meta.data

  new_na <- colnames(res)[is.na(res$clonotype_id)]
  old_na <- colnames(test_vdj)[is.na(test_vdj$clonotype_id)]
  expect_identical(old_na, new_na)                     # cells w/o VDJ data
})

# test_that("Check Seurat return with new_col", {
#   res <- filter_vdj(
#     sobj_in = test_vdj,
#     filt    = any(chains %in% c("IGH", "IGK", "IGL")),
#     new_col = "NEW_COL"
#   )
#
#   expect_true("Seurat" %in% class(res))                # class
#   expect_identical(colnames(res), colnames(test_vdj))  # cells in object
#   expect_true("NEW_COL" %in% colnames(res@meta.data))  # meta.data
#
#   new_na <- colnames(res)[is.na(res$clonotype_id)]
#   old_na <- colnames(test_vdj)[is.na(test_vdj$clonotype_id)]
#   expect_identical(old_na, new_na)                     # cells w/o VDJ data
# })

test_that("Check data.frame return without filtering", {
  res <- filter_vdj(
    sobj_in       = test_vdj,
    filt          = any(chains %in% c("IGH", "IGK", "IGL"))
  )

  res <- res@meta.data

  expect_identical(rownames(res), colnames(test_vdj))                        # cells in object
  expect_identical(sort(colnames(res)), sort(colnames(test_vdj@meta.data)))  # meta.data columns are not in the same order

  new_na <- rownames(res)[is.na(res$clonotype_id)]
  old_na <- colnames(test_vdj)[is.na(test_vdj$clonotype_id)]
  expect_identical(old_na, new_na)                                           # cells w/o VDJ data
})

test_that("Check Seurat return for only cells with VDJ data", {
  res <- filter_vdj(
    sobj_in       = test_vdj,
    filt          = any(chains %in% c("IGH", "IGK", "IGL")),
    clonotype_col = NULL,
    filter_cells  = TRUE
  )

  old_na <- colnames(test_vdj)[!is.na(test_vdj$clonotype_id)]

  expect_true("Seurat" %in% class(res))    # class
  expect_identical(colnames(res), old_na)  # cells w/ VDJ data
})


