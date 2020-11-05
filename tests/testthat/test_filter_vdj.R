
# Check Seurat return without filtering
test_that("filter_vdj Seurat out no filt", {
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

# Check data.frame return without filtering
test_that("filter_vdj df out no filt", {
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

# Check Seurat return for only VDJ cells
test_that("filter_vdj Seurat return VDJ cells", {
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


