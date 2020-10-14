
# Test Seurat return without filtering
test_that("filter_vdj return Seurat", {
  res <- filter_vdj(
    sobj_in       = test_vdj,
    filt          = orig.ident == orig.ident,
  )

  expect_true("Seurat" %in% class(res))                # class
  expect_identical(colnames(res), colnames(test_vdj))  # cells in object
  expect_identical(res@meta.data, test_vdj@meta.data)  # meta.data
})

# Test Seurat return with new_col
test_that("filter_vdj return Seurat", {
  res <- filter_vdj(
    sobj_in       = test_vdj,
    filt          = orig.ident == orig.ident,
    new_col       = "NEW_COL"
  )

  expect_true("Seurat" %in% class(res))                # class
  expect_identical(colnames(res), colnames(test_vdj))  # cells in object
  expect_true("NEW_COL" %in% colnames(res@meta.data))  # meta.data
})

# Test data.frame return without filtering
test_that("filter_vdj return data.frame", {
  res <- filter_vdj(
    sobj_in       = test_vdj,
    filt          = orig.ident == orig.ident,
    return_seurat = FALSE
  )

  expect_true("data.frame" %in% class(res))
  expect_identical(rownames(res), colnames(test_vdj))
  expect_identical(sort(colnames(res)), sort(colnames(test_vdj@meta.data)))  # meta.data columns are not in the same order

  new_na <- rownames(res)[is.na(res$clonotype_id)]
  old_na <- colnames(test_vdj)[is.na(test_vdj@meta.data$clonotype_id)]
  expect_identical(old_na, new_na)
})










