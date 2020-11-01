
test_that("Check Seurat output", {
  res <- tiny_vdj %>%
    calc_diversity(
      clonotype_col = "clonotype_id",
      return_seurat = TRUE
    )

  expect_s4_class(res, "Seurat")                       # class
  expect_identical(colnames(res), colnames(tiny_vdj))  # cells in object
})

test_that("Check tibble output", {
  res <- tiny_vdj %>%
    calc_diversity(
      clonotype_col = "clonotype_id",
      cluster_col   = "seurat_clusters",
      return_seurat = FALSE
    )

  old_nms <- tiny_vdj@meta.data %>%
    filter(!is.na(clonotype_id)) %>%
    pull(seurat_clusters) %>%
    unique()

  expect_s3_class(res, "tbl")
  expect_identical(sort(res$seurat_clusters), sort(old_nms))
})

test_that("Check numeric output", {
  res <- tiny_vdj %>%
    calc_diversity(
      clonotype_col = "clonotype_id",
      return_seurat = FALSE,
      method        = abdiv::simpson
    )

  expect_type(res, "double")
  expect_identical(names(res), "simpson")
})

test_that("Check methods", {
  fns <- map(abdiv::alpha_diversities, ~ {
    eval(parse(text = paste0("abdiv::", .x)))
  })

  names(fns) <- abdiv::alpha_diversities

  res <- tiny_vdj %>%
    calc_diversity(
      clonotype_col = "clonotype_id",
      method = fns
    )

  expect_true(all(abdiv::alpha_diversities %in% names(res@meta.data)))
})
