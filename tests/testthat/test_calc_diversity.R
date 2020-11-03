
test_that("Check all calc_diversity arguments", {
  mets <- abdiv::alpha_diversities %>%
    map(~ eval(parse(text = paste0("abdiv::", .x))))

  args_df <- list(
    cluster_col   = list(NULL, "seurat_clusters"),
    method        = mets,
    prefix        = c("", "X"),
    return_seurat = c(TRUE, FALSE)
  ) %>%
    expand.grid()

  res <- pmap(
    args_df,
    check_args,
    .fn           = calc_diversity,
    sobj_in       = tiny_vdj,
    clonotype_col = "cdr3"
  )
})

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
