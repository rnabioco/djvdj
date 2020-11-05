
# Check all calc_diversity arguments
mets <- abdiv::alpha_diversities %>%
  map(~ eval(parse(text = paste0("abdiv::", .x))))

lst <- list(
  sobj_in       = list(tiny_vdj),
  clonotype_col = "cdr3",
  cluster_col   = list(NULL, "seurat_clusters"),
  method        = mets,
  prefix        = c("", "X"),
  return_seurat = c(TRUE, FALSE)
)

test_all_args(
  lst    = lst,
  .fn    = calc_diversity,
  ttl    = "calc_diversity args",
  chk_fn = expect_silent
)

# Check diversity calculation
mets <- list(
  shannon = abdiv::shannon,
  simpson = abdiv::simpson
)

test_div <- tiny_vdj@meta.data %>%
  as_tibble(rownames = ".cell_id") %>%
  filter(!is.na(cdr3)) %>%
  group_by(cdr3, seurat_clusters) %>%
  summarize(n = n_distinct(.cell_id), .groups = "drop") %>%
  group_by(seurat_clusters) %>%
  summarize(
    shannon = abdiv::shannon(n),
    simpson = abdiv::simpson(n),
    .groups = "drop"
  ) %>%
  arrange(seurat_clusters)

test_that("div calc Seurat out", {
  res <- tiny_vdj %>%
    calc_diversity(
      clonotype_col = "cdr3",
      cluster_col   = "seurat_clusters",
      return_seurat = TRUE,
      method        = mets
    )

  res <- res@meta.data %>%
    as_tibble() %>%
    filter(!is.na(cdr3)) %>%
    select(seurat_clusters, all_of(names(mets))) %>%
    distinct() %>%
    arrange(seurat_clusters)

  expect_identical(res, test_div)
})

test_that("div calc tbl out", {
  res <- tiny_vdj %>%
    calc_diversity(
      clonotype_col = "cdr3",
      cluster_col   = "seurat_clusters",
      return_seurat = FALSE,
      method        = mets
    ) %>%
    arrange(seurat_clusters)

  expect_s3_class(res, "tbl")
  expect_identical(res, test_div)
})

# Check Seurat output
test_that("calc_diversity Seurat out", {
  fns <- abdiv::alpha_diversities %>%
    map(~ eval(parse(text = paste0("abdiv::", .x))))

  names(fns) <- abdiv::alpha_diversities

  res <- tiny_vdj %>%
    calc_diversity(
      clonotype_col = "cdr3",
      method = fns
    )

  res@meta.data <- res@meta.data %>%
    select(-all_of(names(fns)))

  expect_identical(res, tiny_vdj)
})

# Check single method no cluster_col
test_that("single met no cluster_col", {
  res <- tiny_vdj %>%
    calc_diversity(
      clonotype_col = "cdr3",
      return_seurat = FALSE,
      method        = abdiv::simpson
    )

  expect_type(res, "double")
  expect_identical(names(res), "simpson")
})

