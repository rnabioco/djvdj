
# Test inputs
test_clsts <- tiny_vdj@meta.data %>%
  na.omit() %>%
  pull(seurat_clusters) %>%
  unique()

# Check all calc_similarity arguments
mets <- abdiv::beta_diversities %>%
  map(~ eval(parse(text = paste0("abdiv::", .x))))

lst <- list(
  sobj_in       = list(tiny_vdj),
  clonotype_col = "cdr3",
  cluster_col   = "seurat_clusters",
  method        = mets,
  return_seurat = c(TRUE, FALSE)
)

test_all_args(
  lst    = lst,
  .fn    = calc_similarity,
  ttl    = "calc_similarity args",
  chk_fn = expect_silent
)

# Check similarity calculation
test_sim <- tiny_vdj@meta.data %>%
  as_tibble(rownames = ".cell_id") %>%
  filter(!is.na(cdr3)) %>%
  group_by(cdr3, seurat_clusters) %>%
  summarize(n = n_distinct(.cell_id), .groups = "drop") %>%
  pivot_wider(names_from = "seurat_clusters", values_from = "n") %>%
  mutate(across(all_of(test_clsts), replace_na, 0))

test_sim <- map_dfr(test_clsts, ~ {
  x <- test_sim$`6`
  y <- pull(test_sim, .x)

  tibble(
    seurat_clusters = .x,
    x6 = abdiv::binomial_deviance(x, y)
  )
}) %>%
  arrange(seurat_clusters) %>%
  mutate(x6 = if_else(seurat_clusters == "6", 1, x6))

test_that("sim calc", {
  res <- tiny_vdj %>%
    calc_similarity(
      clonotype_col = "cdr3",
      cluster_col   = "seurat_clusters",
      method        = abdiv::binomial_deviance,
      return_seurat = TRUE,
      prefix        = "x"
    )

  res <- res@meta.data %>%
    as_tibble() %>%
    select(seurat_clusters, x6) %>%
    filter(!is.na(x6)) %>%
    distinct() %>%
    arrange(seurat_clusters)

  expect_identical(res, test_sim)
})

# Check Seurat output
test_that("calc_similarity Seurat out", {
  res <- tiny_vdj %>%
    calc_similarity(
      clonotype_col = "cdr3",
      cluster_col   = "seurat_clusters",
      method        = abdiv::binomial_deviance,
      return_seurat = TRUE,
      prefix        = "x"
    )

  expect_s4_class(res, "Seurat")

  res@meta.data <- res@meta.data %>%
    select(-all_of(paste0("x", test_clsts)))

  expect_identical(res, tiny_vdj)
})

# Check matrix output
test_that("calc_similarity mat out", {
  res <- tiny_vdj %>%
    calc_similarity(
      clonotype_col = "cdr3",
      cluster_col   = "seurat_clusters",
      method        = abdiv::binomial_deviance,
      return_seurat = FALSE,
      prefix        = "x"
    )

  expect_type(res, "double")

  res <- unique(c(rownames(res), colnames(res)))

  expect_true(all(sort(test_clsts) == sort(res)))
})
