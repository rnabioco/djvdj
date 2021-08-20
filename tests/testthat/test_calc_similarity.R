
# Test inputs
test_clsts <- tiny_vdj@meta.data %>%
  na.omit() %>%
  pull(seurat_clusters) %>%
  unique()

# Check all calc_similarity arguments
mets <- abdiv::beta_diversities %>%
  map(~ eval(parse(text = paste0("abdiv::", .x))))

arg_lst <- list(
  input         = list(tiny_vdj, tiny_vdj@meta.data),
  clonotype_col = "cdr3",
  cluster_col   = "seurat_clusters",
  method        = mets,
  return_matrix = c(TRUE, FALSE)
)

test_all_args(
  arg_lst = arg_lst,
  .fn     = calc_similarity,
  desc    = "calc_similarity args",
  chk     = expect_silent
)

# Check similarity calculation
test_sim <- tiny_vdj@meta.data %>%
  as_tibble(rownames = ".cell_id") %>%
  filter(!is.na(cdr3)) %>%
  group_by(cdr3, seurat_clusters) %>%
  summarize(n = n_distinct(.cell_id), .groups = "drop") %>%
  pivot_wider(names_from = "seurat_clusters", values_from = "n") %>%
  mutate(across(all_of(test_clsts), replace_na, 0))

clst <- test_clsts[1]
nm   <- paste0("x", clst)

test_sim <- map_dfr(test_clsts, ~ {
  x <- pull(test_sim, clst)
  y <- pull(test_sim, .x)

  tibble(
    seurat_clusters = .x,
    !!sym(nm) := abdiv::binomial_deviance(x, y)
  )
}) %>%
  arrange(seurat_clusters) %>%
  mutate(!!sym(nm) := if_else(seurat_clusters == clst, 1, !!sym(nm)))

test_that("sim calc", {
  res <- tiny_vdj %>%
    calc_similarity(
      clonotype_col = "cdr3",
      cluster_col   = "seurat_clusters",
      method        = abdiv::binomial_deviance,
      return_matrix = FALSE,
      prefix        = "x"
    )

  res <- res@meta.data %>%
    as_tibble() %>%
    select(seurat_clusters, all_of(nm)) %>%
    filter(!is.na(!!sym(nm))) %>%
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
      return_matrix = FALSE,
      prefix        = "x"
    )

  expect_s4_class(res, "Seurat")

  res@meta.data <- res@meta.data %>%
    select(-all_of(paste0("x", test_clsts)))

  expect_identical(res, tiny_vdj)
})

# Check data.frame output
test_that("calc_similarity df out", {
  res <- tiny_vdj@meta.data %>%
    calc_similarity(
      clonotype_col = "cdr3",
      cluster_col   = "seurat_clusters",
      method        = abdiv::binomial_deviance,
      return_matrix = FALSE,
      prefix        = "x"
    )

  expect_s3_class(res, "data.frame")

  res_2 <- res %>%
    select(-all_of(paste0("x", test_clsts)))

  expect_identical(res_2, tiny_vdj@meta.data)

  res <- res %>%
    as_tibble() %>%
    select(seurat_clusters, all_of(nm)) %>%
    filter(!is.na(!!sym(nm))) %>%
    distinct() %>%
    arrange(seurat_clusters)

  expect_identical(res, test_sim)
})

# Check data.frame output
test_that("calc_similarity df out", {
  res <- tiny_vdj@meta.data %>%
    calc_similarity(
      clonotype_col = "cdr3",
      cluster_col   = "seurat_clusters",
      method        = abdiv::binomial_deviance,
      return_matrix = FALSE,
      prefix        = "x"
    )

  expect_s3_class(res, "data.frame")

  res_2 <- res %>%
    select(-all_of(paste0("x", test_clsts)))

  expect_identical(res_2, tiny_vdj@meta.data)

  res <- res %>%
    as_tibble() %>%
    select(seurat_clusters, all_of(nm)) %>%
    filter(!is.na(!!sym(nm))) %>%
    distinct() %>%
    arrange(seurat_clusters)

  expect_identical(res, test_sim)
})

# Check matrix output
test_that("calc_similarity mat out", {
  res <- tiny_vdj %>%
    calc_similarity(
      clonotype_col = "cdr3",
      cluster_col   = "seurat_clusters",
      method        = abdiv::binomial_deviance,
      return_matrix = TRUE,
      prefix        = "x"
    )

  expect_type(res, "double")

  res <- unique(c(rownames(res), colnames(res)))

  expect_true(all(sort(test_clsts) == sort(res)))
})

