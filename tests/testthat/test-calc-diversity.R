# Test data
df_1 <- vdj_so@meta.data

df_2 <- vdj_so@meta.data %>%
  as_tibble(rownames = ".cell_id")

# Check all calc_diversity arguments
mets <- abdiv::alpha_diversities %>%
  map(~ eval(parse(text = paste0("abdiv::", .x))))

names(mets) <- abdiv::alpha_diversities

arg_lst <- list(
  input         = list(vdj_so, vdj_sce, df_1, df_2),
  clonotype_col = "cdr3",
  cluster_col   = list(NULL, "seurat_clusters"),
  method        = append(append(abdiv::shannon, mets), list(mets)),
  prefix        = c("", "X"),
  return_df     = c(TRUE, FALSE)
)

test_all_args(
  arg_lst = arg_lst,
  .fn     = calc_diversity,
  desc    = "calc_diversity args",
  chk     = expect_silent
)

# Check diversity calculation
mets <- list(
  shannon = abdiv::shannon,
  simpson = abdiv::simpson
)

test_div <- vdj_so@meta.data %>%
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
  arrange(seurat_clusters) %>%
  as.data.frame()

test_that("div calc Seurat out", {
  res <- vdj_so %>%
    calc_diversity(
      clonotype_col = "cdr3",
      cluster_col   = "seurat_clusters",
      return_df     = FALSE,
      method        = mets
    )

  res <- res@meta.data %>%
    as_tibble() %>%
    filter(!is.na(cdr3)) %>%
    select(seurat_clusters, all_of(names(mets))) %>%
    distinct() %>%
    arrange(seurat_clusters) %>%
    as.data.frame()

  expect_identical(res, test_div)
})

test_that("div calc df out", {
  res <- vdj_so %>%
    calc_diversity(
      clonotype_col = "cdr3",
      cluster_col   = "seurat_clusters",
      return_df     = TRUE,
      method        = mets
    ) %>%
    arrange(seurat_clusters) %>%
    distinct(seurat_clusters, !!!syms(names(mets))) %>%
    filter(across(all_of(names(mets)), ~ !is.na(.x))) %>%
    remove_rownames()

  expect_s3_class(res, "data.frame")
  expect_identical(res, test_div)
})

# Check Seurat output
test_that("calc_diversity Seurat out", {
  fns <- abdiv::alpha_diversities %>%
    map(~ eval(parse(text = paste0("abdiv::", .x))))

  names(fns) <- abdiv::alpha_diversities

  res <- vdj_so %>%
    calc_diversity(
      clonotype_col = "cdr3",
      method        = fns
    )

  res@meta.data <- res@meta.data %>%
    select(-all_of(names(fns)))

  expect_identical(res, vdj_so)
})

# Check data.frame input
test_that("calc_diversity df in", {
  res <- vdj_so@meta.data %>%
    calc_diversity(
      clonotype_col = "cdr3",
      cluster_col   = "seurat_clusters",
      method        = mets
    ) %>%
    arrange(seurat_clusters) %>%
    distinct(seurat_clusters, !!!syms(names(mets))) %>%
    filter(across(all_of(names(mets)), ~ !is.na(.x))) %>%
    remove_rownames()

  expect_s3_class(res, "data.frame")
  expect_identical(res, test_div)
})

# Check bad method list
test_that("calc_diversity bad method list", {
  expect_error(
    vdj_so %>% calc_diversity(method = list(abdiv::simpson, abdiv::mcintosh_d)),
    "Must include names if using a list of methods"
  )

  res <- vdj_so %>%
    calc_diversity(method = abdiv::simpson)

  expect_true("simpson" %in% colnames(res@meta.data))
})
