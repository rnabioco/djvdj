# Test data
df_1 <- vdj_so@meta.data

df_2 <- vdj_so@meta.data |>
  as_tibble(rownames = ".cell_id")

test_clsts <- df_1 |>
  na.omit() |>
  pull(seurat_clusters) |>
  unique()

# Check all calc_similarity arguments
mets <- abdiv::beta_diversities |>
  purrr::set_names() |>
  map(~ eval(parse(text = paste0("abdiv::", .x))))

arg_lst <- list(
  input       = list(vdj_so, vdj_sce, df_1, df_2),
  data_col    = "cdr3",
  cluster_col = "seurat_clusters",
  method      = mets,
  return_mat  = c(TRUE, FALSE),
  prefix      = "TEST_"
)

test_all_args(
  arg_lst = arg_lst,
  .fn     = calc_similarity,
  desc    = "calc_similarity args",
  chk     = expect_silent
)

# Check similarity calculation
# calculate similarity independently and compare to calc_similarity results
test_sim <- vdj_so@meta.data |>
  as_tibble(rownames = ".cell_id") |>
  filter(!is.na(cdr3)) |>
  group_by(cdr3, seurat_clusters) |>
  summarize(n = n_distinct(.cell_id), .groups = "drop") |>
  tidyr::pivot_wider(names_from = "seurat_clusters", values_from = "n") |>
  mutate(across(all_of(test_clsts), tidyr::replace_na, 0))

clst <- test_clsts[1]
nm   <- paste0("x", clst)

clst_dat <- pull(test_sim, clst)

get_sim_res <- function(method) {
  map_dfr(test_clsts, ~ {
    x <- pull(test_sim, clst)
    y <- pull(test_sim, .x)

    tibble(
      seurat_clusters = .x,
      !!sym(nm) := method(x, y)
    )
  }) |>
    arrange(seurat_clusters) |>
    mutate(
      !!sym(nm) := if_else(
        seurat_clusters == clst,
        method(clst_dat, clst_dat),
        !!sym(nm)
      )
    )
}

purrr::walk(abdiv::beta_diversities, ~ {
  fn <- paste0("abdiv::", .x)

  test_res <- get_sim_res(eval(parse(text = fn)))

  res <- vdj_so |>
    calc_similarity(
      data_col    = "cdr3",
      cluster_col = "seurat_clusters",
      method      = eval(parse(text = fn)),
      return_mat  = FALSE,
      prefix      = "x"
    )

  res <- res@meta.data |>
    as_tibble() |>
    select(seurat_clusters, all_of(nm)) |>
    filter(!is.na(seurat_clusters)) |>
    distinct() |>
    arrange(seurat_clusters)

  test_that(paste0("sim calc ", .x), {
    expect_identical(res, test_res)
  })
})

# Check Seurat output
test_that("calc_similarity Seurat out", {
  res <- vdj_so |>
    calc_similarity(
      data_col    = "cdr3",
      cluster_col = "seurat_clusters",
      method      = abdiv::binomial_deviance,
      return_mat  = FALSE,
      prefix      = "x"
    )

  expect_s4_class(res, "Seurat")

  res@meta.data <- res@meta.data |>
    select(-all_of(paste0("x", test_clsts)))

  expect_identical(res, vdj_so)
})

# Check SCE output
test_that("calc_similarity Seurat out", {
  res <- vdj_sce |>
    calc_similarity(
      data_col    = "cdr3",
      cluster_col = "seurat_clusters",
      method      = abdiv::binomial_deviance,
      return_mat  = FALSE,
      prefix      = "x"
    )

  expect_s4_class(res, "SingleCellExperiment")

  new_clmns <- colnames(res@colData)
  new_clmns <- new_clmns[!new_clmns %in% paste0("x", test_clsts)]

  res@colData <- res@colData[, new_clmns]

  expect_identical(res, vdj_sce)
})

# Check data.frame output
test_that("calc_similarity df out", {
  res <- vdj_so@meta.data |>
    calc_similarity(
      data_col    = "cdr3",
      cluster_col = "seurat_clusters",
      method      = abdiv::binomial_deviance,
      return_mat  = FALSE,
      prefix      = "x"
    )

  expect_s3_class(res, "data.frame")

  res_2 <- res |>
    select(-all_of(paste0("x", test_clsts)))

  expect_identical(res_2, vdj_so@meta.data)
})

# Check matrix output
test_that("calc_similarity mat out", {
  res <- vdj_so |>
    calc_similarity(
      data_col    = "cdr3",
      cluster_col = "seurat_clusters",
      method      = abdiv::binomial_deviance,
      return_mat  = TRUE,
      prefix      = "x"
    )

  expect_type(res, "double")

  res <- unique(c(rownames(res), colnames(res)))

  expect_true(all(sort(test_clsts) == sort(res)))
})

