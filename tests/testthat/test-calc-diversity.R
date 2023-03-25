# Test data
df_1 <- vdj_so@meta.data

df_2 <- vdj_so@meta.data |>
  as_tibble(rownames = ".cell_id")

# Check all calc_diversity arguments
mets <- abdiv::alpha_diversities |>
  map(~ eval(parse(text = paste0("abdiv::", .x))))

names(mets) <- abdiv::alpha_diversities

arg_lst <- list(
  input       = list(vdj_so),
  data_col    = "cdr3",
  chain       = list(NULL, "IGK"),
  cluster_col = list(NULL, "seurat_clusters"),
  downsample  = c(TRUE, FALSE),
  method      = list(mets, list(simpson = abdiv::simpson)),
  return_df   = FALSE,
  n_boots     = 0
)

test_all_args(
  arg_lst = arg_lst,
  .fn     = calc_diversity,
  desc    = "calc_diversity Seurat args",
  chk     = expr(expect_s4_class(.res, "Seurat"))
)

arg_lst$input <- list(vdj_sce)

test_all_args(
  arg_lst = arg_lst,
  .fn     = calc_diversity,
  desc    = "calc_diversity SCE args",
  chk     = expr(expect_s4_class(.res, "SingleCellExperiment"))
)

arg_lst$input <- list(df_1)

test_all_args(
  arg_lst = arg_lst,
  .fn     = calc_diversity,
  desc    = "calc_diversity data.frame args",
  chk     = expr(expect_true(is.data.frame(.res)))
)

arg_lst$input <- list(vdj_so, vdj_sce, df_1, df_2)
arg_lst$return_df <- TRUE

test_all_args(
  arg_lst = arg_lst,
  .fn     = calc_diversity,
  desc    = "calc_diversity return_df",
  chk     = expr(expect_true(is.data.frame(.res)))
)

# Test n_boots
arg_lst$input     <- list(test_so)
arg_lst$n_boots   <- 10
arg_lst$return_df <- FALSE

set.seed(42)

test_all_args(
  arg_lst = arg_lst,
  .fn     = calc_diversity,
  desc    = "calc_diversity n_boots",
  chk     = expr(expect_s4_class(.res, "Seurat"))
)

# Check diversity column prefix
test_that("calc_diversity prefix", {
  res <- vdj_so |>
    calc_diversity(
      data_col    = "cdr3",
      cluster_col = "seurat_clusters",
      return_df   = FALSE,
      method      = mets,
      prefix      = "X_"
    )

  nms <- paste0("X_", names(mets), "_diversity")

  expect_true(all(nms %in% names(res@meta.data)))

  res <- res |>
    mutate_meta(select, -starts_with("X_"))

  expect_identical(res, vdj_so)
})

# Check diversity calculation
# only check the clusters that have the fewest cells since each cluster will be
# downsampled to same number of cells
mets <- list(
  shannon = abdiv::shannon,
  simpson = abdiv::simpson
)

nms <- paste0("cdr3_", names(mets), "_diversity")

test_div <- vdj_so@meta.data |>
  as_tibble(rownames = ".cell_id") |>
  filter(!is.na(cdr3))

sm_clst <- table(test_div$seurat_clusters)
sm_clst <- names(sm_clst[sm_clst == min(sm_clst)])

test_div <- test_div |>
  filter(seurat_clusters %in% sm_clst) |>
  group_by(cdr3, seurat_clusters) |>
  summarize(n = n_distinct(.cell_id), .groups = "drop") |>
  group_by(seurat_clusters) |>
  summarize(
    cdr3_shannon_diversity = abdiv::shannon(n),
    cdr3_simpson_diversity = abdiv::simpson(n),
    .groups = "drop"
  ) |>
  arrange(seurat_clusters) |>
  as.data.frame()


test_that("div calc Seurat out", {
  set.seed(42)

  res <- vdj_so |>
    calc_diversity(
      data_col    = "cdr3",
      cluster_col = "seurat_clusters",
      return_df   = FALSE,
      method      = mets
    )

  res <- res@meta.data |>
    as_tibble() |>
    filter(!is.na(cdr3)) |>
    select(seurat_clusters, all_of(nms)) |>
    distinct() |>
    arrange(seurat_clusters) |>
    as.data.frame() |>
    filter(seurat_clusters %in% sm_clst)

  expect_identical(res, test_div)
})

test_that("div calc df out", {
  res <- vdj_so |>
    calc_diversity(
      data_col = "cdr3",
      cluster_col   = "seurat_clusters",
      return_df     = TRUE,
      method        = mets
    ) |>
    arrange(seurat_clusters) |>
    distinct(seurat_clusters, !!!syms(nms)) |>
    filter(if_all(all_of(nms), ~ !is.na(.x))) |>
    tibble::remove_rownames() |>
    filter(seurat_clusters %in% sm_clst)

  expect_s3_class(res, "data.frame")
  expect_identical(res, test_div)
})

# Check Seurat output
test_that("calc_diversity Seurat out", {
  fns <- abdiv::alpha_diversities |>
    map(~ eval(parse(text = paste0("abdiv::", .x))))

  names(fns) <- abdiv::alpha_diversities

  res <- vdj_so |>
    calc_diversity(
      data_col = "cdr3",
      method   = fns
    )

  res@meta.data <- res@meta.data |>
    select(-all_of(paste0("cdr3_", names(fns), "_diversity")))

  expect_identical(res, vdj_so)
})

# Check data.frame input
test_that("calc_diversity df in", {
  res <- vdj_so@meta.data |>
    calc_diversity(
      data_col    = "cdr3",
      cluster_col = "seurat_clusters",
      method      = mets
    ) |>
    arrange(seurat_clusters) |>
    distinct(seurat_clusters, !!!syms(nms)) |>
    filter(if_all(all_of(nms), ~ !is.na(.x))) |>
    tibble::remove_rownames() |>
    filter(seurat_clusters %in% sm_clst)

  expect_s3_class(res, "data.frame")
  expect_identical(res, test_div)
})

# Check bad method list
test_that("calc_diversity bad method list", {
  expect_error(
    vdj_so |>
      calc_diversity(
        data_col = "clonotype_id",
        method = list(abdiv::simpson, abdiv::mcintosh_d)
      ),
    "Names must be included"
  )

  res <- vdj_so |>
    calc_diversity(data_col = "cdr3", method = abdiv::simpson)

  nms <- paste0("cdr3_simpson_", c("diversity"))

  expect_true(all(nms %in% colnames(res@meta.data)))
})

