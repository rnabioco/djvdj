# Test data
df_1 <- vdj_so@meta.data

df_2 <- vdj_so@meta.data |>
  as_tibble(rownames = ".cell_id")

# Check fetch_vdj arguments
arg_lst <- list(
  input         = list(vdj_so, vdj_sce, df_1, df_2),
  data_cols      = list(NULL, c("umis", "reads")),
  clonotype_col = list(NULL, "clonotype_id"),
  unnest        = c(TRUE, FALSE),
  sep           = c(NULL, ";")
)

test_all_args(
  arg_lst = arg_lst,
  .fn     = fetch_vdj,
  desc    = "fetch_vdj all args",
  chk     = expr(expect_true(is.data.frame(.res)))
)

# Check sep
test_that("fetch_vdj sep", {
  res <- vdj_so |>
    fetch_vdj(
      data_cols      = c("umis", "reads"),
      clonotype_col = "clonotype_id",
      filter_cells  = FALSE,
      sep           = NULL
    )

  expect_identical(res, as_tibble(vdj_so@meta.data, rownames = ".cell_id"))

  expect_warning(
    res <- vdj_so |>
      fetch_vdj(
        data_cols      = c("umis", "reads"),
        clonotype_col = "clonotype_id",
        filter_cells  = FALSE,
        sep           = "BAD"
      ),
    "was not found in umis or reads"
  )
})

# Check filter_cells
test_that("fetch_vdj filter_cells", {
  expect_error(
    vdj_so |>
      fetch_vdj(
        # data_cols     = c("umis", "reads"),
        filter_cells = TRUE
      ),
    "`clonotype_col` must be provided"
  )

  res <- vdj_so |>
    fetch_vdj(
      data_cols      = c("umis", "reads"),
      clonotype_col = "clonotype_id",
      filter_cells  = FALSE
    )

  expect_identical(unique(res$.cell_id), colnames(vdj_so))

  res <- vdj_so |>
    fetch_vdj(
      data_cols      = c("umis", "reads"),
      clonotype_col = "clonotype_id",
      filter_cells  = TRUE
    )

  no_na <- vdj_so@meta.data |>
    filter(!is.na(clonotype_id)) |>
    rownames()

  expect_identical(no_na, unique(res$.cell_id))
})

# Check .filter_chains
test_that(".filter_chains", {
  dat <- vdj_so |>
    fetch_vdj(unnest = FALSE)

  res <- dat |>
    .filter_chains(
      data_cols = c("umis", "reads"),
      chain = NULL
    )

  expect_identical(res, dat)

  expect_error(
    res <- dat |>
      .filter_chains(
        data_cols = c("umis", "nCount_RNA"),
        chain = "IGK"
      ),
    "Some columns do not contain per-chain"
  )

  dat <- vdj_so |>
    fetch_vdj(unnest = TRUE)

  res <- dat |>
    .filter_chains(
      data_cols = c("umis", "reads"),
      chain = "IGK"
    )

  expect_identical(res, filter(dat, chains %in% "IGK"))

  res <- dat |>
    .filter_chains(
      data_cols = "chains",
      chain    = "IGK"
    )

  new_chains <- res$chains |>
    reduce(c) |>
    unique() |>
    na.omit() |>
    as.character()

  expect_identical(new_chains, "IGK")

  res <- dat |>
    .filter_chains(
      data_cols = "chains",
      chain    = c("IGH", "IGL")
    )

  new_chains <- res$chains |>
    reduce(c) |>
    unique() |>
    na.omit() |>
    as.character() |>
    sort()

  expect_identical(new_chains, c("IGH", "IGL"))
})

# Check .filter_chains
test_that(".filter_chains", {
  dat <- vdj_so |>
    fetch_vdj(unnest = FALSE)

  res <- dat |>
    .filter_chains(
      data_cols = c("umis", "reads"),
      chain = NULL
    )

  expect_identical(res, dat)

  expect_error(
    res <- dat |>
      .filter_chains(
        data_cols = c("umis", "nCount_RNA"),
        chain = "IGK"
      ),
    "Some columns do not contain per-chain"
  )

  dat <- vdj_so |>
    fetch_vdj(unnest = TRUE)

  res <- dat |>
    .filter_chains(
      data_cols = c("umis", "reads"),
      chain = "IGK"
    )

  expect_identical(res, filter(dat, chains %in% "IGK"))

  res <- dat |>
    .filter_chains(
      data_cols = "chains",
      chain    = "IGK"
    )

  new_chains <- res$chains |>
    reduce(c) |>
    unique() |>
    na.omit() |>
    as.character()

  expect_identical(new_chains, "IGK")

  res <- dat |>
    .filter_chains(
      data_cols = "chains",
      chain    = c("IGH", "IGL")
    )

  new_chains <- res$chains |>
    reduce(c) |>
    unique() |>
    na.omit() |>
    as.character() |>
    sort()

  expect_identical(new_chains, c("IGH", "IGL"))
})

# Check .get_vdj_cols
test_that(".get_vdj_cols", {
  clmns <- c("reads", "nCount_RNA")

  res <- vdj_so@meta.data |>
    .get_vdj_cols(
      clone_col = NULL,
      cols_in = clmns,
      sep = ";"
    )

  expect_identical(res$vdj, clmns)
  expect_identical(res$sep, "reads")

  res <- vdj_so@meta.data |>
    .get_vdj_cols(
      clone_col = NULL,
      cols_in = clmns,
      sep = "ZZZ"
    )

  expect_identical(res$vdj, clmns)
  expect_null(res$sep)
})

