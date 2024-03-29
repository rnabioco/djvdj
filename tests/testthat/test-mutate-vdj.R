# Check mutate_meta with function
test_that("mutate_meta function input", {
  res <- vdj_sce |>
    mutate_meta(select, .cell_id, orig.ident)

  expect_true(colnames(res@colData) == "orig.ident")
  expect_s4_class(res, "SingleCellExperiment")
})

# Check mutate_meta with formula
test_that("mutate_meta formula input", {
  res <- vdj_sce |>
    mutate_meta(~ select(.x, .cell_id, orig.ident))

  expect_true(colnames(res@colData) == "orig.ident")
  expect_s4_class(res, "SingleCellExperiment")
})

# Check mutate_meta overwrite column
test_that("mutate_meta column overwrite", {
  res <- vdj_sce |>
    mutate_meta(mutate, orig.ident = "HELLO")

  expect_true(all(res$orig.ident == "HELLO"))
  expect_s4_class(res, "SingleCellExperiment")
})

# Check mutate_meta bad function
test_that("mutate_meta bad fn", {
  err <- "`fn` must be a function or formula"

  expect_error(
    mutate_meta(tiny_so, "BAD_FUN"),
    err
  )

  expect_error(
    mutate_meta(tiny_so, 0),
    err
  )

  expect_error(
    mutate_meta(tiny_sce, as.factor("BAD_FUN")),
    err
  )

  # Do not allow user to remove cells from meta.data
  expect_error(
    vdj_sce |>
      mutate_meta(dplyr::filter, nCount_RNA == 1),
    "does not contain the same cells"
  )
})
# Check SingleCellExperiment output
test_that("mutate_vdj SingleCellExperiment out", {
  res <- vdj_sce |>
    mutate_vdj(NEW = paste0(unique(chains), collapse = "_"))

  expect_s4_class(res, "SingleCellExperiment")                     # class
  expect_identical(colnames(res), colnames(vdj_sce))  # cells in object

  expect_identical(
    select(as.data.frame(res@colData), -NEW),
    as.data.frame(vdj_sce@colData)
  )

  res <- vdj_sce |>
    mutate_vdj(NEW = stringr::str_c(chains, collapse = ";"))

  expect_identical(vdj_sce$chains, res$NEW)
})

# Check data.frame output
test_that("mutate_vdj df out", {
  res <- vdj_sce |>
    mutate_vdj(
      NEW       = paste0(unique(chains), collapse = "_"),
      return_df = TRUE
    )

  expect_s3_class(res, "data.frame")                 # class
  expect_identical(rownames(res), colnames(vdj_sce))  # cells in object

  expect_identical(
    select(res, -NEW),
    as.data.frame(vdj_sce@colData)
  )
})

# Check data.frame input
test_that("mutate_vdj df in", {
  res <- vdj_sce@colData |>
    mutate_vdj(NEW = paste0(unique(chains), collapse = "_"))

  expect_s3_class(res, "data.frame")                 # class
  expect_identical(rownames(res), colnames(vdj_sce))  # cells in object

  expect_identical(
    select(res, -NEW), as.data.frame(vdj_sce@colData)
  )
})

test_that("Default separator", {
  res <- vdj_sce |>
    mutate_vdj(
      NEW = ifelse(all(is.na(chains)), NA, paste0(chains, collapse = "_")),
      sep = ";"
    )

  old_nms <- unname(gsub(";", "_", vdj_sce$chains))

  expect_s4_class(res, "SingleCellExperiment")                     # class
  expect_identical(colnames(res), colnames(vdj_sce))  # cells in object
  expect_identical(unname(res$NEW), old_nms)         # new column
})

# Check NULL separator
test_that("mutate_vdj NULL sep", {
  res <- vdj_sce |>
    mutate_vdj(
      NEW = paste0(chains, nCount_RNA, sep = "_"),
      sep = NULL
    )

  old_nms <- paste0(vdj_sce$chains, vdj_sce$nCount_RNA, sep = "_")

  expect_s4_class(res, "SingleCellExperiment")                       # class
  expect_identical(colnames(res), colnames(vdj_sce))    # cells in object
  expect_identical(unname(res$NEW), old_nms)           # new column
})

# Check bad command
test_that("mutate_vdj bad command", {
  expect_error(
    vdj_sce |>
      mutate_vdj(
        NEW = paste0(chains, nCount_RNA, sep = "_"),
        sep = ";"
      )
  )
})

# Check mutate_vdj calculations
test_that("mutate_vdj calcs", {
  res <- vdj_sce@colData |>
    mutate_vdj(NEW_umis = mean(umis))

  x <- vdj_sce$umis |>
    strsplit(";") |>
    purrr::map_dbl(~ mean(as.numeric(.x))) |>
    unname()

  expect_identical(x, res$NEW_umis)

  res <- vdj_sce@colData |>
    mutate_vdj(
      NEW1 = ifelse(all(is.na(v_gene)), NA, paste0(v_gene, collapse = "_")),
      NEW2 = stringr::str_c(v_gene, collapse = "_")
    )

  x <- vdj_sce$v_gene |>
    strsplit(";") |>
    purrr::map_chr(~ ifelse(all(is.na(.x)), NA, paste0(.x, collapse = "_"))) |>
    unname()

  expect_identical(x, res$NEW1)
  expect_identical(x, res$NEW2)

  # Chcek that function and formula inputs return same result
  x <- vdj_sce |>
    mutate_meta(
      mutate,
      NEW = nCount_RNA + nFeature_RNA
    )

  y <- vdj_sce |>
    mutate_meta(~ mutate(.x, NEW = nCount_RNA + nFeature_RNA))

  expect_identical(x, y)
  expect_identical(x$NEW, x$nCount_RNA + x$nFeature_RNA)
})

# Check SCE input
test_that("summarize_vdj SCE out", {
  res <- vdj_sce |>
    summarize_vdj(
      data_cols = "umis",
      col_names = "NEW_{.col}"
    )

  expect_s4_class(res, "SingleCellExperiment")        # class
  expect_identical(colnames(res), colnames(vdj_sce))  # cells in object
  expect_identical(res@colData[, -ncol(res@colData)], vdj_sce@colData)
})

# Check data.frame output
test_that("summarize_vdj df out", {
  res <- vdj_sce |>
    summarize_vdj(
      data_cols = "umis",
      col_names = "NEW_{.col}",
      return_df = TRUE
    )

  expect_s3_class(res, "data.frame")                 # class
  expect_identical(rownames(res), colnames(vdj_sce))  # cells in object
  expect_identical(
    select(res, -starts_with("NEW_")), as.data.frame(vdj_sce@colData)
  )
})

# Check data.frame input
test_that("summarize_vdj df in", {
  res <- vdj_sce@colData |>
    summarize_vdj(
      data_cols = "umis",
      col_names = "NEW_{.col}"
    )

  expect_s3_class(res, "data.frame")                 # class
  expect_identical(rownames(res), colnames(vdj_sce))  # cells in object
  expect_identical(
    select(res, -starts_with("NEW_")), as.data.frame(vdj_sce@colData)
  )
})

# Check calculations
test_that("summarize_vdj calcs", {
  res <- vdj_sce@colData |>
    summarize_vdj(
      data_cols = "umis",
      col_names = "NEW_{.col}"
    )

  x <- vdj_sce$umis |>
    strsplit(";") |>
    purrr::map_dbl(~ mean(as.numeric(.x))) |>
    unname()

  expect_identical(x, res$NEW_umis)

  res <- vdj_sce@colData |>
    summarize_vdj(
      data_cols  = "v_gene",
      fn        = ~ ifelse(all(is.na(.x)), NA, paste0(.x, collapse = "_")),
      col_names = "NEW_{.col}"
    )

  x <- vdj_sce$v_gene |>
    strsplit(";") |>
    purrr::map_chr(~ ifelse(all(is.na(.x)), NA, paste0(.x, collapse = "_"))) |>
    unname()

  expect_identical(x, res$NEW_v_gene)

  # Check chain filtering
  res <- vdj_sce@colData |>
    summarize_vdj(
      data_cols  = "umis",
      col_names = "NEW_{.col}",
      chain     = "IGK"
    )

  x <- strsplit(vdj_sce$umis, ";") |>
    purrr::map(as.numeric)

  y <- strsplit(vdj_sce$chains, ";")

  z <- purrr::map2_dbl(
    x, y,
    ~ ifelse(all(is.na(.x)), NA, mean(.x[.y == "IGK"]))
  ) |>
    unname()

  expect_equal(z, res$NEW_umis)
  expect_identical(res$chains, unname(vdj_sce$chains))
})

# Check bad command
test_that("summarize_vdj return length > 1", {
  res_1 <- vdj_sce |>
    summarize_vdj(
      data_cols = "umis",
      fn = ~ .x[.x > 30]
    )

  res_2 <- vdj_sce |>
    summarize_vdj(
      data_cols = "umis",
      fn = ~ ifelse(length(.x[.x > 30]) > 0, paste0(.x[.x > 30], collapse = ";"), NA)
    )

  expect_identical(res_1, res_2)

  res_1 <- vdj_sce |>
    summarize_vdj(
      data_cols = "chains",
      fn = ~ .x[.x == "IGH"]
    )

  res_2 <- vdj_sce |>
    summarize_vdj(
      data_cols = "chains",
      fn = ~ ifelse(length(.x[.x == "IGH"]) > 0, paste0(.x[.x == "IGH"], collapse = ";"), NA)
    )

  expect_identical(res_1, res_2)
})

# Check bad column chain filtering
test_that("summarize_vdj bad column chain filtering", {
  expect_warning(
    res <- vdj_sce |>
      summarize_vdj(
        data_cols = "nCount_RNA",
        chain = "IGK"
      ),
    "columns do not contain per-chain"
  )

  expect_identical(res, vdj_sce)
})

# Check .prepare_meta
test_that(".prepare_meta", {
  res <- vdj_sce |>
    fetch_vdj(unnest = FALSE)

  expect_error(
    .prepare_meta(vdj_sce, res, ".cell_id"),
    "meta.data cannot include list-cols"
  )

  res <- vdj_sce@colData |>
    as.data.frame() |>
    filter(orig.ident == "avid_1")

  expect_error(
    .prepare_meta(vdj_sce, res, ".cell_id"),
    "meta.data does not contain the same cells"
  )
})

