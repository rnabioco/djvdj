# Check mutate_meta with function
test_that("mutate_meta function input", {
  res <- vdj_so %>%
    mutate_meta(select, .cell_id, orig.ident)

  expect_true(colnames(res@meta.data) == "orig.ident")
  expect_s4_class(res, "Seurat")
})

# Check mutate_meta with formula
test_that("mutate_meta formula input", {
  res <- vdj_so %>%
    mutate_meta(~ select(.x, .cell_id, orig.ident))

  expect_true(colnames(res@meta.data) == "orig.ident")
  expect_s4_class(res, "Seurat")
})

# Check mutate_meta overwrite column
test_that("mutate_meta column overwrite", {
  res <- vdj_so %>%
    mutate_meta(mutate, orig.ident = "HELLO")

  expect_true(all(res$orig.ident == "HELLO"))
  expect_s4_class(res, "Seurat")
})

# Check mutate_meta bad function
test_that("mutate_meta bad fn", {
  expect_error(
    tiny_so %>%
      mutate_meta("BAD_FUN")
  )

  expect_error(
    tiny_so %>%
      mutate_meta(0)
  )

  expect_error(
    tiny_so %>%
      mutate_meta(as.factor("BAD_FUN"))
  )
})

# Check Seurat output
test_that("mutate_vdj Seurat out", {
  res <- vdj_so %>%
    mutate_vdj(NEW = paste0(unique(chains), collapse = "_"))

  expect_s4_class(res, "Seurat")                     # class
  expect_identical(colnames(res), colnames(vdj_so))  # cells in object

  expect_identical(select(res@meta.data, -NEW), vdj_so@meta.data)
})

# Check data.frame output
test_that("mutate_vdj df out", {
  res <- vdj_so %>%
    mutate_vdj(
      NEW       = paste0(unique(chains), collapse = "_"),
      return_df = TRUE
    )

  expect_s3_class(res, "data.frame")                 # class
  expect_identical(rownames(res), colnames(vdj_so))  # cells in object

  expect_identical(select(res, -NEW), vdj_so@meta.data)
})

# Check data.frame input
test_that("mutate_vdj df in", {
  res <- vdj_so@meta.data %>%
    mutate_vdj(NEW = paste0(unique(chains), collapse = "_"))

  expect_s3_class(res, "data.frame")                 # class
  expect_identical(rownames(res), colnames(vdj_so))  # cells in object

  expect_identical(select(res, -NEW), vdj_so@meta.data)
})

test_that("Default separator", {
  res <- vdj_so %>%
    mutate_vdj(
      NEW = ifelse(all(is.na(chains)), NA, paste0(chains, collapse = "_")),
      sep = ";"
    )

  old_nms <- unname(gsub(";", "_", vdj_so$chains))

  expect_s4_class(res, "Seurat")                     # class
  expect_identical(colnames(res), colnames(vdj_so))  # cells in object
  expect_identical(unname(res$NEW), old_nms)         # new column
})

# Check NULL separator
test_that("mutate_vdj NULL sep", {
  res <- vdj_so %>%
    mutate_vdj(
      NEW = paste0(chains, nCount_RNA, sep = "_"),
      sep = NULL
    )

  old_nms <- paste0(vdj_so$chains, vdj_so$nCount_RNA, sep = "_")

  expect_s4_class(res, "Seurat")                       # class
  expect_identical(colnames(res), colnames(vdj_so))    # cells in object
  expect_identical(unname(res$NEW), old_nms)           # new column
})

# Check bad command
test_that("mutate_vdj bad command", {
  expect_error(
    vdj_so %>%
      mutate_vdj(
        NEW = paste0(chains, nCount_RNA, sep = "_"),
        sep = ";"
      )
  )
})





# Check Seurat output
test_that("summarize_vdj Seurat out", {
  res <- vdj_so %>%
    summarize_vdj(
      vdj_cols = "umis",
      col_names = "NEW_{.col}"
    )

  expect_s4_class(res, "Seurat")                     # class
  expect_identical(colnames(res), colnames(vdj_so))  # cells in object
  expect_identical(select(res@meta.data, -starts_with("NEW_")), vdj_so@meta.data)
})

# Check SCE input
test_that("summarize_vdj SCE out", {
  res <- vdj_sce %>%
    summarize_vdj(
      vdj_cols = "umis",
      col_names = "NEW_{.col}"
    )

  expect_s4_class(res, "SingleCellExperiment")        # class
  expect_identical(colnames(res), colnames(vdj_sce))  # cells in object
  expect_identical(res@colData[, -ncol(res@colData)], vdj_sce@colData)
})

# Check data.frame output
test_that("summarize_vdj df out", {
  res <- vdj_so %>%
    summarize_vdj(
      vdj_cols = "umis",
      col_names = "NEW_{.col}",
      return_df = TRUE
    )

  expect_s3_class(res, "data.frame")                 # class
  expect_identical(rownames(res), colnames(vdj_so))  # cells in object
  expect_identical(select(res, -starts_with("NEW_")), vdj_so@meta.data)
})

# Check data.frame input
test_that("summarize_vdj df in", {
  res <- vdj_so@meta.data %>%
    summarize_vdj(
      vdj_cols = "umis",
      col_names = "NEW_{.col}"
    )

  expect_s3_class(res, "data.frame")                 # class
  expect_identical(rownames(res), colnames(vdj_so))  # cells in object
  expect_identical(select(res, -starts_with("NEW_")), vdj_so@meta.data)
})

# Check calculations
test_that("summarize_vdj calcs", {
  res <- vdj_so@meta.data %>%
    summarize_vdj(
      vdj_cols = "umis",
      col_names = "NEW_{.col}"
    )

  x <- vdj_so$umis %>%
    strsplit(";") %>%
    purrr::map_dbl(~ mean(as.numeric(.x))) %>%
    unname()

  expect_identical(x, res$NEW_umis)

  res <- vdj_so@meta.data %>%
    summarize_vdj(
      vdj_cols  = "v_gene",
      fn        = ~ ifelse(all(is.na(.x)), NA, paste0(.x, collapse = "_")),
      col_names = "NEW_{.col}"
    )

  x <- vdj_so$v_gene %>%
    strsplit(";") %>%
    purrr::map_chr(~ ifelse(all(is.na(.x)), NA, paste0(.x, collapse = "_"))) %>%
    unname()

  expect_identical(x, res$NEW_v_gene)

  # Check chain filtering
  res <- vdj_so@meta.data %>%
    summarize_vdj(
      vdj_cols  = "umis",
      col_names = "NEW_{.col}",
      chain     = "IGK"
    )

  x <- strsplit(vdj_so$umis, ";") %>%
    purrr::map(as.numeric)

  y <- strsplit(vdj_so$chains, ";")

  z <- purrr::map2_dbl(x, y, ~ ifelse(all(is.na(.x)), NA, mean(.x[.y == "IGK"]))) %>%
    unname()

  expect_equal(z, res$NEW_umis)
  expect_identical(res$chains, unname(vdj_so$chains))
})

# Check bad command
test_that("summarize_vdj bad command", {
  expect_error(
    res <- vdj_so %>%
      summarize_vdj(
        vdj_cols = "umis",
        fn = ~ .x[chains == "IGK"]
      )
  )
})
