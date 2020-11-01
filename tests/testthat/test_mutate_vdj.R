
test_that("Check Seurat output", {
  res <- tiny_vdj %>%
    mutate_vdj(NEW = paste0(unique(chains), collapse = "_"))

  expect_s4_class(res, "Seurat")                       # class
  expect_identical(colnames(res), colnames(tiny_vdj))  # cells in object
})

test_that("Default separator", {
  res <- tiny_vdj %>%
    mutate_vdj(
      NEW = stringr::str_c(chains, collapse = "_"),
      sep = ";"
    )

  old_nms <- stringr::str_replace_all(tiny_vdj$chains, ";", "_")

  expect_s4_class(res, "Seurat")                       # class
  expect_identical(colnames(res), colnames(tiny_vdj))  # cells in object
  expect_identical(unname(res$NEW), old_nms)           # new column
})

test_that("NULL separator", {
  res <- tiny_vdj %>%
    mutate_vdj(
      NEW = paste0(chains, nCount_RNA, sep = "_"),
      sep = NULL
    )

  old_nms <- paste0(tiny_vdj$chains, tiny_vdj$nCount_RNA, sep = "_")

  expect_s4_class(res, "Seurat")                       # class
  expect_identical(colnames(res), colnames(tiny_vdj))  # cells in object
  expect_identical(unname(res$NEW), old_nms)           # new column
})

test_that("Default separator with bad command", {
  fn <- function() {
    tiny_vdj %>%
      mutate_vdj(
        NEW = paste0(chains, nCount_RNA, sep = "_"),
        sep = ";"
      )
  }

  expect_error(fn())
})
