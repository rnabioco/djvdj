
# Check mutate_meta with function
test_that("mutate_meta function input", {
  res <- tiny_vdj %>%
    mutate_meta(select, .cell_id, orig.ident)

  expect_true(colnames(res@meta.data) == "orig.ident")
  expect_s4_class(res, "Seurat")
})

# Check mutate_meta with formula
test_that("mutate_meta formula input", {
  res <- tiny_vdj %>%
    mutate_meta(~ select(.x, .cell_id, orig.ident))

  expect_true(colnames(res@meta.data) == "orig.ident")
  expect_s4_class(res, "Seurat")
})

# Check mutate_meta overwrite column
test_that("mutate_meta column overwrite", {
  res <- tiny_vdj %>%
    mutate_meta(mutate, orig.ident = "HELLO")

  expect_true(all(res$orig.ident == "HELLO"))
  expect_s4_class(res, "Seurat")
})

# Check mutate_meta bad sobj_in
test_that("mutate_meta bad sobj_in", {
  expect_error(
    tiny_dat %>%
      mutate_meta(select, -orig.ident)
  )
})

# Check mutate_meta bad function
test_that("mutate_meta bad .fun", {
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
  res <- tiny_vdj %>%
    mutate_vdj(NEW = paste0(unique(chains), collapse = "_"))

  expect_s4_class(res, "Seurat")                       # class
  expect_identical(colnames(res), colnames(tiny_vdj))  # cells in object
})

# test_that("Default separator", {
#   res <- tiny_vdj %>%
#     mutate_vdj(
#       NEW = str_c(chains, collapse = "_"),
#       sep = ";"
#     )
#
#   old_nms <- str_replace_all(tiny_vdj$chains, ";", "_")
#
#   expect_s4_class(res, "Seurat")                       # class
#   expect_identical(colnames(res), colnames(tiny_vdj))  # cells in object
#   expect_identical(unname(res$NEW), old_nms)           # new column
# })

# Check NULL separator
test_that("mutate_vdj NULL sep", {
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

# Check bad command
test_that("mutate_vdj bad command", {
  expect_error(
    tiny_vdj %>%
      mutate_vdj(
        NEW = paste0(chains, nCount_RNA, sep = "_"),
        sep = ";"
      )
  )
})
