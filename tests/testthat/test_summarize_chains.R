
test_that("Check all summarize_chains arguments", {
  args_df <- list(
    data_cols    = list("reads", "umis", c("reads", "umis")),
    fn           = list(mean, median),
    chain_col    = list(NULL, "chains"),
    include_cols = list(NULL, "seurat_clusters")
  ) %>%
    expand.grid()

  res <- pmap(
    args_df,
    check_args,
    .fn     = summarize_chains,
    sobj_in = tiny_vdj
  )

  walk(res, expect_s3_class, "tbl")
})
