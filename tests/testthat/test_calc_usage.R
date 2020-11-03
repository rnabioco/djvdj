
test_that("Check all calc_usage arguments", {
  args_df <- list(
    gene_cols   = list("v_gene", "d_gene", "j_gene", "c_gene", c("v_gene", "j_gene")),
    cluster_col = list(NULL, "seurat_clusters"),
    chain       = list(NULL, "IGH", "IGL", "IGK")
  ) %>%
    expand.grid()

  res <- pmap(
    args_df,
    check_args,
    .fn       = calc_usage,
    sobj_in   = tiny_vdj,
    chain_col = "chains"
  )

  walk(res, expect_s3_class, "tbl")
})
