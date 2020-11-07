
# Check all calc_usage arguments
arg_lst <- list(
  sobj_in     = list(tiny_vdj),
  gene_cols   = list("v_gene", "d_gene", "j_gene", "c_gene", c("v_gene", "j_gene")),
  cluster_col = list(NULL, "seurat_clusters"),
  chain       = list(NULL, "IGH", "IGL", "IGK"),
  chain_col   = "chains"
)

test_all_args(
  arg_lst = arg_lst,
  .fn     = calc_usage,
  ttl     = "calc_usage args",
  chk_fn  = expect_s3_class,
  chk_arg = "tbl"
)
