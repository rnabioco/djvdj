
# Check all calc_usage arguments
arg_lst <- list(
  input       = list(tiny_vdj),
  gene_cols   = list("v_gene", "d_gene", "j_gene", "c_gene", c("v_gene", "j_gene")),
  cluster_col = list(NULL, "seurat_clusters"),
  chain       = list(NULL, "IGH", "IGL", "IGK"),
  chain_col   = "chains"
)

test_all_args(
  arg_lst = arg_lst,
  .fn     = calc_usage,
  desc    = "calc_usage args",
  chk     = expr(expect_s3_class(.res, "tbl"))
)
