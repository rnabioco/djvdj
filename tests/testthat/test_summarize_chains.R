
# Check all summarize_chains arguments
args_lst <- list(
  sobj_in      = list(tiny_vdj),
  data_cols    = list("reads", "umis", c("reads", "umis")),
  fn           = list(mean, median),
  chain_col    = list(NULL, "chains"),
  include_cols = list(NULL, "seurat_clusters")
)

test_all_args(
  lst     = args_lst,
  .fn     = summarize_chains,
  ttl     = "summarize_chains args",
  chk_fn  = expect_s3_class,
  chk_arg = "tbl"
)
