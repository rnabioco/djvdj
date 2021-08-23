# Test data
df_1 <- vdj_so@meta.data

df_2 <- vdj_so@meta.data %>%
  as_tibble(rownames = ".cell_id")

# Check all summarize_chains arguments
arg_lst <- list(
  input        = list(vdj_so, vdj_sce, df_1, df_2),
  data_cols    = list("reads", "umis", c("reads", "umis")),
  fn           = list(mean, median),
  chain_col    = list(NULL, "chains"),
  include_cols = list(NULL, "seurat_clusters")
)

test_all_args(
  arg_lst = arg_lst,
  .fn     = summarize_chains,
  desc    = "summarize_chains args",
  chk     = expr(expect_s3_class(.res, "tbl"))
)
