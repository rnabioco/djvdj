
# Check all cluster_seqs args
arg_lst <- list(
  input       = list(vdj_so, vdj_sce),
  data_col    = "cdr3",
  chain       = "IGH",
  method      = c("louvain", "leiden"),
  k           = c(5, 10),
  resolution  = c(0.1, 5),
  return_df   = c(TRUE, FALSE)
)

test_all_args(
  arg_lst = arg_lst,
  .fn     = cluster_seqs,
  desc    = "cluster_seqs args",
  chk     = expect_silent
)

# Check all plot_seq_motifs args
arg_lst <- list(
  input       = list(vdj_so, vdj_sce),
  data_col    = "cdr3",
  chain       = "IGH",
  cluster_col = list(NULL, "seurat_clusters"),
  plot_colors = list(NULL, test_cols),
  plot_lvls   = list(NULL, test_lvls),
  width       = c(0.5, 10),
  align_end   = c("3", "5"),
  facet_rows  = list(NULL, 2)
)

test_all_args(
  arg_lst = arg_lst,
  .fn     = plot_seq_motifs,
  desc    = "plot_seq_motifs args",
  chk     = expect_silent
)

