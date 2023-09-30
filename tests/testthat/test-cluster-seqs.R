
test_cols <- c(
  "#E69F00", "#56B4E9", "#009E73",
  "#F0E442", "#d7301f", "#0072B2",
  "#D55E00", "#6A51A3", "#CC79A7",
  "#999999", "#875C04", "#000000"
)

test_lvls <- unique(vdj_sce$seurat_clusters) |>
  as.character() |>
  rev()

# Check all cluster_seqs args
arg_lst <- list(
  input       = list(vdj_sce),
  data_col    = "cdr3",
  chain       = list(NULL, "IGK"),
  method      = c("louvain", "leiden"),
  k           = 10,
  resolution  = list(5, c(0.1, 5)),
  return_df   = c(TRUE, FALSE)
)

test_all_args(
  arg_lst = arg_lst,
  .fn     = cluster_sequences,
  desc    = "cluster_seqs args",
  chk     = expect_silent
)

# Check all plot_seq_motifs args
arg_lst <- list(
  input       = list(vdj_sce),
  data_col    = "cdr3",
  chain       = "IGH",
  cluster_col = list(NULL, "seurat_clusters"),
  plot_colors = list(NULL, test_cols),
  plot_lvls   = list(NULL),
  width       = c(2),
  align_end   = c("3", "5"),
  panel_nrow  = list(NULL, 2)
)

test_all_args(
  arg_lst = arg_lst,
  .fn     = plot_motifs,
  desc    = "plot_seq_motifs args",
  chk     = expect_silent
)

