
# Test inputs
test_cols <- c(
  "#E69F00", "#56B4E9", "#009E73",
  "#F0E442", "#d7301f", "#0072B2",
  "#D55E00", "#6A51A3", "#CC79A7",
  "#999999", "#875C04", "#000000"
)

test_lvls <- unique(tiny_vdj$seurat_clusters)

tiny_dat <- tiny_vdj@meta.data %>%
  as_tibble(rownames = ".cell_id")

# Check all plot_features arguments
lst <- list(
  sobj_in     = list(tiny_vdj, tiny_dat),
  feature     = c("seurat_clusters", "nCount_RNA"),
  data_slot   = c("data", "counts"),
  plot_colors = list(NULL, test_cols),
  feat_lvls   = list(NULL, test_lvls),
  facet_id    = list(NULL, "orig.ident"),
  facet_lvls  = list(NULL, "AVID-seq"),
  min_pct     = list(NULL, 0.05),
  max_pct     = list(NULL, 0.95),
  lm_line     = c(TRUE, FALSE)
)

test_all_args(
  lst     = lst,
  .fn     = plot_features,
  ttl     = "plot_features args",
  chk_fn  = expect_s3_class,
  chk_arg = "ggplot"
)

# Check all plot_cell_count arguments
lst <- list(
  sobj_in     = list(tiny_vdj, tiny_dat),
  x           = "orig.ident",
  fill_col    = list(NULL, "seurat_clusters"),
  facet_col   = list(NULL, "orig.ident"),
  yaxis       = c("fraction", "count"),
  plot_colors = list(NULL, test_cols),
  plot_lvls   = list(NULL, "AVID-seq"),
  order_count = c(TRUE, FALSE),
  n_label     = c(TRUE, FALSE),
  label_aes   = list(list(), list(size = 2))
)

test_all_args(
  lst     = lst,
  .fn     = plot_cell_count,
  ttl     = "plot_cell_count args",
  chk_fn  = expect_s3_class,
  chk_arg = "ggplot"
)

# Check all plot_reads arguments
lst <- list(
  sobj_in     = list(tiny_vdj),
  data_cols   = list("reads", "umis", c("reads", "umis")),
  chain_col   = list(NULL, "chains"),
  cluster_col = list(NULL, "seurat_clusters"),
  type        = c("violin", "histogram", "density"),
  plot_colors = list(NULL, test_cols),
  plot_lvls   = list(NULL, test_lvls)
)

test_all_args(
  lst     = lst,
  .fn     = plot_reads,
  ttl     = "plot_reads args",
  chk_fn  = expect_s3_class,
  chk_arg = "ggplot"
)

# Check all plot_abundance arguments for line plot
ab_lst <- list(
  sobj_in       = list(tiny_vdj),
  clonotype_col = "cdr3",
  cluster_col   = list(NULL, "seurat_clusters"),
  type          = "line",
  label_col     = list(NULL, "cdr3"),
  yaxis         = c("percent", "frequency"),
  plot_colors   = list(NULL, test_cols),
  plot_lvls     = list(NULL, test_lvls),
  label_aes     = list(list(), list(size = 2))
)

test_all_args(
  lst     = ab_lst,
  .fn     = plot_abundance,
  ttl     = "plot_abundance args",
  chk_fn  = expect_s3_class,
  chk_arg = "ggplot"
)

# Check all plot_abundance arguments for bar plot
ab_lst$type <- "bar"
ab_lst$label_col <- "cdr3"

test_all_args(
  lst     = ab_lst,
  .fn     = plot_abundance,
  ttl     = "plot_abundance args",
  chk_fn  = expect_s3_class,
  chk_arg = "ggplot"
)

# Check all plot_diversity arguments
mets <- abdiv::alpha_diversities %>%
  map(~ eval(parse(text = paste0("abdiv::", .x))))

lst <- list(
  sobj_in       = list(tiny_vdj),
  clonotype_col = "cdr3",
  cluster_col   = list(NULL, "seurat_clusters"),
  method        = mets,
  plot_colors   = list(NULL, test_cols),
  plot_lvls     = list(NULL, test_lvls)
)

test_all_args(
  lst     = lst,
  .fn     = plot_diversity,
  ttl     = "plot_diversity args",
  chk_fn  = expect_s3_class,
  chk_arg = "ggplot"
)

# Check all plot_similarity arguments
mets <- abdiv::beta_diversities %>%
  map(~ eval(parse(text = paste0("abdiv::", .x))))

lst <- list(
  sobj_in       = list(tiny_vdj),
  clonotype_col = "cdr3",
  cluster_col   = "seurat_clusters",
  method        = mets,
  plot_colors   = list(NULL, test_cols)
)

test_all_args(
  lst     = lst,
  .fn     = plot_similarity,
  ttl     = "plot_similarity args",
  chk_fn  = expect_s3_class,
  chk_arg = "ggplot"
)

# Check all plot_usage arguments
lst <- list(
  sobj_in     = list(tiny_vdj),
  gene_cols   = list("v_gene", "d_gene", "j_gene", "c_gene", c("v_gene", "j_gene")),
  cluster_col = list(NULL, "seurat_clusters"),
  chain       = list(NULL, "IGH", "IGL", "IGK"),
  chain_col   = "chains",
  type        = c("heatmap", "bar"),
  plot_colors = list(NULL, test_cols),
  plot_lvls   = list(NULL, test_lvls),
  yaxis       = c("percent", "frequency")
)

test_all_args(
  lst = lst,
  .fn = plot_usage,
  ttl = "plot_usage args"
)
