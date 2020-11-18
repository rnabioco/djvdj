
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

# Check all plot_features arguments for character feature
arg_lst <- list(
  x           = list("UMAP_1", c(x = "UMAP_1")),
  y           = list("UMAP_2", c(y = "UMAP_2")),
  sobj_in     = list(tiny_vdj, tiny_dat),
  feature     = list("seurat_clusters", c(clust = "seurat_clusters")),
  data_slot   = c("data", "counts"),
  plot_colors = list(NULL, test_cols),
  feat_lvls   = list(NULL, test_lvls),
  facet_col   = list(NULL, "orig.ident", c("orig.ident", "seurat_clusters")),
  facet_lvls  = list(NULL, "AVID-seq"),
  min_pct     = list(NULL, 0.05),
  max_pct     = list(NULL, 0.95),
  lm_line     = c(TRUE, FALSE)
)

test_all_args(
  arg_lst = arg_lst,
  .fn     = plot_features,
  ttl     = "plot_features args chr feat",
  chk     = expr(expect_s3_class(.res, "ggplot"))
)

# Check all plot_features arguments for numeric feature
arg_lst$feature   <- "nCount_RNA"
arg_lst$feat_lvls <- NULL

test_all_args(
  arg_lst = arg_lst,
  .fn     = plot_features,
  ttl     = "plot_features args num feat",
  chk     = expr(expect_s3_class(.res, "ggplot"))
)

# Check plot_features warning for numeric feature
test_that("plot_features warning num feat", {
  expect_warning(
    plot_features(
      sobj_in   = tiny_vdj,
      feature   = "nCount_RNA",
      feat_lvls = test_lvls
    )
  )
})

# Check plot_features error for same x and y
test_that("plot_features error same x y", {
  expect_error(
    plot_features(
      sobj_in   = tiny_vdj,
      x         = "UMAP_1",
      y         = "UMAP_1",
      feature   = "nCount_RNA",
      feat_lvls = test_lvls
    )
  )
})

# Check plot_features feature not found
arg_lst <- list(
  sobj_in = list(tiny_vdj, tiny_dat),
  feature = "BAD_FEATURE"
)

test_all_args(
  arg_lst = arg_lst,
  .fn     = plot_features,
  ttl     = "plot_features no feat found",
  chk     = expect_error
)

# Check all plot_cell_count arguments
arg_lst <- list(
  sobj_in     = list(tiny_vdj, tiny_dat),
  x           = "orig.ident",
  fill_col    = list(NULL, "seurat_clusters"),
  facet_col   = list(NULL, "orig.ident"),
  yaxis       = c("fraction", "counts"),
  plot_colors = list(NULL, test_cols),
  plot_lvls   = list(NULL, "AVID-seq"),
  n_label     = c(TRUE, FALSE),
  label_aes   = list(list(), list(size = 2))
)

test_all_args(
  arg_lst = arg_lst,
  .fn     = plot_cell_count,
  ttl     = "plot_cell_count args",
  chk     = expr(expect_s3_class(.res, "ggplot"))
)

# Check plot_cell_count bad yaxis
test_that("plot_cell_count bad yaxis", {
  expect_error(
    plot_cell_count(
      sobj_in = tiny_vdj,
      x       = "orig.ident",
      yaxis   = "BAD"
    )
  )
})

# Check all plot_reads arguments
arg_lst <- list(
  sobj_in     = list(tiny_vdj),
  data_cols   = list("reads", "umis", c("reads", "umis")),
  chain_col   = list(NULL, "chains"),
  cluster_col = list(NULL, "seurat_clusters"),
  type        = c("violin", "histogram", "density"),
  plot_colors = list(NULL, test_cols),
  plot_lvls   = list(NULL, test_lvls)
)

test_all_args(
  arg_lst = arg_lst,
  .fn     = plot_reads,
  ttl     = "plot_reads args",
  chk     = expr(expect_s3_class(.res, "ggplot"))
)

# Check plot_reads bad type
test_that("plot_reads bad type", {
  expect_error(
    plot_reads(
      sobj_in = tiny_vdj,
      type    = "BAD"
    )
  )
})

# Check all plot_abundance arguments for line plot
arg_lst <- list(
  sobj_in       = list(tiny_vdj),
  clonotype_col = "cdr3_nt",
  cluster_col   = list(NULL, "seurat_clusters"),
  type          = "line",
  label_col     = list(NULL, "cdr3"),
  color_col     = list(NULL, "seurat_clusters"),
  yaxis         = c("percent", "frequency"),
  plot_colors   = list(NULL, test_cols),
  plot_lvls     = list(NULL, test_lvls),
  label_aes     = list(list(), list(size = 2))
)

test_all_args(
  arg_lst = arg_lst,
  .fn     = plot_abundance,
  ttl     = "plot_abundance bar args",
  chk     = expr(expect_s3_class(.res, "ggplot"))
)

# Check all plot_abundance arguments for bar plot
arg_lst$type <- "bar"
arg_lst$label_col <- "cdr3"

test_all_args(
  arg_lst = arg_lst,
  .fn     = plot_abundance,
  ttl     = "plot_abundance line args",
  chk     = expr(expect_s3_class(.res, "ggplot"))
)

# Check plot_abundance axis labels
arg_lst <- list(
  sobj_in       = list(tiny_vdj),
  clonotype_col = "cdr3_nt",
  label_col     = "cdr3",
  yaxis         = "percent",
  type          = c("bar", "line")
)

test_all_args(
  arg_lst = arg_lst,
  .fn     = plot_abundance,
  ttl     = "plot_abundance axis labels",
  chk     = expr(expect_true(.res$label$y == "percent"))
)

arg_lst$yaxis <- "frequency"

test_all_args(
  arg_lst = arg_lst,
  .fn     = plot_abundance,
  ttl     = "plot_abundance axis labels",
  chk     = expr(expect_true(.res$label$y == "frequency"))
)

# Check plot_abundance bad yaxis
test_that("plot_abundance bad yaxis", {
  expect_error(
    plot_abundance(
      sobj_in       = tiny_vdj,
      type          = "line",
      clonotype_col = "cdr3_nt",
      yaxis         = "BAD"
    )
  )
})

# Check plot_abundance bad type
test_that("plot_abundance bad type", {
  expect_error(
    plot_abundance(
      sobj_in       = tiny_vdj,
      type          = "BAD",
      clonotype_col = "cdr3_nt"
    )
  )
})

# Check plot_abundance bad label_col
test_that("plot_abundance bad label_col", {
  expect_error(
    plot_abundance(
      sobj_in       = tiny_vdj,
      type          = "bar",
      clonotype_col = "cdr3_nt",
      label_col     = NULL
    )
  )
})

# Check all plot_diversity arguments
mets <- abdiv::alpha_diversities %>%
  map(~ eval(parse(text = paste0("abdiv::", .x))))

names(mets) <- abdiv::alpha_diversities

arg_lst <- list(
  sobj_in       = list(tiny_vdj),
  clonotype_col = "cdr3_nt",
  cluster_col   = list(NULL, "seurat_clusters"),
  method        = append(mets, list(mets)),
  plot_colors   = list(NULL, test_cols),
  plot_lvls     = list(NULL, test_lvls)
)

test_all_args(
  arg_lst = arg_lst,
  .fn     = plot_diversity,
  ttl     = "plot_diversity args",
  chk     = expr(expect_s3_class(.res, "ggplot"))
)

# Check plot_diversity bad method names
mets <- abdiv::alpha_diversities %>%
  map(~ eval(parse(text = paste0("abdiv::", .x))))

test_that("plot_diversity bad names", {
  expect_error(
    plot_diversity(
      sobj_in       = tiny_vdj,
      clonotype_col = "cdr3_nt",
      method        = mets
    )
  )
})

# Check all plot_similarity arguments
mets <- abdiv::beta_diversities %>%
  map(~ eval(parse(text = paste0("abdiv::", .x))))

arg_lst <- list(
  sobj_in       = list(tiny_vdj),
  clonotype_col = "cdr3_nt",
  cluster_col   = "seurat_clusters",
  method        = mets,
  plot_colors   = list(NULL, test_cols)
)

test_all_args(
  arg_lst = arg_lst,
  .fn     = plot_similarity,
  ttl     = "plot_similarity args",
  chk     = expr(expect_s3_class(.res, "ggplot"))
)

# Check plot_similarity bad clonotype_col
test_that("plot_similarity bad clonotype col", {
  expect_error(
    plot_similarity(
      sobj_in       = tiny_vdj,
      clonotype_col = NULL,
      cluster_col   = "orig.ident"
    )
  )
})

# Check plot_similarity bad cluster_col
test_that("plot_similarity bad cluster col", {
  expect_error(
    plot_similarity(
      sobj_in       = tiny_vdj,
      cluster_col   = NULL,
      clonotype_col = "cdr3_nt"
    )
  )
})

# Check all plot_usage arguments return ggplot
test_genes <- c(
  "IGHV1",    "IGHJ3", "IGHM",
  "IGKV4-53", "IGKJ2", "IGKC",
  "IGLV1",    "IGLJ1", "IGLC1"
)

arg_lst <- list(
  sobj_in     = list(tiny_vdj),
  gene_cols   = list("v_gene", "j_gene", c("v_gene", "j_gene")),
  cluster_col = NULL,
  chain       = list(NULL, "IGH", "IGL", "IGK"),
  chain_col   = "chains",
  type        = c("heatmap", "bar"),
  n_genes     = list(NULL, 10),
  plot_genes  = list(NULL, test_genes),
  plot_colors = list(NULL, test_cols),
  plot_lvls   = list(NULL, test_lvls),
  yaxis       = c("percent", "frequency")
)

test_all_args(
  arg_lst = arg_lst,
  .fn     = plot_usage,
  ttl     = "plot_usage args",
  chk     = expr(expect_s3_class(.res, "ggplot"))
)

# Check all plot_usage arguments return ggplot list
arg_lst$cluster_col <- "seurat_clusters"

test_all_args(
  arg_lst = arg_lst,
  .fn     = plot_usage,
  ttl     = "plot_usage args",
  chk     = expr(expect_type(.res, "list"))
)

# Check plot_usage bad type
test_that("plot_usage bad type", {
  expect_error(
    plot_usage(
      sobj_in   = tiny_vdj,
      gene_cols = "v_gene",
      type      = "BAD"
    )
  )
})

# Check plot_usage bad yaxis
test_that("plot_usage bad type", {
  expect_error(
    plot_usage(
      sobj_in   = tiny_vdj,
      gene_cols = "v_gene",
      yaxis     = "BAD"
    )
  )
})

# Check plot_usage bad gene_cols
test_that("plot_usage bad gene_cols", {
  expect_error(
    plot_usage(
      sobj_in   = tiny_vdj,
      gene_cols = c("v_gene", "d_gene", "j_gene")
    )
  )
})

# Check .set_lvls levels
test_that(".set_lvls args", {
  res <- .set_lvls(
    df_in = tiny_dat,
    clmn  = "seurat_clusters",
    lvls  = test_lvls
  )

  expect_identical(test_lvls, levels(res$seurat_clusters))
})

# Check .set_lvls bad levels
test_that(".set_lvls bad lvls", {
  lvls <- test_lvls[2:length(test_lvls)]

  expect_error(
    .set_lvls(
      df_in = tiny_dat,
      clmn  = "seurat_clusters",
      lvls  = lvls
    )
  )
})
