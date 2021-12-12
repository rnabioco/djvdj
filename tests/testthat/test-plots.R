# Test data
test_cols <- c(
  "#E69F00", "#56B4E9", "#009E73",
  "#F0E442", "#d7301f", "#0072B2",
  "#D55E00", "#6A51A3", "#CC79A7",
  "#999999", "#875C04", "#000000"
)

test_lvls <- unique(vdj_so$seurat_clusters) %>%
  as.character() %>%
  rev()

df_1 <- vdj_so@meta.data

df_2 <- vdj_so@meta.data %>%
  as_tibble(rownames = ".cell_id")

# Check all plot_features arguments except data_slot
arg_lst <- list(
  x           = list("UMAP_1", c(x = "UMAP_1")),
  y           = list("UMAP_2", c(y = "UMAP_2")),
  input       = list(vdj_so, vdj_sce, df_1),
  feature     = list("seurat_clusters", c(clust = "seurat_clusters")),
  plot_colors = list(NULL, test_cols),
  plot_lvls   = list(NULL, test_lvls),
  min_q       = list(NULL, 0.05),
  max_q       = list(NULL, 0.95)
)

test_all_args(
  arg_lst = arg_lst,
  .fn     = plot_features,
  desc    = "plot_features args chr feat",
  chk     = expr(expect_s3_class(.res, "ggplot"))
)

# Check all plot_features arguments with data_slot
arg_lst$input     <- list(vdj_so)
arg_lst$data_slot <- c("data", "counts")

test_all_args(
  arg_lst = arg_lst,
  .fn     = plot_features,
  desc    = "plot_features args chr feat",
  chk     = expr(expect_s3_class(.res, "ggplot"))
)

# Check all plot_features arguments for numeric feature
arg_lst$feature   <- "nCount_RNA"
arg_lst$plot_lvls <- list(NULL)

test_all_args(
  arg_lst = arg_lst,
  .fn     = plot_features,
  desc    = "plot_features args num feat",
  chk     = expr(expect_s3_class(.res, "ggplot"))
)

# Check plot_features warning for numeric feature
test_that("plot_features warning num feat", {
  expect_warning(
    vdj_so %>%
      plot_features(
        feature   = "nCount_RNA",
        plot_lvls = test_lvls
      )
  )
})

# Check plot_features error for same x and y
test_that("plot_features error same x y", {
  expect_error(
    vdj_so %>%
      plot_features(
        x         = "UMAP_1",
        y         = "UMAP_1",
        feature   = "nCount_RNA",
        plot_lvls = test_lvls
      )
  )
})

# Check plot_features feature not found
arg_lst <- list(
  input   = list(vdj_so, vdj_sce, df_1),
  feature = "BAD_FEATURE"
)

test_all_args(
  arg_lst = arg_lst,
  .fn     = plot_features,
  desc    = "plot_features no feat found",
  chk     = expect_error
)

# Check all plot_reads arguments
arg_lst <- list(
  input       = list(vdj_so, vdj_sce),
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
  desc    = "plot_reads args",
  chk     = expr(expect_s3_class(.res, "ggplot"))
)

# Check plot_reads bad type
test_that("plot_reads bad type", {
  expect_error(
    vdj_so %>%
      plot_reads(type = "BAD")
  )
})

# Check all plot_abundance arguments for line plot
arg_lst <- list(
  input         = list(vdj_so, vdj_sce),
  clonotype_col = "cdr3_nt",
  cluster_col   = list(NULL, "seurat_clusters"),
  type          = "line",
  label_col     = c("cdr3", "orig.ident"),
  color_col     = list(NULL, "seurat_clusters"),
  yaxis         = c("percent", "frequency"),
  plot_colors   = list(NULL, test_cols),
  plot_lvls     = list(NULL, test_lvls),
  label_aes     = list(list(), list(size = 2)),
  n_clonotypes  = c(0, 5)
)

test_all_args(
  arg_lst = arg_lst,
  .fn     = plot_abundance,
  desc    = "plot_abundance line args",
  chk     = expr(expect_s3_class(.res, "ggplot"))
)

# Check all plot_abundance arguments for bar plot
arg_lst$type <- "bar"
arg_lst$n_clonotypes <- 5

test_all_args(
  arg_lst = arg_lst,
  .fn     = plot_abundance,
  desc    = "plot_abundance bar args",
  chk     = expr(expect_s3_class(.res, "ggplot"))
)

# Check plot_abundance axis labels
arg_lst <- list(
  input         = list(vdj_so, vdj_sce),
  clonotype_col = "cdr3_nt",
  label_col     = "cdr3",
  yaxis         = "percent",
  type          = c("bar", "line")
)

test_all_args(
  arg_lst = arg_lst,
  .fn     = plot_abundance,
  desc    = "plot_abundance axis labels",
  chk     = expr(expect_true(.res$label$y == "percent"))
)

arg_lst$yaxis <- "frequency"

test_all_args(
  arg_lst = arg_lst,
  .fn     = plot_abundance,
  desc    = "plot_abundance axis labels",
  chk     = expr(expect_true(.res$label$y == "frequency"))
)

# Check plot_abundance bad yaxis
test_that("plot_abundance bad yaxis", {
  expect_error(
    vdj_so %>%
      plot_abundance(
        type          = "line",
        clonotype_col = "cdr3_nt",
        yaxis         = "BAD"
      )
  )
})

# Check plot_abundance bad type
test_that("plot_abundance bad type", {
  expect_error(
    vdj_so %>%
      plot_abundance(
        type          = "BAD",
        clonotype_col = "cdr3_nt"
      )
  )
})

# Check plot_abundance bad n_clonotypes
test_that("plot_abundance bad n_clonotypes", {
  expect_error(
    vdj_so %>%
      plot_abundance(
        type          = "bar",
        clonotype_col = "cdr3_nt",
        n_clonotypes  = 0
      )
  )
})

# Check plot_abundance bad label_col
test_that("plot_abundance bad label_col", {
  expect_error(
    vdj_so %>%
      plot_abundance(
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
  input         = list(vdj_so, vdj_sce),
  clonotype_col = "cdr3_nt",
  cluster_col   = list(NULL, "seurat_clusters"),
  method        = append(mets, list(mets)),
  plot_colors   = list(NULL, test_cols),
  plot_lvls     = list(NULL, test_lvls)
)

test_all_args(
  arg_lst = arg_lst,
  .fn     = plot_diversity,
  desc    = "plot_diversity args",
  chk     = expr(expect_s3_class(.res, "ggplot"))
)

# Check plot_diversity bad method names
mets <- abdiv::alpha_diversities %>%
  map(~ eval(parse(text = paste0("abdiv::", .x))))

test_that("plot_diversity bad names", {
  expect_error(
    vdj_so %>%
      plot_diversity(
        clonotype_col = "cdr3_nt",
        method        = mets
      )
  )
})

# Check all plot_similarity arguments
mets <- abdiv::beta_diversities %>%
  map(~ eval(parse(text = paste0("abdiv::", .x))))

arg_lst <- list(
  input         = list(vdj_so, vdj_sce),
  clonotype_col = "cdr3_nt",
  cluster_col   = "seurat_clusters",
  method        = mets,
  plot_colors   = list(NULL, test_cols)
)

test_all_args(
  arg_lst = arg_lst,
  .fn     = plot_similarity,
  desc    = "plot_similarity args",
  chk     = expr(expect_s3_class(.res, "ggplot"))
)

# Check plot_similarity bad clonotype_col
test_that("plot_similarity bad clonotype col", {
  expect_error(
    vdj_so %>%
      plot_similarity(
        clonotype_col = NULL,
        cluster_col   = "orig.ident"
      )
  )
})

# Check plot_similarity bad cluster_col
test_that("plot_similarity bad cluster col", {
  expect_error(
    vdj_so %>%
      plot_similarity(
        cluster_col   = NULL,
        clonotype_col = "cdr3_nt"
      )
  )
})

# Check all plot_usage arguments return ggplot
arg_lst <- list(
  input       = list(vdj_so, vdj_sce),
  gene_cols   = list("v_gene", "j_gene", c("v_gene", "j_gene")),
  chain       = list(NULL, "IGH", "IGL", "IGK"),
  type        = c("heatmap", "bar"),
  plot_colors = list(NULL, test_cols),
  plot_lvls   = list(NULL, test_lvls),
  yaxis       = c("percent", "frequency")
)

test_all_args(
  arg_lst = arg_lst,
  .fn     = plot_usage,
  desc    = "plot_usage args",
  chk     = expr(expect_s3_class(.res, "ggplot"))
)

# Check all plot_usage arguments return ggplot list
arg_lst$cluster_col <- "seurat_clusters"

test_all_args(
  arg_lst = arg_lst,
  .fn     = plot_usage,
  desc    = "plot_usage args",
  chk     = expr(expect_type(.res, "list"))
)

# Check plot_usage plot_genes
test_genes <- c(
  "IGHV1-26", "IGHV1-72", "IGHV1-64",
  "IGHV3-6",  "IGHV1-82", "IGHV1-53"
)

arg_lst <- list(
  input       = list(vdj_so, vdj_sce),
  gene_cols   = list("v_gene", c("v_gene", "j_gene")),
  plot_genes  = list(test_genes),
  type        = c("heatmap", "bar"),
  plot_colors = list(NULL, test_cols),
  plot_lvls   = list(NULL, test_lvls),
  yaxis       = c("percent", "frequency")
)

test_all_args(
  arg_lst = arg_lst,
  .fn     = plot_usage,
  desc    = "plot_usage args",
  chk     = expr(expect_s3_class(.res, "ggplot"))
)

# Check plot_usage bad plot_genes
test_that("plot_usage bad plot_genes", {
  expect_error(
    plot_usage(vdj_so, gene_cols = "v_gene", plot_genes = "BAD"),
    "None of the provided genes were found"
  )

  expect_warning(
    plot_usage(vdj_so, gene_cols = "v_gene", plot_genes = c(test_genes, "BAD")),
    "Some genes not found: "
  )
})

# Check plot_usage bad type
test_that("plot_usage bad type", {
  expect_error(
    vdj_so %>%
      plot_usage(
        gene_cols = "v_gene",
        type      = "BAD"
      )
  )
})

# Check plot_usage bad yaxis
test_that("plot_usage bad type", {
  expect_error(
    vdj_so %>%
      plot_usage(
        gene_cols = "v_gene",
        yaxis     = "BAD"
      )
  )
})

# Check plot_usage bad gene_cols
test_that("plot_usage bad gene_cols", {
  expect_error(
    vdj_so %>%
      plot_usage(gene_cols = c("v_gene", "d_gene", "j_gene"))
  )
})

# Check all plot_cdr3_length arguments
arg_lst <- list(
  input       = list(vdj_so, vdj_sce, df_1),
  length_col  = c("cdr3_length", "cdr3_nt_length"),
  cluster_col = list(NULL, "seurat_clusters"),
  chain       = list(NULL, "IGH", c("IGH", "IGK")),
  type        = c("histogram", "density", "violin"),
  yaxis       = c("frequency", "percent"),
  plot_colors = list(NULL, test_cols),
  plot_lvls   = list(NULL, test_lvls)
)

test_all_args(
  arg_lst = arg_lst,
  .fn     = plot_cdr3_length,
  desc    = "plot_cdr3_length args",
  chk     = expr(expect_s3_class(.res, "ggplot"))
)

# Check .set_lvls levels
arg_lst <- list(df_in = list(df_1, df_2))

pwalk(arg_lst, ~ {
  test_that(".set_lvls args", {
    res <- .set_lvls(
      ...,
      clmn  = "seurat_clusters",
      lvls  = test_lvls
    )

    expect_identical(test_lvls, levels(res$seurat_clusters))
  })
})

# Check .set_lvls bad levels
lvls <- test_lvls[2:length(test_lvls)]

arg_lst <- list(
  df_in = list(df_1, df_2),
  clmn  = "seurat_clusters",
  lvls  = lvls
)

test_all_args(
  arg_lst = arg_lst,
  .fn     = .set_lvls,
  desc    = ".set_lvls bad lvls",
  chk     = expect_error
)

