# Test data
test_cols <- c(
  "#E69F00", "#56B4E9", "#009E73",
  "#F0E442", "#d7301f", "#0072B2",
  "#D55E00", "#6A51A3", "#CC79A7",
  "#999999", "#875C04", "#000000"
)

test_lvls <- unique(vdj_so$seurat_clusters) |>
  as.character() |>
  rev()

test_lvls_2 <- unique(test_so$sample) |>
  rev()

df_1 <- vdj_so@meta.data

df_2 <- vdj_so@meta.data |>
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
arg_lst$trans     <- c("identity", "log10")

test_all_args(
  arg_lst = arg_lst,
  .fn     = plot_features,
  desc    = "plot_features args num feat",
  chk     = expr(expect_s3_class(.res, "ggplot"))
)

# Check plot_features warning for numeric feature
test_that("plot_features warning num feat", {
  expect_warning(
    vdj_so |>
      plot_features(
        feature   = "nCount_RNA",
        plot_lvls = test_lvls
      )
  )
})

# Check plot_features error for same x and y
test_that("plot_features error same x y", {
  expect_error(
    vdj_so |>
      plot_features(
        x         = "UMAP_1",
        y         = "UMAP_1",
        feature   = "nCount_RNA",
        plot_lvls = test_lvls
      ),
    "`x` and `y` must be different"
  )
})

# Check plot_features feature not found
test_that("plot_fetures bad feature", {
  fn <- function(obj) {
    plot_features(obj, feature = "BAD_FEATURE")
  }

  expect_error(
    expect_warning(fn(vdj_so), "requested variables were not found"),
    "not present in input"
  )

  expect_error(fn(vdj_sce), "not present in input")
  expect_error(fn(df_1), "not present in input")
  expect_error(fn(df_2), "not present in input")
})

# Check all plot_vdj_feature arguments
arg_lst <- list(
  input       = list(vdj_so, vdj_sce, df_1),
  data_col    = "umis",
  chain       = list(NULL, "IGH", c("IGH", "IGK")),
  plot_colors = list(NULL, test_cols),
  min_q       = list(NULL, 0.05),
  max_q       = list(NULL, 0.95)
)

test_all_args(
  arg_lst = arg_lst,
  .fn     = plot_vdj_feature,
  desc    = "plot_vdj_feature args",
  chk     = expr(expect_s3_class(.res, "ggplot"))
)

arg_lst$data_col <- "chains"

chain_lvls <- vdj_so@meta.data |>
  pull(chains) |>
  unique() |>
  c("IGH;IGH", "IGH") |>
  list()

arg_lst$plot_lvls <- chain_lvls

test_all_args(
  arg_lst = arg_lst,
  .fn     = plot_vdj_feature,
  desc    = "plot_vdj_feature args",
  chk     = expr(expect_s3_class(.res, "ggplot"))
)

# plot_vdj_feature bad chain filtering
test_that("plot_vdj_feature bad chain filtering", {

  expect_warning(
    vdj_so |>
      plot_vdj_feature(
        data_col = "nCount_RNA",
        chain = "IGK"
      ),
    "Some columns do not contain per-chain"
  )
})

# Check all plot_vdj arguments
arg_lst <- list(
  input       = list(vdj_so, vdj_sce, df_1),
  data_col    = "umis",
  per_cell    = c(FALSE, TRUE),
  cluster_col = list(NULL, "seurat_clusters"),
  chain       = list(NULL, "IGH", c("IGH", "IGK")),
  method      = c("histogram", "density", "violin", "boxplot"),
  units       = c("frequency", "percent"),
  plot_colors = list(NULL, test_cols),
  trans       = c("identity", "log10")
)

test_all_args(
  arg_lst = arg_lst,
  .fn     = plot_vdj,
  desc    = "plot_vdj args",
  chk     = expr(expect_s3_class(.res, "ggplot"))
)

# Check all plot_clone_frequency arguments for line plot
arg_lst <- list(
  input        = list(vdj_so, vdj_sce),
  data_col     = "cdr3_nt",
  cluster_col  = list(NULL, "seurat_clusters"),
  method       = "line",
  units        = c("percent", "frequency"),
  plot_colors  = list(NULL, test_cols),
  plot_lvls    = list(NULL, test_lvls),
  label_params = list(list(), list(size = 2)),
  n_clones     = c(0, 5)
)

test_all_args(
  arg_lst = arg_lst,
  .fn     = plot_clone_frequency,
  desc    = "plot_clone_frequency line args",
  chk     = expr(expect_s3_class(.res, "ggplot"))
)

# Check all plot_clone_frequency arguments for bar plot
arg_lst$method <- "bar"
arg_lst$n_clones <- 5

test_all_args(
  arg_lst = arg_lst,
  .fn     = plot_clone_frequency,
  desc    = "plot_clone_frequency bar args",
  chk     = expr(expect_s3_class(.res, "ggplot"))
)

# Check plot_clone_frequency axis labels
arg_lst <- list(
  input  = list(vdj_so, vdj_sce),
  units  = "percent",
  method = c("bar", "line")
)

test_all_args(
  arg_lst = arg_lst,
  .fn     = plot_clone_frequency,
  desc    = "plot_clone_frequency axis labels",
  chk     = expr(expect_true(.res$labels$y == "% of cells"))
)

arg_lst$units <- "frequency"

test_all_args(
  arg_lst = arg_lst,
  .fn     = plot_clone_frequency,
  desc    = "plot_clone_frequency axis labels",
  chk     = expr(expect_true(.res$label$y == "number of cells"))
)

# Check plot_clone_frequency bad units
test_that("plot_clone_frequency bad units", {
  expect_error(
    vdj_so |>
      plot_clone_frequency(
        method   = "line",
        data_col = "cdr3_nt",
        units    = "BAD"
      )
  )
})

# Check plot_clone_frequency bad method
test_that("plot_clone_frequency bad method", {
  expect_error(
    vdj_so |>
      plot_clone_frequency(
        method          = "BAD",
        data_col = "cdr3_nt"
      )
  )
})

# Check plot_clone_frequency bad n_clonotypes
test_that("plot_clone_frequency bad n_clonotypes", {
  expect_error(
    vdj_so |>
      plot_clone_frequency(
        method          = "bar",
        data_col = "cdr3_nt",
        n_clones  = 0
      )
  )
})

# Check all plot_diversity arguments
mets <- abdiv::alpha_diversities |>
  map(~ eval(parse(text = paste0("abdiv::", .x))))

names(mets) <- abdiv::alpha_diversities

arg_lst <- list(
  input       = list(vdj_so, vdj_sce),
  data_col    = "cdr3_nt",
  cluster_col = list(NULL, "seurat_clusters"),
  method      = list(mets),
  chain       = list(NULL, "IGK"),
  plot_colors = list(NULL, test_cols),
  plot_lvls   = list(NULL, test_lvls)
)

test_all_args(
  arg_lst = arg_lst,
  .fn     = plot_diversity,
  desc    = "plot_diversity args",
  chk     = expr(expect_s3_class(.res, "ggplot"))
)

# Check plot_diversity bad method names
mets <- abdiv::alpha_diversities |>
  map(~ eval(parse(text = paste0("abdiv::", .x))))

test_that("plot_diversity bad names", {
  expect_error(
    vdj_so |>
      plot_diversity(
        data_col = "cdr3_nt",
        method   = unname(mets)
      ),
    "Names must be included"
  )
})

# Check all plot_rarefaction arguments
arg_lst <- list(
  input        = list(vdj_so, vdj_sce),
  data_col     = "cdr3_nt",
  cluster_col  = list(NULL, "seurat_clusters"),
  method       = c("richness", "shannon", "invsimpson"),
  chain        = list(NULL, "IGK"),
  plot_colors  = list(NULL, test_cols),
  plot_lvls    = list(NULL, test_lvls),
  panel_nrow   = list(NULL, 2),
  panel_scales = c("free", "fixed"),
  ci_alpha     = 0.5,
  n_boots      = 0
)

test_all_args(
  arg_lst = arg_lst,
  .fn     = plot_rarefaction,
  desc    = "plot_rarefaction args",
  chk     = expr(expect_s3_class(.res, "ggplot"))
)

# Check all plot_frequency arguments
arg_lst <- list(
  input       = list(test_so),
  data_col    = "cell_type",
  cluster_col = list(NULL, "sample"),
  units       = c("percent", "frequency"),
  stack       = c(TRUE, FALSE),
  plot_colors = list(NULL, test_cols),
  plot_lvls   = list(NULL, test_lvls_2),
  n_top       = list(NULL, 4),
  n_label     = c(TRUE, FALSE),
  label_params = list(list(size = 4), list())
)

test_all_args(
  arg_lst = arg_lst,
  .fn     = plot_frequency,
  desc    = "plot_frequency args",
  chk     = expr(expect_s3_class(.res, "ggplot"))
)

arg_lst$group_col   <- "orig.ident"
arg_lst$cluster_col <- "sample"
arg_lst$plot_lvls   <- list(NULL, c("BL6", "MD4"))

test_all_args(
  arg_lst = arg_lst,
  .fn     = plot_frequency,
  desc    = "plot_frequency args 2",
  chk     = expr(expect_s3_class(.res, "ggplot"))
)

# Check plot_diversity bad method names
mets <- abdiv::alpha_diversities |>
  map(~ eval(parse(text = paste0("abdiv::", .x))))

test_that("plot_diversity bad names", {
  expect_error(
    vdj_so |>
      plot_diversity(
        data_col = "cdr3_nt",
        method   = unname(mets)
      ),
    "Names must be included"
  )
})

# Check all plot_similarity arguments
# exclude methods that produce -Inf
exclude <- c(
  "morisita", "kullback_leibler_divergence",
  "kulczynski_first", "jaccard_turnover"
)

mets <- abdiv::beta_diversities
mets <- mets[!mets %in% exclude]

mets <- purrr::set_names(mets) |>
  map(~ eval(parse(text = paste0("abdiv::", .x))))

arg_lst <- list(
  input       = list(vdj_so, vdj_sce),
  data_col    = "cdr3_nt",
  cluster_col = "seurat_clusters",
  chain       = list(NULL, "IGH"),
  # method      = list(list("METHOD" = abdiv::jaccard), "count"),
  plot_colors = list(NULL, test_cols, "blue")
)

test_all_args(
  arg_lst = arg_lst,
  .fn     = plot_similarity,
  desc    = "plot_similarity args",
  chk     = expect_silent
)

# Check plot_similarity circos plot
arg_lst$chain <- NULL
arg_lst$method <- "circos"

test_all_args(
  arg_lst = arg_lst,
  .fn     = plot_similarity,
  desc    = "plot_similarity circos args",
  chk     = expect_silent
)

# Check plot_similarity bad data_col
test_that("plot_similarity bad clonotype col", {
  expect_error(
    vdj_so |>
      plot_similarity(
        data_col    = "BAD",
        cluster_col = "orig.ident"
      )
  )
})

# Check plot_similarity bad cluster_col
test_that("plot_similarity bad cluster col", {
  expect_error(
    vdj_so |>
      plot_similarity(
        cluster_col = NULL,
        data_col    = "cdr3_nt"
      )
  )
})

# Check calc_mds
arg_lst <- list(
  input       = list(test_so),
  data_col    = "cdr3_nt",
  cluster_col = "sample",
  method      = c("jaccard", "horn_morisita"),
  chain       = list(NULL, "IGK"),
  prefix      = c("", "TEST")
)

test_all_args(
  arg_lst = arg_lst,
  .fn     = calc_mds,
  desc    = "calc_mds Seurat args",
  chk     = expr(expect_s4_class(.res, "Seurat"))
)

# Check plot_mds
arg_lst <- list(
  input       = list(test_so),
  data_col    = "cdr3_nt",
  cluster_col = "sample",
  method      = c("jaccard", "horn_morisita"),
  chain       = list(NULL, "IGK"),
  plot_colors = list(NULL, c("red1", "red2", "red3", "blue1", "blue2", "blue3")),
  label_points = c(TRUE, FALSE),
  plot_lvls    = list(NULL, test_lvls_2)
)

test_all_args(
  arg_lst = arg_lst,
  .fn     = plot_mds,
  desc    = "plot_mds Seurat args",
  chk     = expr(expect_s3_class(.res, "ggplot"))
)

# Check all plot_gene_usage arguments for one gene
arg_lst <- list(
  input       = list(vdj_so, vdj_sce),
  data_cols   = list("v_gene", "j_gene"),
  chain       = list(NULL, "IGH", "IGL", "IGK"),
  cluster_col = list(NULL, "seurat_clusters"),
  method      = c("heatmap", "bar"),
  plot_colors = list(NULL, test_cols),
  plot_lvls   = list(NULL, test_lvls),
  units       = c("percent", "frequency"),
  trans       = c("identity", "log10")
)

test_all_args(
  arg_lst = arg_lst,
  .fn     = plot_gene_usage,
  desc    = "plot_gene_usage single gene args",
  chk     = expr(expect_s3_class(.res, "ggplot"))
)

# Check all plot_gene_usage arguments for two genes
# when circos option is included, the first circos test prints message which
# causes test to fail, this does not happen when running interactively, not
# sure what is going here
arg_lst$data_cols   <- list(c("v_gene", "j_gene"))
arg_lst$chain       <- list(NULL, "IGH")
arg_lst$method      <- "heatmap"
arg_lst$plot_colors <- list(NULL, rep(test_cols, 3))
arg_lst$plot_lvls   <- list(NULL, test_lvls)
arg_lst$trans       <- c("identity", "log1p")
arg_lst$return_list <- TRUE

test_all_args(
  arg_lst = arg_lst,
  .fn     = plot_gene_usage,
  desc    = "plot_gene_usage paired gene args",
  chk     = expect_silent
)

# Check plot_gene_usage vdj_genes single column
test_genes <- vdj_so |>
  fetch_vdj() |>
  pull(v_gene) |>
  na.omit() |>
  head(10)

arg_lst <- list(
  input       = list(vdj_so, vdj_sce),
  data_cols   = "v_gene",
  vdj_genes   = list(test_genes),
  method        = c("heatmap", "bar"),
  plot_colors = list(NULL, test_cols),
  plot_lvls   = list(NULL, test_lvls),
  units       = c("percent", "frequency"),
  trans       = c("identity", "log10")
)

test_all_args(
  arg_lst = arg_lst,
  .fn     = plot_gene_usage,
  desc    = "plot_gene_usage vdj_genes single args",
  chk     = expr(expect_s3_class(.res, "ggplot"))
)

# Check plot_gene_usage vdj_genes multiple columns
# when circos option is included, the first circos test fails on github actions
arg_lst$method    <- "heatmap"
arg_lst$data_cols <- list(c("v_gene", "j_gene"))
arg_lst$plot_lvls <- NULL
arg_lst$trans     <- c("identity", "log1p")

test_all_args(
  arg_lst = arg_lst,
  .fn     = plot_gene_usage,
  desc    = "plot_gene_usage vdj_genes paired args",
  chk     = expect_silent
)

# Check plot_gene_usage bad plot_genes
test_that("plot_gene_usage bad plot_genes", {
  expect_error(
    plot_gene_usage(vdj_so, data_cols = "v_gene", vdj_genes = "BAD"),
    "None of the provided genes were found"
  )

  expect_warning(
    plot_gene_usage(vdj_so, data_cols = "v_gene", vdj_genes = c(test_genes, "BAD")),
    "The following genes were not found"
  )
})

# Check plot_gene_usage bad method
test_that("plot_gene_usage bad method", {
  expect_error(
    vdj_so |>
      plot_gene_usage(
        data_cols = "v_gene",
        method      = "BAD"
      )
  )
})

# Check plot_gene_usage bad units
test_that("plot_gene_usage bad method", {
  expect_error(
    vdj_so |>
      plot_gene_usage(
        data_cols = "v_gene",
        units     = "BAD"
      )
  )
})

# Check plot_gene_usage bad data_cols
test_that("plot_gene_usage bad data_cols", {
  expect_error(
    vdj_so |>
      plot_gene_usage(data_cols = c("v_gene", "d_gene", "j_gene"))
  )
})

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

