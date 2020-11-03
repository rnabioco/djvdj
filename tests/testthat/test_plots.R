
test_cols <- c(
  "#E69F00", "#56B4E9", "#009E73",
  "#F0E442", "#d7301f", "#0072B2",
  "#D55E00", "#6A51A3", "#CC79A7",
  "#999999", "#875C04"
)

test_lvls <- unique(tiny_vdj$seurat_clusters)

test_that("Check plot_features", {
  args_df <- list(
    feature     = c("seurat_clusters", "nCount_RNA"),
    data_slot   = c("data", "counts"),
    plot_colors = list(NULL, test_cols),
    feat_lvls   = list(NULL, test_lvls),
    facet_id    = list(NULL, "orig.ident"),
    facet_lvls  = list(NULL, "orig.ident"),
    min_pct     = list(NULL, 0.05),
    max_pct     = list(NULL, 0.95),
    lm_line     = c(TRUE, FALSE)
  ) %>%
    expand.grid() %>%
    mutate(
      across(where(is.factor), as.character),
      n = rownames(.)
    )

  res <- pmap(
    args_df,
    check_args,
    .fn     = plot_features,
    obj_in  = tiny_vdj
  )

  walk(res, expect_s3_class, "ggplot")
})

test_that("Check plot_cell_count", {
  args_df <- list(
    x           = "orig.ident",
    fill_col    = list(NULL, "seurat_clusters"),
    facet_col   = list(NULL, "orig.ident"),
    yaxis       = c("fraction", "count"),
    plot_colors = list(NULL, test_cols),
    plot_lvls   = "orig.ident",
    order_count = c(TRUE, FALSE),
    n_label     = c(TRUE, FALSE),
    label_aes   = list(list(), list(size = 2))
  ) %>%
    expand.grid() %>%
    mutate(
      across(where(is.factor), as.character),
      n = rownames(.)
    )

  res <- pmap(
    args_df,
    check_args,
    .fn     = plot_cell_count,
    sobj_in = tiny_vdj
  )

  walk(res, expect_s3_class, "ggplot")
})

test_that("Check plot_reads", {
  args_df <- list(
    data_cols   = list("reads", "umis", c("reads", "umis")),
    chain_col   = list(NULL, "chains"),
    cluster_col = list(NULL, "seurat_clusters"),
    type        = c("violin", "histogram", "density"),
    plot_colors = list(NULL, test_cols),
    plot_lvls   = list(NULL, test_lvls)
  ) %>%
    expand.grid() %>%
    mutate(n = rownames(.))

  res <- pmap(
    args_df,
    check_args,
    .fn     = plot_reads,
    sobj_in = tiny_vdj
  )

  walk(res, expect_s3_class, "ggplot")
})

test_that("Check plot_abundance line plot", {
  args_df <- list(
    cluster_col = list(NULL, "seurat_clusters"),
    label_col   = list(NULL, "cdr3"),
    yaxis       = c("percent", "frequency"),
    plot_colors = list(NULL, test_cols),
    plot_lvls   = list(NULL, test_lvls),
    label_aes   = list(list(), list(size = 2))
  ) %>%
    expand.grid() %>%
    mutate(n = rownames(.))

  res <- pmap(
    args_df,
    check_args,
    .fn           = plot_abundance,
    sobj_in       = tiny_vdj,
    clonotype_col = "cdr3",
    type          = "line"
  )

  walk(res, expect_s3_class, "ggplot")
})

test_that("Check plot_abundance bar plot", {
  args_df <- list(
    cluster_col = list(NULL, "seurat_clusters"),
    yaxis       = c("percent", "frequency"),
    plot_colors = list(NULL, test_cols),
    plot_lvls   = list(NULL, test_lvls),
    label_aes   = list(list(), list(size = 2))
  ) %>%
    expand.grid() %>%
    mutate(n = rownames(.))

  res <- pmap(
    args_df,
    check_args,
    .fn           = plot_abundance,
    sobj_in       = tiny_vdj,
    clonotype_col = "cdr3",
    type          = "bar",
    label_col     = "cdr3"
  )

  walk(res, expect_s3_class, "ggplot")
})

test_that("Check plot_diversity", {
  mets <- abdiv::alpha_diversities %>%
    map(~ eval(parse(text = paste0("abdiv::", .x))))

  args_df <- list(
    cluster_col = list(NULL, "seurat_clusters"),
    method      = mets,
    plot_colors = list(NULL, test_cols),
    plot_lvls   = list(NULL, test_lvls)
  ) %>%
    expand.grid() %>%
    mutate(n = rownames(.))

  res <- pmap(
    args_df,
    check_args,
    .fn           = plot_diversity,
    sobj_in       = tiny_vdj,
    clonotype_col = "cdr3"
  )

  walk(res, expect_s3_class, "ggplot")
})

test_that("Check plot_similarity", {
  mets <- abdiv::beta_diversities %>%
    map(~ eval(parse(text = paste0("abdiv::", .x))))

  args_df <- list(
    method      = mets,
    plot_colors = list(NULL, test_cols)
  ) %>%
    expand.grid() %>%
    mutate(n = rownames(.))

  res <- pmap(
    args_df,
    check_args,
    .fn           = plot_similarity,
    sobj_in       = tiny_vdj,
    clonotype_col = "cdr3",
    cluster_col   = "seurat_clusters"
  )

  walk(res, expect_s3_class, "ggplot")
})

test_that("Check plot_usage", {
  args_df <- list(
    gene_cols   = list("v_gene", "d_gene", "j_gene", "c_gene", c("v_gene", "j_gene")),
    cluster_col = list(NULL, "seurat_clusters"),
    chain       = list(NULL, "IGH", "IGL", "IGK"),
    type        = c("heatmap", "bar"),
    plot_colors = list(NULL, test_cols),
    plot_lvls   = list(NULL, test_lvls),
    yaxis       = c("percent", "frequency")
  ) %>%
    expand.grid() %>%
    mutate(n = rownames(.))

  res <- pmap(
    args_df,
    check_args,
    .fn       = plot_usage,
    sobj_in   = tiny_vdj,
    chain_col = "chains"
  )
})
