# Environment that holds global variables
djvdj_global <- new.env(parent = emptyenv())

# Column name to use for storing cell barcodes
djvdj_global$cell_col <- ".cell_id"

# Default argument classes to use with .check_args()
djvdj_global$arg_classes <- list(
  list(arg = "data_col"),
  list(arg = "cluster_col", allow_null = TRUE),
  list(arg = "group_col", allow_null = TRUE),
  list(arg = "clonotype_col", allow_null = TRUE),
  list(arg = "chain_col"),
  list(arg = "downsample", Class = "logical"),
  list(arg = "n_boots", Class = "numeric"),
  list(arg = "chain", len_one = FALSE, allow_null = TRUE),
  list(arg = "prefix", allow_null = TRUE),
  list(arg = "return_df", Class = "logical"),
  list(arg = "sep", allow_null = TRUE),
  list(arg = "plot_colors", len_one = FALSE, allow_null = TRUE),
  list(arg = "plot_lvls", Class = list(c("character", "factor")), len_one = FALSE, allow_null = TRUE),
  list(arg = "panel_nrow", Class = "numeric", allow_null = TRUE),
  list(arg = "panel_scales"),
  list(arg = "n_label", Class = "logical"),
  list(arg = "label_params", Class = "list", len_one = FALSE),
  list(arg = "units"),
  list(arg = "trans"),
  list(arg = "per_cell", Class = "logical"),

  # plot_rarefaction
  list(arg = "ci_alpha", Class = "numeric"),

  # plot_clone_frequency
  list(arg = "n_clones", Class = "numeric", allow_null = TRUE),

  # plot_frequency
  list(arg = "n_top", Class = "numeric", allow_null = TRUE),
  list(arg = "other_label"),
  list(arg = "stack", Class = "logical"),

  # plot_gene_usage
  list(arg = "vdj_genes", len_one = FALSE, allow_null = TRUE),
  list(arg = "n_genes", Class = "numeric"),
  list(arg = "rotate_labels", Class = "logical"),
  list(arg = "return_list", Class = "logical"),

  # calc_similarity
  list(arg = "return_mat", Class = "logical"),

  # plot_similarity
  list(arg = "cluster_heatmap", Class = "logical"),
  list(arg = "remove_upper_triangle", Class = "logical"),
  list(arg = "remove_diagonal", Class = "logical"),

  # cluster_sequences
  list(arg = "resolution", Class = "numeric", len_one = FALSE),
  list(arg = "k", Class = "numeric"),
  list(arg = "dist_method", allow_null = TRUE),
  list(arg = "run_umap", Class = "logical"),

  # plot_motifs
  list(arg = "width", Class = "numeric"),
  list(arg = "align_end"),

  # import_vdj
  list(arg = "vdj_dir", len_one = FALSE, allow_null = TRUE),
  list(arg = "filter_paired", Class = "logical"),
  list(arg = "define_clonotypes", allow_null = TRUE),
  list(arg = "include_mutations", Class = "logical"),
  list(arg = "aggr_dir", allow_null = TRUE),

  # fetch_vdj
  list(arg = "filter_cells", Class = "logical"),
  list(arg = "unnest", Class = "logical"),

  # summarize_vdj
  list(arg = "col_names", allow_null = TRUE),

  # plot_features
  list(arg = "feature", allow_null = TRUE),
  list(arg = "x"),
  list(arg = "y"),
  list(arg = "min_q", Class = "numeric", allow_null = TRUE),
  list(arg = "max_q", Class = "numeric", allow_null = TRUE),
  list(arg = "na_color"),
  list(arg = "data_slot"),

  # Arguments that vary
  data_cols     = list(arg = "data_cols", len_one = FALSE),
  method        = list(arg = "method"),
  filter_chains = list(arg = "filter_chains", Class = "logical")
)
