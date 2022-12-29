# Environment that holds global variables
global <- new.env()

# Column name to use for storing cell barcodes
global$cell_col <- ".cell_id"

# Default chain column
global$chain_col <- "chains"

# Default clonotype column
global$clonotype_col <- "clonotype_id"

# Default chain separator
global$sep <- ";"

# Default argument classes to use with .check_args()
global$arg_classes <- list(
  data_col      = list(),
  cluster_col   = list(allow_null = TRUE),
  group_col     = list(allow_null = TRUE),
  clonotype_col = list(allow_null = TRUE),
  chain_col     = list(),
  downsample    = list(Class = "logical"),
  n_boots       = list(Class = "numeric"),
  chain         = list(len_one = FALSE, allow_null = TRUE),
  prefix        = list(allow_null = TRUE),
  return_df     = list(Class = "logical"),
  sep           = list(allow_null = TRUE),
  plot_colors   = list(len_one = FALSE, allow_null = TRUE),
  plot_lvls     = list(
    Class = list(c("character", "factor")), len_one = FALSE, allow_null = TRUE
  ),
  panel_nrow    = list(Class = "numeric", allow_null = TRUE),
  panel_scales  = list(),
  n_label       = list(Class = "logical"),
  label_params  = list(Class = "list", len_one = FALSE),
  units         = list(),
  trans         = list(),
  per_cell      = list(Class = "logical"),

  # plot_rarefaction
  ci_alpha = list(Class = "numeric"),

  # plot_clone_frequency
  n_clones = list(Class = "numeric", allow_null = TRUE),

  # plot_frequency
  n_top       = list(Class = "numeric", allow_null = TRUE),
  other_label = list(),
  stack       = list(Class = "logical"),

  # plot_gene_usage
  vdj_genes     = list(len_one = FALSE, allow_null = TRUE),
  n_genes       = list(Class = "numeric"),
  rotate_labels = list(Class = "logical"),
  return_list   = list(Class = "logical"),

  # calc_similarity
  return_mat = list(Class = "logical"),

  # plot_similarity
  cluster_heatmap       = list(Class = "logical"),
  remove_upper_triangle = list(Class = "logical"),
  remove_diagonal       = list(Class = "logical"),

  # cluster_sequences
  resolution  = list(Class = "numeric", len_one = FALSE),
  k           = list(Class = "numeric"),
  dist_method = list(allow_null = TRUE),
  run_umap    = list(Class = "logical"),

  # plot_motifs
  width     = list(Class = "numeric"),
  align_end = list(),

  # import_vdj
  vdj_dir           = list(len_one = FALSE, allow_null = TRUE),
  filter_paired     = list(Class = "logical"),
  define_clonotypes = list(allow_null = TRUE),
  include_mutations = list(Class = "logical"),
  aggr_dir          = list(allow_null = TRUE),

  # fetch_vdj
  filtrer_cells = list(Class = "logical"),
  unnest        = list(Class = "logical"),

  # summarize_vdj
  col_names = list(allow_null = TRUE),

  # plot_features
  feature   = list(allow_null = TRUE),
  x         = list(),
  y         = list(),
  min_q     = list(Class = "numeric", allow_null = TRUE),
  max_q     = list(Class = "numeric", allow_null = TRUE),
  na_color  = list(),
  data_slot = list(),

  # Arguments that vary
  data_cols     = list(len_one = FALSE),
  method        = list(),
  filter_chains = list(Class = "logical")
)
