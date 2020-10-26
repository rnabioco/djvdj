
# djvdj <img src="man/figures/djvdj-logo.png" align="right" height="145">

<!-- badges: start -->

[![R build
status](https://github.com/rnabioco/djvdj/workflows/R-CMD-check/badge.svg)](https://github.com/rnabioco/djvdj/actions)
<!-- badges: end -->

The goal of djvdj is to provide tools to analyze AVID-seq signals
alongside single-cell VDJ sequencing data.

<br>

## Installation

You can install the development version of djvdj from
[GitHub](https://github.com/rnabioco/djvdj) with:

``` r
# install.packages("devtools")
devtools::install_github("rnabioco/djvdj")
```

<br>

## Vignette

Write a description here…

![](man/figures/README-rna_umap-1.png)<!-- -->

<br>

### Import VDJ data

`import_vdj` takes the output files from `cellranger vdj` and adds
clonotype information to the meta.data for an existing Seurat object.
For cells with multiple chains, the information for each chain is stored
as a single row, separated by a “;” (or a character specified by `sep`).
For cells that do not have any VDJ sequencing data, NAs will be added to
the meta.data.

If the Seurat object contains data for multiple runs, a vector
containing paths to the VDJ data for each sample can be given. If
multiple paths are provided, cell prefixes should be included as names
for the vector.

``` r
# Create vector of paths for cellranger output
samples <- levels(so_tcr$orig.ident)
paths   <- file.path("data", str_c(samples, "_TCR"))

names(paths) <- str_c(samples, "_GE")

# Import VDJ data
so_tcr <- import_vdj(
  sobj_in        = so_tcr,  # Seurat object
  vdj_dir        = paths,   # Directories containing cellranger output files
  prefix         = "",      # Prefix to add to new meta.data columns
  filter_contigs = TRUE     # Only include chains with at least one productive contig
)

# Take a look at the meta.data
vdj_cols <- c(
  "clonotype_id", "cdr3",
  "chains", "v_gene", 
  "j_gene", "reads",
  "umis"
)

so_tcr@meta.data %>%
  as_tibble() %>%
  filter(!is.na(clonotype_id)) %>%
  select(all_of(vdj_cols))
#> # A tibble: 3,850 x 7
#>    clonotype_id     cdr3             chains  v_gene      j_gene    reads   umis 
#>    <chr>            <chr>            <chr>   <chr>       <chr>     <chr>   <chr>
#>  1 WT_DN3_GE_clono… CAVQGANTEVFF     TRB     TRBV17      TRBJ1-1   42590   72   
#>  2 WT_DN3_GE_clono… CASSHPGQNSGNTLYF TRB     TRBV5       TRBJ1-3   12670   17   
#>  3 WT_DN3_GE_clono… CASSHWGETLYF     TRB     TRBV5       TRBJ2-3   31772   46   
#>  4 WT_DN3_GE_clono… CGARAQGLYNSPLYF  TRB     TRBV20      TRBJ1-6   7124    15   
#>  5 WT_DN3_GE_clono… CTAPAGGQNTEVFF   TRB     TRBV1       TRBJ1-1   4744    8    
#>  6 WT_DN3_GE_clono… CASSQDLDWGGEQFF  TRB     TRBV5       TRBJ2-1   34966   52   
#>  7 WT_DN3_GE_clono… CASRTGGCYEQYF;C… TRB;TRB TRBV13-1;T… TRBJ2-7;… 4166;6… 7;171
#>  8 WT_DN3_GE_clono… CASSPGTENTLYF    TRB     TRBV12-2    TRBJ2-4   9180    18   
#>  9 WT_DN3_GE_clono… CASSLKGARSDYTF   TRB     TRBV12-1    TRBJ1-2   19742   29   
#> 10 WT_DN3_GE_clono… CASRLTGRDSDYTF;… TRB;TRB TRBV15;TRB… TRBJ1-2;… 28860;… 48;22
#> # … with 3,840 more rows
```

<br>

### Quality Control

`djvdj` provides a range of tools to assess data quality and filter
single-cell V(D)J data. `mutate_vdj` can be used to add new cell labels
to the object meta.data, which is helpful for assessing quality.
`filter_vdj` will filter cells based on V(D)J information in the
meta.data or just remove V(D)J data from the object without discarding
cells. This can be used to filter V(D)J data based on read support. To
assess overall read and and UMI support, the function `plot_reads` will
summarize these metrics for each clonotype.

``` r
plot_reads(
  sobj_in      = so_tcr,        # Seurat object
  chain_col    = "chains",      # Column containing chains for each cell
  cluster_col  = "orig.ident",  # Column containing labels to group by
  plot_colors  = ito_cols       # Plot colors
) +
  guides(fill = FALSE, color = FALSE)
```

![](man/figures/README-read_support-1.png)<!-- -->

<br>

### Clonotype Abundance

To identify the top clonotypes in each sample or cluster, clonotype
abundance can be calculated using the `calc_abundance` function. The
corresponding `plot_abundance` function will plot clonotype frequency.

``` r
plot_abundance(
  sobj_in       = so_tcr,        # Seurat object
  clonotype_col = "cdr3",        # meta.data column containing clonotype IDs
  cluster_col   = "orig.ident",  # meta.data column containing cell labels
  
  plot_colors = ito_cols,        # Plot colors
  yaxis       = "percent",       # Units to plot
  label_col   = "cdr3",          # meta.data column containing labels
  n_labels    = 1,               # Number of top clonotypes to label
  size        = 1                # Additional ggplot options
) +
  theme(legend.title = element_blank())
```

![](man/figures/README-abund_plots-1.png)<!-- -->

<br>

### Repertoire Diversity

The function `calc_diversity` will calculate repertoire diversity on a
per-cluster basis. Using the `cluster_col` argument, any meta.data
column containing cell labels can be used for calculations.
`calc_diversity` uses the R package `abdiv` for performing diversity
calculations and any `abdiv` diversity function can be specified using
the `method` argument. The `plot_diversity` function can summarize these
results.

``` r
# Metrics to plot
fns <- list(
  "simpson"     = abdiv::simpson,
  "shannon"     = abdiv::shannon,
  "margalef"    = abdiv::margalef,
  "menhinick"   = abdiv::menhinick,
  "brillouin_d" = abdiv::brillouin_d
)

plot_diversity(
  sobj_in       = so_tcr,        # Seurat object
  clonotype_col = "cdr3",        # meta.data column containing clonotype ids
  cluster_col   = "orig.ident",  # meta.data column containing cell labels
  method        = fns,           # abdiv method to use
  plot_colors   = ito_cols
) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
```

![](man/figures/README-div_plots-1.png)<!-- -->

<br>

### Repertoire Overlap

To compare repertoires for different samples or clusters,
`calc_similarity` can calculate a variety of different similarity
metrics. The `cluster_col` should be used to specify the meta.data
column containing cell labels for comparison. Like `calc_diversity`, an
`abdiv` function can be specified with the `method` argument. These
results can be visualized with `plot_similarity`.

``` r
heat_theme <- theme(
  legend.title = element_blank(),
  legend.text  = element_text(size = 8)
)

# Sample heatmap
ident_heat <- plot_similarity(
  sobj_in       = so_tcr,                 # Seurat object
  clonotype_col = "cdr3",                 # meta.data column containing clonotype IDs
  cluster_col   = "orig.ident",           # meta.data column containing cell labels
  method        = abdiv::jaccard,         # Method to use
  plot_colors   = c("grey90", "#009E73")  # Plot colors
) +
  heat_theme

# Cluster heatmap
clust_heat <- plot_similarity(
  sobj_in       = so_tcr,
  clonotype_col = "cdr3",
  cluster_col   = "seurat_clusters",
  method        = abdiv::jaccard,
  plot_colors   = c("grey90", "#56B4E9"),  
  size          = 0.2,                    # Additional ggplot options
  color         = "white"                 # Additional ggplot options
) +
  heat_theme +
  theme(axis.text.x  = element_text(angle = 0))

# Combine heatmaps
plot_grid(ident_heat, clust_heat, align = "h")
```

![](man/figures/README-sim_plots-1.png)<!-- -->

<br>

### Gene Usage

The V(D)J data imported from Cell Ranger also includes the specific
genes detected for each cell. The function `calc_usage` can be used to
calculate the fraction of cells that express different V(D)J genes. This
function will produce a table summarizing the results. To only include
results for a certain chain, the `chain` and `chain_col` arguments can
be used to specify the meta.data column containing the chains detected
for each cell. These results can be visualized with the `plot_usage`
function.

``` r
plot_usage(
  sobj_in     = so_tcr,                # Seurat object
  gene_cols   = "v_gene",              # meta.data column(s) containing genes
  cluster_col = "orig.ident",          # meta.data column containing cell labels
  type        = "bar",                 # Type of plot
  chain       = "TRB",                 # Chain to use for filtering genes
  chain_col   = "chains",              # meta.data column containing chains
  
  yaxis       = "percent",             # Units to plot
  plot_colors = ito_cols,              # Colors to use for heatmap
  plot_genes  = NULL,                  # A list of genes to plot
  n_genes     = NULL,                  # The number of top genes to plot
  
  size        = 0.2,                   # Additional ggplot options
  color       = "white"                # Additional ggplot options
)
```

![](man/figures/README-usage_plots_1-1.png)<!-- -->

<br>

By passing multiple columns to `gene_cols`, the frequency that different
genes are used together can also be shown.

``` r
ggs <- plot_usage(
  sobj_in     = so_tcr,                 # Seurat object
  gene_cols   = c("v_gene", "j_gene"),  # meta.data column(s) containing genes
  cluster_col = "orig.ident",           # meta.data column containing cell labels
  chain       = "TRB",                  # Chain to use for filtering genes
  chain_col   = "chains",               # meta.data column containing chains identified
  plot_colors = c("grey90", "#6A51A3")  # Colors to use for heatmap
) %>%
  imap(~ .x + ggtitle(.y))

plot_grid(plotlist = ggs)
```

![](man/figures/README-usage_plots_2-1.png)<!-- -->
