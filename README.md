
<!-- README.md is generated from README.Rmd. Please edit that file -->

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

## TCR analysis

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
paths   <- file.path("data", str_c(samples, "_TCR"), "outs")

names(paths) <- str_c(samples, "_GE")

# Import VDJ data
so_tcr <- import_vdj(
  sobj_in         = so_tcr,  # Seurat object
  vdj_dir         = paths,   # Directories containing cellranger output files
  prefix          = "",      # Prefix to add to new meta.data columns
  filter_contigs = TRUE      # Only include chains with at least one productive contig
)

# Take a look at the added meta.data
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
#>    clonotype_id cdr3               chains  v_gene      j_gene     reads    umis 
#>    <chr>        <chr>              <chr>   <chr>       <chr>      <chr>    <chr>
#>  1 clonotype5   CAVQGANTEVFF       TRB     TRBV17      TRBJ1-1    42590    72   
#>  2 clonotype6   CASSHPGQNSGNTLYF   TRB     TRBV5       TRBJ1-3    12670    17   
#>  3 clonotype7   CASSHWGETLYF       TRB     TRBV5       TRBJ2-3    31772    46   
#>  4 clonotype8   CGARAQGLYNSPLYF    TRB     TRBV20      TRBJ1-6    7124     15   
#>  5 clonotype9   CTAPAGGQNTEVFF     TRB     TRBV1       TRBJ1-1    4744     8    
#>  6 clonotype10  CASSQDLDWGGEQFF    TRB     TRBV5       TRBJ2-1    34966    52   
#>  7 clonotype11  CASRTGGCYEQYF;CAS… TRB;TRB TRBV13-1;T… TRBJ2-7;T… 4166;67… 7;171
#>  8 clonotype2   CASSPGTENTLYF      TRB     TRBV12-2    TRBJ2-4    9180     18   
#>  9 clonotype12  CASSLKGARSDYTF     TRB     TRBV12-1    TRBJ1-2    19742    29   
#> 10 clonotype13  CASRLTGRDSDYTF;CA… TRB;TRB TRBV15;TRB… TRBJ1-2;T… 28860;1… 48;22
#> # … with 3,840 more rows
```

<br>

### Filtering

`filter_vdj` allows you to filter a Seurat object using the added
clonotype information or any other columns present in the meta.data.
When filtering, columns with VDJ data will be expanded based on the
delimiter “;” (or a character passed to `sep`). The columns that are
expanded for filtering can be specified with the `vdj_cols` argument. By
default filtering is only performed on cells that include VDJ data.

Filter to only include cells with paired alpha and beta chains.

``` r
so_filt <- filter_vdj(
  sobj_in  = so_tcr,                            # Seurat object
  filt     = all(c("TRA", "TRB") %in% chains),  # Expression for filtering
  vdj_cols = "chains"
)

# Take a look at the meta.data
so_filt@meta.data %>%
  as_tibble() %>%
  filter(!is.na(clonotype_id)) %>%
  select(all_of(vdj_cols))
#> # A tibble: 654 x 7
#>    clonotype_id cdr3             chains   v_gene       j_gene     reads    umis 
#>    <chr>        <chr>            <chr>    <chr>        <chr>      <chr>    <chr>
#>  1 clonotype490 CALSGSNTGYQNFYF… TRA;TRB  TRAV15-2-DV… TRAJ49;TR… 2286;13… 3;26 
#>  2 clonotype7   CAASASANKMIF;CA… TRA;TRB  TRAV14-2;TR… TRAJ47;TR… 19138;2… 13;23
#>  3 clonotype10  CATTGFASALTF;CA… TRA;TRB… TRAV16D-DV1… TRAJ35;TR… 904;639… 2;6;…
#>  4 clonotype11  CALGMNYNQGKLIF;… TRA;TRB  TRAV6N-7;TR… TRAJ23;TR… 8336;84… 8;8  
#>  5 clonotype16  CALISGSFNKLTF;C… TRA;TRB  TRAV12D-2;T… TRAJ4;TRB… 2190;82… 1;10 
#>  6 clonotype25  CAMRGGEGSWQLIF;… TRA;TRB  TRAV16N;TRB… TRAJ22;TR… 11096;1… 8;11 
#>  7 clonotype26  CAAYNYAQGLTF;CA… TRA;TRB  TRAV14N-3;T… TRAJ26;TR… 10850;1… 7;15 
#>  8 clonotype27  CALVMNYNQGKLIF;… TRA;TRB… TRAV13-1;TR… TRAJ23;TR… 16886;1… 13;1…
#>  9 clonotype29  CAVSNNNNAPRF;CA… TRA;TRB  TRAV3-3;TRB… TRAJ43;TR… 11874;2… 8;25 
#> 10 clonotype30  CAGHYNVLYF;CASS… TRA;TRB  TRAV7-3;TRB… TRAJ21;TR… 4144;59… 2;6  
#> # … with 644 more rows
```

<br>

Instead of filtering, `filter_vdj` can also add a new column to the
object meta.data using the `new_col` argument. The values present in
`new_col` will be based on the values (or expression) passed to the
`true` and `false` arguments.

In this example a new column is added indicating whether the cell has a
paired alpha and beta chain. This is useful for generating new cell
labels for plotting.

``` r
so_tcr <- filter_vdj(
  sobj_in = so_tcr,                            # Seurat object
  filt    = all(c("TRA", "TRB") %in% chains),  # Condition to use for filtering
  new_col = "Paired",                          # Name of new column
  true    = "paired",                          # Value when condition is TRUE
  false   = "unpaired"                         # Value when condition is FALSE
)

# Take a look at the meta.data
vdj_cols <- c(vdj_cols, "Paired")

so_tcr@meta.data %>%
  as_tibble() %>%
  filter(!is.na(clonotype_id)) %>%
  select(all_of(vdj_cols))
#> # A tibble: 3,850 x 8
#>    clonotype_id cdr3           chains  v_gene     j_gene    reads   umis  Paired
#>    <chr>        <chr>          <chr>   <chr>      <chr>     <chr>   <chr> <chr> 
#>  1 clonotype5   CAVQGANTEVFF   TRB     TRBV17     TRBJ1-1   42590   72    unpai…
#>  2 clonotype6   CASSHPGQNSGNT… TRB     TRBV5      TRBJ1-3   12670   17    unpai…
#>  3 clonotype7   CASSHWGETLYF   TRB     TRBV5      TRBJ2-3   31772   46    unpai…
#>  4 clonotype8   CGARAQGLYNSPL… TRB     TRBV20     TRBJ1-6   7124    15    unpai…
#>  5 clonotype9   CTAPAGGQNTEVFF TRB     TRBV1      TRBJ1-1   4744    8     unpai…
#>  6 clonotype10  CASSQDLDWGGEQ… TRB     TRBV5      TRBJ2-1   34966   52    unpai…
#>  7 clonotype11  CASRTGGCYEQYF… TRB;TRB TRBV13-1;… TRBJ2-7;… 4166;6… 7;171 unpai…
#>  8 clonotype2   CASSPGTENTLYF  TRB     TRBV12-2   TRBJ2-4   9180    18    unpai…
#>  9 clonotype12  CASSLKGARSDYTF TRB     TRBV12-1   TRBJ1-2   19742   29    unpai…
#> 10 clonotype13  CASRLTGRDSDYT… TRB;TRB TRBV15;TR… TRBJ1-2;… 28860;… 48;22 unpai…
#> # … with 3,840 more rows
```

![](man/figures/README-pair_umap-1.png)<!-- -->

<br>

More complicated statements referring to meta.data columns can be used
for the `filt`, `true`, and `false` arguments. For more detailed
analysis of the chains detected for each cell, a new cell label can be
created for only the unique chains detected in each cell.

``` r
so_tcr <- filter_vdj(
  sobj_in = so_tcr,                                # Seurat object
  new_col = "uniq_chains",                         # Name of new column
  true    = str_c(unique(chains), collapse = "_")  # Value when condition is TRUE
)

# Take a look at the meta.data
vdj_cols <- c(vdj_cols, "uniq_chains")

so_tcr@meta.data %>%
  as_tibble() %>%
  filter(!is.na(clonotype_id), n_chains > 2) %>%
  select(all_of(vdj_cols))
#> # A tibble: 101 x 9
#>    clonotype_id cdr3     chains  v_gene   j_gene  reads umis  Paired uniq_chains
#>    <chr>        <chr>    <chr>   <chr>    <chr>   <chr> <chr> <chr>  <chr>      
#>  1 clonotype10  CATTGFA… TRA;TR… TRAV16D… TRAJ35… 904;… 2;6;… paired TRA_TRB    
#>  2 clonotype27  CALVMNY… TRA;TR… TRAV13-… TRAJ23… 1688… 13;1… paired TRA_TRB    
#>  3 clonotype36  CAAYSGG… TRA;TR… TRAV14-… TRAJ53… 1210… 14;8… paired TRA_TRB    
#>  4 clonotype85  CVVVDLP… TRA;TR… TRAV11;… TRAJ28… 6974… 7;11… paired TRA_TRB    
#>  5 clonotype136 CAARRGS… TRA;TR… TRAV14D… TRAJ33… 2014… 13;3… paired TRA_TRB    
#>  6 clonotype138 CVAHNNA… TRA;TR… TRAV11D… TRAJ39… 7230… 5;10… paired TRA_TRB    
#>  7 clonotype174 CAASTSG… TRA;TR… TRAV7D-… TRAJ22… 1815… 11;8… paired TRA_TRB    
#>  8 clonotype205 CALLASS… TRA;TR… TRAV12-… TRAJ50… 1294… 19;5… paired TRA_TRB    
#>  9 clonotype215 CAPGTGG… TRA;TR… TRAV13-… TRAJ12… 1188… 13;4… paired TRA_TRB    
#> 10 clonotype282 CAASALR… TRA;TR… TRAV14N… TRAJ13… 1337… 8;15… paired TRA_TRB    
#> # … with 91 more rows
```

![](man/figures/README-chains_umap-1.png)<!-- -->

<br>

### Clonotype Abundance

To identify the top clonotypes in each sample or cluster, clonotype
abundance can be calculated using the `calc_abundance` function.

``` r
x <- calc_abundance(
  sobj_in       = so_tcr,        # Seurat object
  clonotype_col = "cdr3",        # meta.data column containing clonotype IDs
  cluster_col   = "orig.ident",  # meta.data column containing cell labels
  prefix        = ""             # Prefix to add to new meta.data columns
)
```

<br>

### Repertoire Diversity

The function `calc_diversity` will calculate repertoire diversity on a
per-cluster basis. Using the `cluster_col` argument, any meta.data
column containing cell labels can be used for calculations.
`calc_diversity` uses the R package `abdiv` for performing diversity
calculations and any `abdiv` diversity function can be specified using
the `method` argument.

Possible methods for calculating diversity include:

    #>  [1] "berger_parker_d"  "brillouin_d"      "dominance"        "heip_e"          
    #>  [5] "invsimpson"       "kempton_taylor_q" "margalef"         "mcintosh_d"      
    #>  [9] "mcintosh_e"       "menhinick"        "pielou_e"         "richness"        
    #> [13] "shannon"          "simpson"          "simpson_e"        "strong"

<br>

In this example we are calculating the Shannon index of diversity for
each sample in the orig.ident meta.data column.

``` r
so_tcr <- calc_diversity(
  sobj_in       = so_tcr,         # Seurat object
  clonotype_col = "cdr3",         # meta.data column containing clonotype ids
  cluster_col   = "orig.ident",   # meta.data column containing cell labels
  method        = abdiv::shannon  # abdiv method to use
)
```

<br>

For each ‘calculation’ function provided by `djvdj`, there is a matching
`plot` function that will generate a summary plot. For `plot_diversity`
this is a bar graph.

``` r
plot_diversity(
  sobj_in       = so_tcr,          # Seurat object
  clonotype_col = "cdr3",          # meta.data column containing clonotype ids
  cluster_col   = "orig.ident",    # meta.data column containing cell labels
  method        = abdiv::shannon,  # abdiv method to use
  plot_colors   = ito_cols 
) +
  theme_cowplot() +
  theme(
    legend.position = "none",
    axis.title.x    = element_blank(),
    axis.text.x     = element_text(angle = 45, hjust = 1)
  )
```

![](man/figures/README-div_plots-1.png)<!-- -->

<br>

### Repertoire Overlap

To compare repertoires for different samples or clusters,
`calc_similarity` can calculate a variety of different similarity
metrics. The `cluster_col` should be used to specify the meta.data
column containing cell labels for comparison. Like `calc_diversity`, an
`abdiv` function can be specified with the `method` argument.

Possible methods for calculating repertoire similarity include:

    #>  [1] "binomial_deviance"                  "bray_curtis"                       
    #>  [3] "bray_curtis_balanced"               "bray_curtis_gradient"              
    #>  [5] "canberra"                           "chebyshev"                         
    #>  [7] "chord"                              "clark_coefficient_of_divergence"   
    #>  [9] "correlation_distance"               "cosine_distance"                   
    #> [11] "cy_dissimilarity"                   "euclidean"                         
    #> [13] "geodesic_metric"                    "hamming"                           
    #> [15] "hellinger"                          "horn_morisita"                     
    #> [17] "jaccard"                            "jaccard_nestedness"                
    #> [19] "jaccard_turnover"                   "kulczynski_first"                  
    #> [21] "kulczynski_second"                  "kullback_leibler_divergence"       
    #> [23] "manhattan"                          "mean_character_difference"         
    #> [25] "minkowski"                          "modified_mean_character_difference"
    #> [27] "morisita"                           "rms_distance"                      
    #> [29] "rogers_tanimoto"                    "russel_rao"                        
    #> [31] "ruzicka"                            "ruzicka_balanced"                  
    #> [33] "ruzicka_gradient"                   "sokal_michener"                    
    #> [35] "sokal_sneath"                       "sorenson"                          
    #> [37] "sorenson_nestedness"                "sorenson_turnover"                 
    #> [39] "weighted_kulczynski_second"         "yule_dissimilarity"

<br>

By default `calc_similarity` will add a new meta.data column for each
comparison. In this example we are calculating the jaccard dissimilarity
index for all combinations of cell labels present in the `orig.ident`
column.

``` r
so_tcr <- calc_similarity(
  sobj_in       = so_tcr,          # Seurat object
  clonotype_col = "cdr3",          # meta.data column containing clonotype ids
  cluster_col   = "orig.ident",    # meta.data column containing cell labels
  method        = abdiv::jaccard,  # abdiv method to use
  prefix        = "jcrd_",         # Prefix to add to new meta.data columns 
  return_seurat = TRUE             # Return Seurat object with results added to meta.data
)

# Take a look at the meta.data
so_tcr@meta.data %>%
  as_tibble() %>%
  filter(!is.na(clonotype_id), n_chains > 2) %>%
  select(all_of(vdj_cols), starts_with("jcrd_"))
#> # A tibble: 101 x 13
#>    clonotype_id cdr3  chains v_gene j_gene reads umis  Paired uniq_chains
#>    <chr>        <chr> <chr>  <chr>  <chr>  <chr> <chr> <chr>  <chr>      
#>  1 clonotype10  CATT… TRA;T… TRAV1… TRAJ3… 904;… 2;6;… paired TRA_TRB    
#>  2 clonotype27  CALV… TRA;T… TRAV1… TRAJ2… 1688… 13;1… paired TRA_TRB    
#>  3 clonotype36  CAAY… TRA;T… TRAV1… TRAJ5… 1210… 14;8… paired TRA_TRB    
#>  4 clonotype85  CVVV… TRA;T… TRAV1… TRAJ2… 6974… 7;11… paired TRA_TRB    
#>  5 clonotype136 CAAR… TRA;T… TRAV1… TRAJ3… 2014… 13;3… paired TRA_TRB    
#>  6 clonotype138 CVAH… TRA;T… TRAV1… TRAJ3… 7230… 5;10… paired TRA_TRB    
#>  7 clonotype174 CAAS… TRA;T… TRAV7… TRAJ2… 1815… 11;8… paired TRA_TRB    
#>  8 clonotype205 CALL… TRA;T… TRAV1… TRAJ5… 1294… 19;5… paired TRA_TRB    
#>  9 clonotype215 CAPG… TRA;T… TRAV1… TRAJ1… 1188… 13;4… paired TRA_TRB    
#> 10 clonotype282 CAAS… TRA;T… TRAV1… TRAJ1… 1337… 8;15… paired TRA_TRB    
#> # … with 91 more rows, and 4 more variables: jcrd_KI_DN3 <dbl>,
#> #   jcrd_KI_DN4 <dbl>, jcrd_WT_DN3 <dbl>, jcrd_WT_DN4 <dbl>
```

<br>

Alternatively, `calc_similarity` can output a matrix

``` r
calc_similarity(
  sobj_in       = so_tcr,          # Seurat object
  clonotype_col = "cdr3",          # meta.data column containing clonotype ids
  cluster_col   = "orig.ident",    # meta.data column containing cell labels
  method        = abdiv::jaccard,  # abdiv method to use
  return_seurat = FALSE            # Return Seurat object with results added to meta.data
)
#>           KI_DN3    KI_DN4 WT_DN3
#> KI_DN4 0.9987069        NA     NA
#> WT_DN3 1.0000000 0.9986326     NA
#> WT_DN4 0.9969174 0.9980024      1
```

<br>

A heatmap summarizing the results can be generated using the
`plot_similarity` function. Here we are creating two heatmaps, one to
compare the different samples and one to compare cell clusters.

``` r
# Sample heatmap
ident_heat <- plot_similarity(
  sobj_in       = so_tcr,
  clonotype_col = "cdr3",
  cluster_col   = "orig.ident",
  method        = abdiv::jaccard,
  plot_colors   = c("grey90", "#56B4E9")
) +
  theme(legend.title = element_blank())

# Cluster heatmap
clust_heat <- plot_similarity(
  sobj_in       = so_tcr,
  clonotype_col = "cdr3",
  cluster_col   = "seurat_clusters",
  method        = abdiv::jaccard,
  plot_colors   = c("grey90", "#56B4E9"),
  
  size          = 1,       # Additional options to pass to ggplot
  color         = "white"  # Additional options to pass to ggplot
) +
  theme(legend.title = element_blank())

# Combine heatmaps
plot_grid(ident_heat, clust_heat)
```

![](man/figures/README-sim_plots-1.png)<!-- -->

<br>

### VDJ Gene Usage

The V(D)J data imported from Cell Ranger also includes the specific
genes detected for each cell. The function `calc_usage` can be used to
calculate the fraction of cells that express different V(D)J genes. This
function will produce a table summarizing the results. To only include
results for a certain chain, the `chain` and `chain_col` arguments can
be used to specify the meta.data column containing the chains detected
for each cell. By default the results for all chains will be included.

In this example we are summarizing the usage of different V genes for
the TRB chain

``` r
calc_usage(
  sobj_in     = so_tcr,        # Seurat object
  gene_cols   = "v_gene",      # meta.data column containing genes
  cluster_col = "orig.ident",  # meta.data column containing cell labels
  chain       = "TRB",         # Chain to use for filtering genes
  chain_col   = "chains"       # meta.data column containing chains identified for each cell
)
#> # A tibble: 92 x 5
#>    v_gene orig.ident n_cells  freq   pct
#>    <chr>  <chr>        <dbl> <int> <dbl>
#>  1 None   WT_DN3         597     0 0    
#>  2 None   WT_DN4         909    12 1.32 
#>  3 None   KI_DN3         728     0 0    
#>  4 None   KI_DN4        1616     8 0.495
#>  5 TRBV1  WT_DN3         597    23 3.85 
#>  6 TRBV1  WT_DN4         909    38 4.18 
#>  7 TRBV1  KI_DN3         728    34 4.67 
#>  8 TRBV1  KI_DN4        1616    66 4.08 
#>  9 TRBV10 WT_DN3         597    31 5.19 
#> 10 TRBV10 WT_DN4         909    59 6.49 
#> # … with 82 more rows
```

<br>

The companion function `plot_usage` can be used to create a heatmap
summarizing these results. Using the `yaxis` argument, the percent or
absolute count (frequency) can be used for plotting. The genes plotting
can also be selected using the `plot_genes` argument, or the total
number of genes to plot can be passed to `n_genes`.

In this example we are only plotting the top 10 genes that have the
highest average usage.

``` r
plot_usage(
  sobj_in     = so_tcr,                    # Seurat object
  gene_cols   = "v_gene",                  # meta.data column(s) containing genes
  cluster_col = "orig.ident",              # meta.data column containing cell labels
  chain       = "TRB",                     # Chain to use for filtering genes
  chain_col   = "chains",                  # meta.data column containing chains
  
  plot_colors = c("grey90", ito_cols[5]),  # Colors to use for heatmap
  plot_genes  = NULL,                      # A list of genes to plot
  n_genes     = 10,                        # The number of top genes to plot
  yaxis       = "percent"                  # Units to plot
)
```

![](man/figures/README-usage_plots_1-1.png)<!-- -->

<br>

By passing multiple columns to `gene_cols`, the frequency that different
genes are used together can also be summarized.

``` r
calc_usage(
  sobj_in     = so_tcr,                 # Seurat object
  gene_cols   = c("v_gene", "j_gene"),  # meta.data column(s) containing genes
  cluster_col = "orig.ident",           # meta.data column containing cell labels
  chain       = "TRB",                  # Chain to use for filtering genes
  chain_col   = "chains"                # meta.data column containing chains
)
#> # A tibble: 1,016 x 6
#>    v_gene j_gene  orig.ident n_cells  freq   pct
#>    <chr>  <chr>   <chr>        <dbl> <int> <dbl>
#>  1 None   None    WT_DN3         597     0 0    
#>  2 None   None    WT_DN4         909    12 1.32 
#>  3 None   None    KI_DN3         728     0 0    
#>  4 None   None    KI_DN4        1616     8 0.495
#>  5 TRBV1  TRBJ1-1 WT_DN3         597     3 0.503
#>  6 TRBV1  TRBJ1-1 WT_DN4         909     8 0.880
#>  7 TRBV1  TRBJ1-1 KI_DN3         728     4 0.549
#>  8 TRBV1  TRBJ1-1 KI_DN4        1616    10 0.619
#>  9 TRBV1  TRBJ1-2 WT_DN3         597     1 0.168
#> 10 TRBV1  TRBJ1-2 WT_DN4         909     3 0.330
#> # … with 1,006 more rows
```

<br>

When multiple gene columns are passed to `plot_usage`, a list of plots
will be returned, one for each cell label in the `cluster_col` column.

``` r
gg <- plot_usage(
  sobj_in     = so_tcr,                   # Seurat object
  gene_cols   = c("v_gene", "j_gene"),    # meta.data column(s) containing genes
  cluster_col = "orig.ident",             # meta.data column containing cell labels
  chain       = "TRB",                    # Chain to use for filtering genes
  chain_col   = "chains",                 # meta.data column containing chains identified
  plot_colors = c("grey90", ito_cols[8])  # Colors to use for heatmap
) %>%
  imap(~ .x + ggtitle(.y))

# Combine heatmaps
plot_grid(plotlist = gg)
```

![](man/figures/README-usage_plots_2-1.png)<!-- -->
