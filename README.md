
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

## Vignette

Splenocytes from MD4 transgenic mice which have monoclonal B cells that
all bind hen egg lysozyme (HEL) antigen were mixed with splenocytes from
C57BL/6 mice at a 1:20 ratio. The cells were stained with the HEL
AVID-tag and sequencing libraries were prepared to capture gene
expression, B cell receptor sequences, and AVID-tag signals using the
10x Genomics 5’ immune profiling kit.

<img src="man/figures/README-rna_umap-1.png" width="100%" />

<br>

### Import VDJ data

`import_vdj` takes the output files from `cellranger vdj` and adds
clonotype information to the meta.data for an existing Seurat object.
For cells that do not have any VDJ sequencing data, NAs will be included
in the meta.data. The `filter_contigs` argument will only include chains
that have at least one contig that is full length and productive.

``` r
so_vdj <- import_vdj(
  sobj_in        = so,              # Seurat object                         
  vdj_dir        = params$vdj_dir,  # Directory containing cellranger output files
  prefix         = "",              # Prefix to add to new meta.data columns
  cell_prefix    = "",              # Prefix to add to cell barcodes
  filter_contigs = TRUE             # Only include chains with at least one productive contig
)

vdj_cols <- c(
  "chain",      "cdr3",
  "clone_freq", "clone_prop",
  "n_chains"
)

so_vdj@meta.data %>%
  as_tibble() %>%
  select(orig.ident, nCount_RNA, nFeature_RNA, all_of(vdj_cols))
#> # A tibble: 7,137 x 8
#>    orig.ident nCount_RNA nFeature_RNA chain cdr3  clone_freq clone_prop n_chains
#>    <fct>           <dbl>        <int> <chr> <chr>      <int>      <dbl>    <int>
#>  1 AVID-seq          884          551 IGH;… CVKG…          1   0.000262        2
#>  2 AVID-seq         3061          970 <NA>  <NA>          NA  NA              NA
#>  3 AVID-seq         1297          677 IGH;… CARG…          1   0.000262        2
#>  4 AVID-seq         1570          848 <NA>  <NA>          NA  NA              NA
#>  5 AVID-seq         2277          818 <NA>  <NA>          NA  NA              NA
#>  6 AVID-seq         1320          566 <NA>  <NA>          NA  NA              NA
#>  7 AVID-seq          570          348 <NA>  <NA>          NA  NA              NA
#>  8 AVID-seq          909          489 IGH;… CTVS…          1   0.000262        2
#>  9 AVID-seq         1072          588 IGH;… CARS…          1   0.000262        2
#> 10 AVID-seq         1143          594 IGH;… CARS…          1   0.000262        2
#> # … with 7,127 more rows
```

<br>

### Filtering

`filter_vdj` allows you to filter a Seurat object using the added
clonotype information or any other columns present in the meta.data. For
cells with multiple chains, the information for each chain is stored as
a single row, separated by a “;”. When filtering, columns with VDJ data
will be expanded based on the delimiter “;”. The columns that are
expanded for filtering can be specified with the `split_cols` argument.
By default filtering is only performed on cells that include VDJ data.

Filter to only include cells with paired light and heavy chains.

``` r
so_filt <- filter_vdj(
  sobj_in = so_vdj,                                              # Seurat object
  filt    = "IGH" %in% chain && any(c("IGK", "IGL") %in% chain)  # Expression for filtering
)

so_filt@meta.data %>%
  as_tibble() %>%
  filter(!is.na(clonotype_id)) %>%
  select(all_of(vdj_cols))
#> # A tibble: 3,353 x 5
#>    chain       cdr3                               clone_freq clone_prop n_chains
#>    <chr>       <chr>                                   <int>      <dbl>    <int>
#>  1 IGH;IGK     CVKGYDYDWYFDVW;CLQYDNLWTF                   1   0.000262        2
#>  2 IGH;IGK     CARGRLGYAMDYW;CQHFWSTPWTF                   1   0.000262        2
#>  3 IGH;IGK     CTVSYTKDWYFDVW;CAQNLELPLTF                  1   0.000262        2
#>  4 IGH;IGK     CARSYDYDPLYYAMDYW;CLQSDNLPLTF               1   0.000262        2
#>  5 IGH;IGK     CARSRLAYW;CLQYASSPFTF                       1   0.000262        2
#>  6 IGH;IGK;IGL CAKRGYSNSLDYW;CQHFWSTPYTF;CALWYSN…          1   0.000262        3
#>  7 IGH;IGK     CANPITTAEGWYFDVW;CLQHGESPYTF                1   0.000262        2
#>  8 IGH;IGK     CARSYGYAMDYW;CWQGTHFPYTF                    1   0.000262        2
#>  9 IGH;IGK     CARWVYGSAWFAYW;CMQHLEYPFTF                  1   0.000262        2
#> 10 IGH;IGK     CARSHGYDFYAMDYW;CQHFWGTPRTF                 1   0.000262        2
#> # … with 3,343 more rows
```

<br>

Instead of filtering, `filter_vdj` can also add new cell labels to the
object meta.data using the `new_col` argument. Here a new column is
added indicating whether each cell has a paired heavy and light chain.
This is useful for plotting.

``` r
so_vdj <- filter_vdj(
  sobj_in = so_vdj,                                               # Seurat object
  filt    = "IGH" %in% chain && any(c("IGK", "IGL") %in% chain),  # Condition to use for filtering
  new_col = "Paired",                                             # Name of new column
  true    = "paired",                                             # Value when condition is TRUE
  false   = "unpaired"                                            # Value when condition is FALSE
)

vdj_cols <- c(vdj_cols, "Paired")

so_vdj@meta.data %>%
  as_tibble() %>%
  filter(!is.na(clonotype_id)) %>%
  select(all_of(vdj_cols))
#> # A tibble: 3,820 x 6
#>    chain      cdr3                        clone_freq clone_prop n_chains Paired 
#>    <chr>      <chr>                            <int>      <dbl>    <int> <chr>  
#>  1 IGH;IGK    CVKGYDYDWYFDVW;CLQYDNLWTF            1   0.000262        2 paired 
#>  2 IGH;IGK    CARGRLGYAMDYW;CQHFWSTPWTF            1   0.000262        2 paired 
#>  3 IGH;IGK    CTVSYTKDWYFDVW;CAQNLELPLTF           1   0.000262        2 paired 
#>  4 IGH;IGK    CARSYDYDPLYYAMDYW;CLQSDNLP…          1   0.000262        2 paired 
#>  5 IGH;IGK    CARSRLAYW;CLQYASSPFTF                1   0.000262        2 paired 
#>  6 IGH;IGK;I… CAKRGYSNSLDYW;CQHFWSTPYTF;…          1   0.000262        3 paired 
#>  7 IGK        CQQWSSNPLTF                          3   0.000785        1 unpair…
#>  8 IGH;IGK    CANPITTAEGWYFDVW;CLQHGESPY…          1   0.000262        2 paired 
#>  9 IGH;IGK    CARSYGYAMDYW;CWQGTHFPYTF             1   0.000262        2 paired 
#> 10 IGH;IGK    CARWVYGSAWFAYW;CMQHLEYPFTF           1   0.000262        2 paired 
#> # … with 3,810 more rows
```

<img src="man/figures/README-pair_umap-1.png" width="100%" />

<br>

More complicated statements referring to meta.data columns can be used
for the `filt`, `true`, and `false` arguments. For more detailed
analysis of the chains detected for each cell, a new cell label can be
created for only the unique chains.

``` r
so_vdj <- filter_vdj(
  sobj_in = so_vdj,                                # Seurat object
  filt    = length(unique(chain)) < 4,             # Condition to use for filtering
  new_col = "uniq_chains",                         # Name of new column
  true    = str_c(unique(chain), collapse = "_"),  # Value when condition is TRUE
  false   = "other"                                # Value when condition is FALSE
)

vdj_cols <- c(vdj_cols, "uniq_chains")

so_vdj@meta.data %>%
  as_tibble() %>%
  filter(!is.na(clonotype_id), n_chains > 2) %>%
  select(all_of(vdj_cols))
#> # A tibble: 526 x 7
#>    chain    cdr3               clone_freq clone_prop n_chains Paired uniq_chains
#>    <chr>    <chr>                   <int>      <dbl>    <int> <chr>  <chr>      
#>  1 IGH;IGK… CAKRGYSNSLDYW;CQH…          1   0.000262        3 paired IGH_IGK_IGL
#>  2 IGH;IGH… CARGDYW;CTTWLRLRS…          1   0.000262        4 paired IGH_IGK    
#>  3 IGH;IGK… CAKPRYYYGSSFYAMDY…          1   0.000262        3 paired IGH_IGK    
#>  4 IGH;IGK… CARGPYYTNGGAMDYW;…          1   0.000262        3 paired IGH_IGK    
#>  5 IGH;IGH… CARSYPYFDYW;CARSS…          1   0.000262        4 paired IGH_IGK    
#>  6 IGH;IGK… CALDSSGFAYW;CQQYW…          1   0.000262        3 paired IGH_IGK    
#>  7 IGH;IGH… CARHDGLPGAMDYW;CA…          1   0.000262        4 paired IGH_IGK    
#>  8 IGH;IGH… CAEGSSNWYFDVW;CAR…          1   0.000262        4 paired IGH_IGK    
#>  9 IGH;IGK… CTSPPYEGYYAMDYW;C…          1   0.000262        3 paired IGH_IGK    
#> 10 IGH;IGH… CTRLLTGYYFDYW;CAR…          1   0.000262        4 paired IGH_IGK    
#> # … with 516 more rows
```

<img src="man/figures/README-chains_umap-1.png" width="100%" />

<br>

MD4 B cells are expected to have IGK chains with the CDR3 amino acid
sequence CQQSNSWPYTF. Using `filter_vdj` we can visualize which cells
have this sequence.

``` r
# Add new cell labels to meta.data
so_vdj <- filter_vdj(
  sobj_in = so_vdj,                                   # Seurat object
  filt    = "CQQSNSWPYTF" %in% cdr3[chain == "IGK"],  # Condition to use for filtering
  new_col = "IGK_seq",                                # Name of new column
  true    = "CQQSNSWPYTF",                            # Value when condition is TRUE
  false   = "other"                                   # Value when condition is FALSE
)

vdj_cols <- c(vdj_cols, "uniq_chains")

so_vdj@meta.data %>%
  as_tibble() %>%
  filter(!is.na(clonotype_id), n_chains > 2) %>%
  select(all_of(vdj_cols))
#> # A tibble: 526 x 7
#>    chain    cdr3               clone_freq clone_prop n_chains Paired uniq_chains
#>    <chr>    <chr>                   <int>      <dbl>    <int> <chr>  <chr>      
#>  1 IGH;IGK… CAKRGYSNSLDYW;CQH…          1   0.000262        3 paired IGH_IGK_IGL
#>  2 IGH;IGH… CARGDYW;CTTWLRLRS…          1   0.000262        4 paired IGH_IGK    
#>  3 IGH;IGK… CAKPRYYYGSSFYAMDY…          1   0.000262        3 paired IGH_IGK    
#>  4 IGH;IGK… CARGPYYTNGGAMDYW;…          1   0.000262        3 paired IGH_IGK    
#>  5 IGH;IGH… CARSYPYFDYW;CARSS…          1   0.000262        4 paired IGH_IGK    
#>  6 IGH;IGK… CALDSSGFAYW;CQQYW…          1   0.000262        3 paired IGH_IGK    
#>  7 IGH;IGH… CARHDGLPGAMDYW;CA…          1   0.000262        4 paired IGH_IGK    
#>  8 IGH;IGH… CAEGSSNWYFDVW;CAR…          1   0.000262        4 paired IGH_IGK    
#>  9 IGH;IGK… CTSPPYEGYYAMDYW;C…          1   0.000262        3 paired IGH_IGK    
#> 10 IGH;IGH… CTRLLTGYYFDYW;CAR…          1   0.000262        4 paired IGH_IGK    
#> # … with 516 more rows
```

<img src="man/figures/README-seq_umap-1.png" width="100%" />

<br>

### Repertoire stats

The functions `calc_diversity` and `calc_jaccard` will calculate
repertoire diversity and repertoire overlap on a per-cluster basis.
These functions can be given any meta.data column containing cell labels
to use for calculations.

Calculate repertoire diversity with `calc_diversity`. The inverse
Simpson index is used to measure diversity for each cluster.

``` r
so_vdj <- calc_diversity(
  sobj_in       = so_vdj,          # Seurat object
  clonotype_col = "clonotype_id",  # meta.data column containing clonotype ids
  cluster_col   = "type_mouse",    # meta.data column containing cell labels
  prefix        = ""               # Prefix to add to new meta.data columns
)
```

<img src="man/figures/README-div_umap-1.png" width="100%" />

<br>

Calculate repertoire overlap with `calc_jaccard` for the cell groups
present in the `cluster_col`. Using the `return_seurat` argument,
`calc_jaccard` can also output a matrix for plotting.

``` r
so_vdj <- calc_jaccard(
  sobj_in       = so_vdj,             # Seurat object
  clonotype_col = "clonotype_id",     # meta.data column containing clonotype ids
  cluster_col   = "seurat_clusters",  # meta.data column containing cell labels
  prefix        = "x",                # Prefix to add to new meta.data columns 
  return_seurat = TRUE                # Return Seurat object with results added to meta.data
)
```

<img src="man/figures/README-jaccard_umap-1.png" width="100%" />
