% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plots.R
\name{plot_reads}
\alias{plot_reads}
\title{Plot read support for chains}
\usage{
plot_reads(
  sobj_in,
  data_cols = c("reads", "umis"),
  chain_col = "chains",
  cluster_col = NULL,
  type = "violin",
  plot_colors = NULL,
  plot_lvls = NULL,
  facet_rows = 1,
  facet_scales = "free_x",
  ...,
  sep = ";"
)
}
\arguments{
\item{sobj_in}{Seurat object}

\item{data_cols}{meta.data columns containing UMI and/or read counts}

\item{chain_col}{meta.data column containing chains. If chain_col is
provided, reads and UMIs will be plotted separately for each chain}

\item{cluster_col}{meta.data column containing cluster IDs. If cluster_col
is provided reads and UMIs will be plotted separately for each cluster.}

\item{type}{Type of plot to create, can be 'violin', 'histogram' or 'density'}

\item{plot_colors}{Character vector containing colors for plotting}

\item{plot_lvls}{Character vector containing levels for ordering}

\item{facet_rows}{The number of facet rows. Use this argument if both
chain_col and cluster_col are provided}

\item{facet_scales}{This argument passes a scales specification to
facet_wrap. Can be "fixed", "free", "free_x", or "free_y"}

\item{...}{Additional arguments to pass to ggplot}

\item{sep}{Separator to use for expanding data_cols}
}
\value{
ggplot object
}
\description{
Plot read support for chains
}
