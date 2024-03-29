% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calc-similarity.R
\name{plot_mds}
\alias{plot_mds}
\title{Create MDS plot}
\usage{
plot_mds(
  input,
  data_col,
  cluster_col,
  method = "jaccard",
  chain = NULL,
  chain_col = global$chain_col,
  plot_colors = NULL,
  plot_lvls = names(plot_colors),
  label_points = TRUE,
  n_label = "none",
  label_params = list(),
  sep = global$sep,
  ...
)
}
\arguments{
\item{input}{Single cell object or data.frame containing V(D)J data. If a
data.frame is provided, the cell barcodes should be stored as row names.}

\item{data_col}{meta.data column containing values to use for
calculating pairwise similarity between clusters, e.g. 'clonotype_id'}

\item{cluster_col}{meta.data column containing cluster IDs to use for
calculating overlap}

\item{method}{Method to use for comparing clusters and calculating MDS
coordinates, available methods include:
\itemize{
\item 'jaccard', Jaccard dissimilarity index implemented with \code{\link[abdiv:jaccard]{abdiv::jaccard()}}
\item 'horn_morisita', Horn-Morisita index implemented with
\code{\link[abdiv:morisita]{abdiv::horn_morisita()}}
}}

\item{chain}{Chain to use for comparing clusters. To perform calculations
using a single chain, the column passed to \code{data_col} must contain
per-chain data such as CDR3 sequences. Set to \code{NULL} to include all chains.}

\item{chain_col}{meta.data column containing chains for each cell}

\item{plot_colors}{Character vector containing colors for plotting}

\item{plot_lvls}{Levels to use for ordering clusters}

\item{label_points}{Label points on plot}

\item{n_label}{Location on plot where n label should be added, this can be
any combination of the following:
\itemize{
\item 'corner', display the total number of points plotted in the top right
corner, the position of the label can be modified by passing \code{x} and \code{y}
specifications with the \code{label_params} argument
\item 'none', do not display the number of points plotted
}}

\item{label_params}{Named list providing additional parameters to modify
n label aesthetics, e.g. list(size = 4, color = "red")}

\item{sep}{Separator used for storing per-chain V(D)J data for each cell}

\item{...}{Additional arguments to pass to \code{\link[ggplot2:geom_point]{ggplot2::geom_point()}}}
}
\value{
ggplot object
}
\description{
Calculate MDS coordinates based on a beta diversity metric and plot results.
}
\examples{
# Calculate MDS coordinates
plot_mds(
  vdj_sce,
  data_col = "clonotype_id",
  cluster_col = "isotype"
)

# Calculate MDS coordinates based on IGK CDR3 sequences
plot_mds(
  vdj_sce,
  data_col    = "cdr3",
  cluster_col = "isotype",
  chain       = "IGK"
)

# Calculate repertoire similarity using the Horn-Morisita index
plot_mds(
  vdj_sce,
  data_col    = "clonotype_id",
  cluster_col = "isotype",
  method      = "horn_morisita"
)

}
\seealso{
\code{\link[=calc_mds]{calc_mds()}}, \code{\link[=calc_similarity]{calc_similarity()}}, \code{\link[=plot_similarity]{plot_similarity()}},
\code{\link[MASS:isoMDS]{MASS::isoMDS()}}
}
