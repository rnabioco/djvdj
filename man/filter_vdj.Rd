% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/filter-vdj.R
\name{filter_vdj}
\alias{filter_vdj}
\title{Filter V(D)J data in object}
\usage{
filter_vdj(
  input,
  filt,
  data_cols = NULL,
  clonotype_col = "clonotype_id",
  sep = global$sep,
  per_chain = TRUE
)
}
\arguments{
\item{input}{Single cell object or data.frame containing V(D)J data. If a
data.frame is provided, cell barcodes should be stored as row names.}

\item{filt}{Condition to use for filtering V(D)J data. To allow for
filtering of per-chain V(D)J data, the data for each cell is converted into
a vector, e.g. 'IGH;IGK' is equivalent to c('IGH', 'IGK'). This allows R
vector operations to be performed on the per-chain values. The filtering
condition must return TRUE/FALSE for each chain or a single TRUE/FALSE for
each cell. Data can be filtered based on cell barcodes by referring to the
'.cell_id' column.}

\item{data_cols}{meta.data columns containing V(D)J data to use for
filtering. If NULL, V(D)J data are automatically selected by identifying
columns that have NAs in the same rows as clonotype_col.}

\item{clonotype_col}{meta.data column containing clonotype IDs. This column
is used to determine which columns contain V(D)J data.}

\item{sep}{Separator used for storing per cell V(D)J data.}

\item{per_chain}{If \code{TRUE} per-chain data will be parsed so values for
individual chains can be filtered, if \code{FALSE} values will be filtered as is.}
}
\value{
Object with filtered meta.data
}
\description{
Remove V(D)J data for chains/cells that do not satisfy the provided
condition
}
\examples{
# Only include V(D)J data for productive chains
res <- filter_vdj(vdj_sce, productive)

# Only include V(D)J data for cells with paired chains
res <- filter_vdj(vdj_sce, paired)

# Only include V(D)J data for cells with at least one heavy and one light
# chain
res <- filter_vdj(
  vdj_sce,
  "IGH" \%in\% chains && any(c("IGK", "IGL") \%in\% chains)
)

# Only include V(D)J data for cells that have an IGH, IGK, and IGL chain
res <- filter_vdj(
  vdj_sce,
  all(c("IGH", "IGK", "IGL") \%in\% chains)
)

# Only include V(D)J data for heavy chains
res <- filter_vdj(vdj_sce, chains == "IGH")

# Remove chains that only have 1 UMI for support
res <- filter_vdj(vdj_sce, umis > 1)

# Filter based on cell barcode
res <- filter_vdj(vdj_sce, .cell_id == "1_ACGGAGACATGCTGGC-1")

}
