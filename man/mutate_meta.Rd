% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{mutate_meta}
\alias{mutate_meta}
\title{Manipulate Seurat object meta.data}
\usage{
mutate_meta(sobj_in, .fun, ...)
}
\arguments{
\item{sobj_in}{Seurat object}

\item{.fun}{Function or formula to use for modifying the meta.data. If a
formula is provided, use .x to refer to the meta.data table.}

\item{...}{Arguments to pass to the provided function}
}
\value{
Seurat object
}
\description{
Manipulate Seurat object meta.data
}