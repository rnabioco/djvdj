% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plots.R
\name{.add_lm}
\alias{.add_lm}
\title{Add regression line to ggplot}
\usage{
.add_lm(gg_in, lab_pos = NULL, lab_size = 3.5, ...)
}
\arguments{
\item{gg_in}{ggplot object}

\item{lab_pos}{Position of correlation coefficient label. Set to NULL to
omit label}

\item{lab_size}{Size of label}

\item{...}{Additional arguments to pass to geom_smooth}
}
\value{
ggplot object with added regression line
}
\description{
Add regression line to ggplot
}
