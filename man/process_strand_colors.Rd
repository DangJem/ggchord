% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helpers_params.R
\name{process_strand_colors}
\alias{process_strand_colors}
\title{Process gene color parameters in strand mode}
\usage{
process_strand_colors(gene_colors)
}
\arguments{
\item{gene_colors}{Color vector (can be NULL, single value, vector of length 2, named vector with "+"/"-")}
}
\value{
Named vector (names are "+"/"-"), standardized color values (default "+" is red, "-" is blue)
}
\description{
Standardizes gene color parameters in strand mode (color by strand direction) into a named vector with "+"/"-"
}
\keyword{internal}
