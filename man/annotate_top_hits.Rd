% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/annotate.R
\name{annotate_top_hits}
\alias{annotate_top_hits}
\title{Pull top hits from a gwas object}
\usage{
annotate_top_hits(gwas, threshold = 5e-08)
}
\arguments{
\item{gwas}{A gwas object containing the data to plot.}

\item{threshold}{The p-value threshold to filter the top hits. Default is 5e-8.}
}
\value{
A data frame with the top hits.
}
\description{
Pull top hits from a gwas object
}
