% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/annotate.R
\name{annotate_with_centromere}
\alias{annotate_with_centromere}
\alias{annotate_with_centromere.data.frame}
\alias{annotate_with_centromere.tbl_df}
\alias{annotate_with_centromere.GWASFormatter}
\title{Annotate data with centromere information}
\usage{
annotate_with_centromere(x, ...)

\method{annotate_with_centromere}{data.frame}(x, chrom_col = "CHROM", pos_col = "POS", ...)

\method{annotate_with_centromere}{tbl_df}(x, ...)

\method{annotate_with_centromere}{GWASFormatter}(x, ...)
}
\arguments{
\item{x}{A data frame or tibble containing variant data.}

\item{...}{Additional arguments passed to methods.}

\item{chrom_col}{Name of the column containing chromosome information. Default is "CHROM".}

\item{pos_col}{Name of the column containing position information. Default is "POS".}
}
\value{
A data frame with the centromere information.
}
\description{
Annotate data with centromere information
}
\section{Methods (by class)}{
\itemize{
\item \code{annotate_with_centromere(data.frame)}: Annotate a data frame or tibble with centromere information

\item \code{annotate_with_centromere(tbl_df)}: Alias for the data.frame method

\item \code{annotate_with_centromere(GWASFormatter)}: Centromere annotation method for GWASFormatter objects

}}
