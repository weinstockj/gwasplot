\name{meta_analyze_fe}
\alias{meta_analyze_fe}
\alias{meta_analyze_fe.data.frame}
\alias{meta_analyze_fe.tbl_df}
\alias{meta_analyze_fe.GWASFormatter}
\title{Fixed-effects meta-analysis of two GWAS datasets}
\usage{
meta_analyze_fe(x, y, ...)

\method{meta_analyze_fe}{data.frame}(
  x,
  y,
  keep = c("overlap"),
  allow_swaps = TRUE,
  ...
)

\method{meta_analyze_fe}{tbl_df}(
  x,
  y,
  keep = c("overlap"),
  allow_swaps = TRUE,
  ...
)

\method{meta_analyze_fe}{GWASFormatter}(
  x,
  y,
  keep = c("overlap"),
  allow_swaps = TRUE,
  ...
)
}
\arguments{
\item{x}{A \code{GWASFormatter} object, data.frame, or tibble containing standardized GWAS results.}

\item{y}{A second \code{GWASFormatter} object, data.frame, or tibble containing standardized GWAS results.}

\item{keep}{Output row-set policy. Currently only \code{"overlap"} is supported.}

\item{allow_swaps}{Logical. If \code{TRUE}, simple \code{REF}/\code{ALT} swaps are harmonized by flipping the second study's \code{BETA}.}

\item{...}{Additional arguments passed to methods.}
}
\value{
For \code{GWASFormatter} inputs, a modified \code{GWASFormatter} object whose data contains the meta-analysis result. For eager inputs, a tibble containing the harmonized overlap and fixed-effects summary statistics.
}
\description{
Perform inverse-variance fixed-effects meta-analysis on two GWAS datasets that have been standardized to the package's expected schema.
}
\details{
The current implementation keeps only the harmonized overlap between the two studies. Exact allele matches are used directly. Simple \code{REF}/\code{ALT} swaps are supported by flipping the second study's \code{BETA}. The returned result preserves the package's standard GWAS columns, including \code{CHROM}, \code{POS}, \code{REF}, \code{ALT}, \code{ID}, \code{BETA}, \code{SE}, and \code{PVALUE}, so it can be used by downstream plotting and annotation functions.

When both inputs are \code{GWASFormatter} objects, the harmonization and fixed-effects computation are performed inside DuckDB against the objects' backing tables. Each \code{GWASFormatter} now owns its own DuckDB table rather than sharing a single global \code{summary_stats} table.

The first implementation supports two studies only, inverse-variance fixed-effects meta-analysis only, and the overlap-only row set. Variants that cannot be harmonized by exact allele match or simple \code{REF}/\code{ALT} swap are excluded from the result.
}