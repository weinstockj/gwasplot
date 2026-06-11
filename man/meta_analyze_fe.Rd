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
  keep = c("overlap", "union"),
  allow_swaps = TRUE,
  ...
)

\method{meta_analyze_fe}{tbl_df}(
  x,
  y,
  keep = c("overlap", "union"),
  allow_swaps = TRUE,
  ...
)

\method{meta_analyze_fe}{GWASFormatter}(
  x,
  y,
  keep = c("overlap", "union"),
  allow_swaps = TRUE,
  ...
)
}
\arguments{
\item{x}{A \code{GWASFormatter} object, data.frame, or tibble containing standardized GWAS results.}

\item{y}{A second \code{GWASFormatter} object, data.frame, or tibble containing standardized GWAS results.}

\item{keep}{Output row-set policy. \code{"overlap"} (default) keeps only variants
present in both studies (intersection). \code{"union"} keeps variants present in
either study; a variant found in only one study passes through with that study's
estimate and \code{N_studies = 1} (and \code{matched_by} of \code{"x_only"} or
\code{"y_only"}).}

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
By default (\code{keep = "overlap"}) only the harmonized intersection of the two
studies is returned. With \code{keep = "union"}, variants present in either study
are retained: overlapping variants are meta-analyzed as usual, while a variant
seen in only one study passes through with that study's estimate and
\code{N_studies = 1}. Exact allele matches are used directly. Simple
\code{REF}/\code{ALT} swaps are supported by flipping the second study's
\code{BETA}. The returned result preserves the package's standard GWAS columns,
including \code{CHROM}, \code{POS}, \code{REF}, \code{ALT}, \code{ID},
\code{BETA}, \code{SE}, and \code{PVALUE}, so it can be used by downstream
plotting and annotation functions.

When both inputs are \code{GWASFormatter} objects, the harmonization and fixed-effects computation are performed inside DuckDB against the objects' backing tables. Each \code{GWASFormatter} now owns its own DuckDB table rather than sharing a single global \code{summary_stats} table.

This implementation supports two studies and inverse-variance fixed-effects meta-analysis. For \code{keep = "overlap"}, variants that cannot be harmonized by exact allele match or simple \code{REF}/\code{ALT} swap are excluded from the result.
}