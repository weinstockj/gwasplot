# modified from here: https://stackoverflow.com/questions/27788968/how-would-one-check-the-system-memory-available-using-r-on-a-windows-machine
getFreeMemoryGB <- function() {
  osName <- Sys.info()[["sysname"]]
  if (osName == "Windows") {
    x <- system2(
      "wmic",
      args = "OS get FreePhysicalMemory /Value",
      stdout = TRUE
    )
    x <- x[grepl("FreePhysicalMemory", x)]
    x <- gsub("FreePhysicalMemory=", "", x, fixed = TRUE)
    x <- gsub("\r", "", x, fixed = TRUE)
    return(as.integer(x) / (1024^2))
  } else if (osName == 'Linux') {
    x <- system2('free', args = '-k', stdout = TRUE)
    x <- strsplit(x[2], " +")[[1]][4]
    return(as.integer(x) / (1024^2))
  } else {
    stop("Only supported on Windows and Linux")
  }
}

#' Compute the genomic inflation factor (lambda GC)
#'
#' Calculates lambda GC as the ratio of the median observed chi-squared statistic
#' (derived from p-values) to the expected median under the null (0.4549).
#' A value near 1.0 indicates no inflation; values above ~1.05 suggest
#' population stratification or other confounding.
#'
#' @param x A numeric vector of p-values, a data.frame/tibble with a `PVALUE`
#'   column, or a GWASFormatter object.
#' @param ... Additional arguments (unused).
#' @return A single numeric value: the genomic inflation factor.
#' @export
lambda_gc <- function(x, ...) UseMethod("lambda_gc")

#' @describeIn lambda_gc Method for numeric vectors of p-values
#' @export
lambda_gc.numeric <- function(x, ...) {
  chi2 <- qchisq(x, df = 1, lower.tail = FALSE)
  val <- median(chi2, na.rm = TRUE) / qchisq(0.5, df = 1)
  cli::cli_inform("lambda_GC = {round(val, 4)}")
  val
}

#' @describeIn lambda_gc Method for data.frame/tibble with a PVALUE column
#' @export
lambda_gc.data.frame <- function(x, ...) {
  lambda_gc.numeric(x$PVALUE, ...)
}

#' @describeIn lambda_gc Method for GWASFormatter objects
#' @export
lambda_gc.GWASFormatter <- function(x, ...) {
  pvals <- x$data %>% dplyr::pull(PVALUE)
  lambda_gc.numeric(pvals, ...)
}

#' Take a random sample of rows from the summary stats table
#'
#' @param g A GWASFormatter object.
#' @param rows The number of rows to sample. Default is 100.
compute_sample = function(g, rows = 100L) {
  df = DBI::dbGetQuery(
    g$con,
    glue("SELECT * FROM {g$table_name} USING SAMPLE {rows}")
  )

  df
}
