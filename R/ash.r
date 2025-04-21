
#' Append lfsr results from ashr to a data frame
#' 
#' This function appends lfsr results from the ashr package to a data frame containing GWAS summary statistics.
#' 
#' @param x A data frame containing GWAS summary statistics with columns "BETA" and "SE".
#' @param ... Additional arguments passed to the method.
#' 
#' @return A data frame with additional columns for lfsr, qvalue, pm, and psd.
#' @export
append_ashr_results <- function(x, ...) {
    UseMethod("append_lfsr")
}

#' @export
append_ashr_results.data.frame <- function(x, ...) {

    if (!requireNamespace("ashr", quietly = TRUE)) {
        cli::cli_abort("The 'ashr' package is required for this function. Please install it.")
    }
    
    # Check if the required columns are present
    if (!all(c("BETA", "SE") %in% names(x))) {
        cli::cli_abort("The data frame must contain 'BETA' and 'SE' columns.")
    }
    
    lfsr_result <- ashr::ash(x$BETA, x$SE)

    # Append lfsr to the data frame
    x$lfsr <- ashr::get_lfsr(lfsr_result)
    x$qvalue <- ashr::get_qvalue(lfsr_result)
    x$pm <- ashr::get_pm(lfsr_result)
    x$psd <- ashr::get_psd(lfsr_result)
    
    return(x)
}

#' @export
append_ashr_results.tbl_df <- function(x, ...) {
    append_ashr_results.data.frame(x, ...)# Check if the required columns are present
}

#' @export
append_ashr_results.GWASFormatter <- function(x, ...) {
    cli::cli_abort("The 'append_ashr_results' method is not implemented for GWASFormatter objects.")
}