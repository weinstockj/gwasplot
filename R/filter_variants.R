#' @title Filter variants
#' @description Filter variants in a GWASFormatter object based on a whitelist or blacklist of variants.
#' @param x A GWASFormatter object.
#' @param subset A file path to a whitelist of variants in Parquet format. Only variants in this file will be kept.
#' @export
filter_variants = function(x, subset = NULL, exclude = NULL, ...) {
  UseMethod("filter_variants")
}

#' @export
filter_variants.GWASFormatter = function(x, subset = NULL, exclude = NULL, ...) {

    if(!is.null(subset) && !is.null(exclude)) {
        cli::cli_abort("You cannot specify both subset and exclude arguments.")
    }

    start_time <- Sys.time()
    cli::cli_process_start("Filtering variants")

    if(!is.null(subset)) {

        whitelist = dplyr::tbl(x$con, glue::glue("read_parquet('{subset}')"))

        x$data = dplyr::semi_join(x$data, whitelist, by = c("chrom", "pos", "ref", "alt"))
            dplyr::compute(name = "filtered_variants", temporary = FALSE, overwrite = TRUE)
    }

    if(!is.null(exclude)) {
        
        blacklist = dplyr::tbl(x$con, glue::glue("read_parquet('{exclude}')"))

        x$data = dplyr::anti_join(x$data, blacklist, by = c("chrom", "pos", "ref", "alt"))
            dplyr::compute(name = "filtered_variants", temporary = FALSE, overwrite = TRUE)
    }

    end_time <- Sys.time()
    elapsed <- round(difftime(end_time, start_time, units = "secs"), 2)

    cli::cli_alert_success("Filtered variants")

    return(x)
}
