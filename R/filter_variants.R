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

        x$data = dplyr::semi_join(x$data, whitelist, by = c("chrom", "pos", "ref", "alt")) %>%
            dplyr::compute(name = "filtered_variants", temporary = FALSE, overwrite = TRUE)
    }

    if(!is.null(exclude)) {

        blacklist = dplyr::tbl(x$con, glue::glue("read_parquet('{exclude}')"))

        x$data = dplyr::anti_join(x$data, blacklist, by = c("chrom", "pos", "ref", "alt")) %>%
            dplyr::compute(name = "filtered_variants", temporary = FALSE, overwrite = TRUE)
    }

    end_time <- Sys.time()
    elapsed <- round(difftime(end_time, start_time, units = "secs"), 2)

    cli::cli_alert_success("Filtered variants")

    return(x)
}

#' @title Exclude difficult regions
#' @description Exclude difficult regions from a GWASFormatter object.
#' @param x A GWASFormatter object.
#' @param ... Additional arguments passed to the method.
#' @export
exclude_difficult_regions = function(x, ...) {
    UseMethod("exclude_difficult_regions")
}

#' @title Exclude difficult regions
#' @param beds_exclude A character vector of bed files to exclude. Options are "hg19diff", "UCSC_unusual", "GRC_exclusions", and "GIAB_difficult_regions".
#' @param active_table The name of the active table to exclude difficult regions from. Options are "filtered_variants" and "summary_stats".
#' 
#' @export
exclude_difficult_regions.GWASFormatter = function(x, beds_exclude = c("hg19diff", "UCSC_unusual", "GRC_exclusions"), active_table = "filtered_variants",  ...) {
  start_time <- Sys.time()
  cli::cli_alert_info("Now excluding 'difficult' regions...")

  active_table = match.arg(active_table, c("filtered_variants", "summary_stats"))
  stopifnot(all(beds_exclude %in% c("hg19diff", "UCSC_unusual", "GRC_exclusions", "GIAB_difficult_regions")))
  
  con = db_connect()

  active_tbl = dplyr::tbl(con, active_table)

  for (bed in beds_exclude) {
    data(list = bed)

    bed_tbl = dplyr::copy_to(con, get(bed), name = bed, overwrite = TRUE, temporary = FALSE) %>%
      dplyr::select(chrom, start, end) %>%
      dplyr::mutate(start = start + 1) # convert from 0-based to 1-based

    active_tbl = dplyr::anti_join(active_tbl, bed_tbl, by = dplyr::join_by(chrom, between(x$POS, y$start, y$end))) 
  }

  active_tbl = active_tbl %>%
      dplyr::compute(name = "filtered_variants", temporary = FALSE, overwrite = TRUE) # now executive lazy query

  DBI::dbDisconnect(con)
  end_time <- Sys.time()
  elapsed <- round(difftime(end_time, start_time, units = "secs"), 2)
  cli::cli_alert_info("Done in {elapsed} seconds")

  x$data = active_tbl

  return(x)
}

#' @export
exclude_difficult_regions.data.frame = function(x, beds_exclude = c("hg19diff", "UCSC_unusual", "GRC_exclusions"), ...) {
  start_time <- Sys.time()
  cli::cli_alert_info("Now excluding 'difficult' regions...")

  for (bed in beds_exclude) {
    data(list = bed)

    bed_tbl = get(bed) %>%
      dplyr::select(chrom, start, end) %>%
      dplyr::mutate(start = start + 1) # convert from 0-based to 1-based

    x = x %>%
      dplyr::anti_join(bed_tbl, by = dplyr::join_by(CHROM == chrom, between(x$POS, y$start, y$end)))
  }

  end_time <- Sys.time()
  elapsed <- round(difftime(end_time, start_time, units = "secs"), 2)
  cli::cli_alert_info("Done in {elapsed} seconds")

  return(x)
}

#' @export
exclude_difficult_regions.tbl_df = function(x, ...) {
    exclude_difficult_regions.data.frame(x, ...)
}

