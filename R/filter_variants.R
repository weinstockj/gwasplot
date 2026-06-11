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

      materialize_gwas_tbl_(
        x,
        dplyr::semi_join(x$data, whitelist, by = c("chrom", "pos", "ref", "alt")),
        suffix = "filtered"
      )
    }

    if(!is.null(exclude)) {
        blacklist = dplyr::tbl(x$con, glue::glue("read_parquet('{exclude}')"))

      materialize_gwas_tbl_(
        x,
        dplyr::anti_join(x$data, blacklist, by = c("chrom", "pos", "ref", "alt")),
        suffix = "filtered"
      )
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

#' @export
exclude_difficult_regions.GWASFormatter = function(x, beds_exclude = c("hg19diff", "UCSC_unusual", "GRC_exclusions"), active_table = "filtered_variants",  ...) {
  start_time <- Sys.time()
  cli::cli_alert_info("Now excluding 'difficult' regions...")

  active_table = match.arg(active_table, c("filtered_variants", "summary_stats"))
  stopifnot(all(beds_exclude %in% c("hg19diff", "UCSC_unusual", "GRC_exclusions", "GIAB_difficult_regions")))
  
  con = x$con %||% db_connect()

  active_tbl = resolve_gwas_active_tbl_(x, active_table)

  for (bed in beds_exclude) {
    data(list = bed)

    bed_tbl = copy_to_if_missing(con, get(bed), name = bed, temporary = FALSE) %>%
      dplyr::select(chrom, start, end) %>%
      dplyr::mutate(start = start + 1) # convert from 0-based to 1-based

    active_tbl = dplyr::anti_join(active_tbl, bed_tbl, by = dplyr::join_by(chrom, between(x$POS, y$start, y$end))) 
  }

    materialize_gwas_tbl_(x, active_tbl, suffix = "excluded")

  if (is.null(x$con)) {
    DBI::dbDisconnect(con)
  }
  end_time <- Sys.time()
  elapsed <- round(difftime(end_time, start_time, units = "secs"), 2)
  cli::cli_alert_info("Done in {elapsed} seconds")

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

