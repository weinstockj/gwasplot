#' @title Filter variants
#' @description Filter variants in a GWASFormatter object based on a whitelist or blacklist of variants.
#' @param x A GWASFormatter object.
#' @param subset A file path to a whitelist of variants in Parquet format. Only variants in this file will be kept.
#' @export
filter_variants = function(x, subset = NULL, exclude = NULL, ...) {
  UseMethod("filter_variants")
}

#' @export
filter_variants.GWASFormatter = function(
  x,
  subset = NULL,
  exclude = NULL,
  ...
) {
  if (!is.null(subset) && !is.null(exclude)) {
    cli::cli_abort("You cannot specify both subset and exclude arguments.")
  }

  start_time <- Sys.time()
  cli::cli_process_start("Filtering variants")

  if (!is.null(subset)) {
    whitelist = dplyr::tbl(x$con, glue::glue("read_parquet('{subset}')"))

    materialize_gwas_tbl_(
      x,
      dplyr::semi_join(x$data, whitelist, by = c("chrom", "pos", "ref", "alt")),
      suffix = "filtered"
    )
  }

  if (!is.null(exclude)) {
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
exclude_difficult_regions.GWASFormatter = function(
  x,
  beds_exclude = c("hg19diff", "UCSC_unusual", "GRC_exclusions"),
  active_table = "filtered_variants",
  ...
) {
  start_time <- Sys.time()
  cli::cli_alert_info("Now excluding 'difficult' regions...")

  active_table = match.arg(
    active_table,
    c("filtered_variants", "summary_stats")
  )
  stopifnot(all(
    beds_exclude %in%
      c("hg19diff", "UCSC_unusual", "GRC_exclusions", "GIAB_difficult_regions")
  ))

  con = x$con %||% db_connect()

  active_tbl = resolve_gwas_active_tbl_(x, active_table)

  for (bed in beds_exclude) {
    data(list = bed)

    bed_tbl = copy_to_if_missing(
      con,
      get(bed),
      name = bed,
      temporary = FALSE
    ) %>%
      dplyr::select(chrom, start, end) %>%
      dplyr::mutate(start = start + 1) # convert from 0-based to 1-based

    active_tbl = dplyr::anti_join(
      active_tbl,
      bed_tbl,
      by = dplyr::join_by(chrom, between(x$POS, y$start, y$end))
    )
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
exclude_difficult_regions.data.frame = function(
  x,
  beds_exclude = c("hg19diff", "UCSC_unusual", "GRC_exclusions"),
  ...
) {
  start_time <- Sys.time()
  cli::cli_alert_info("Now excluding 'difficult' regions...")

  for (bed in beds_exclude) {
    data(list = bed)

    bed_tbl = get(bed) %>%
      dplyr::select(chrom, start, end) %>%
      dplyr::mutate(start = start + 1) # convert from 0-based to 1-based

    x = x %>%
      dplyr::anti_join(
        bed_tbl,
        by = dplyr::join_by(CHROM == chrom, between(x$POS, y$start, y$end))
      )
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

# Identify aberrant "lonely peak" lead variants: a strong lead with no nearby
# supporting signal. Returns the subset of `df` rows that are aberrant leads.
detect_aberrant_leads_ = function(
  df,
  lead_pvalue_threshold,
  support_pvalue_threshold,
  window_kb
) {
  window_bp <- window_kb * 1000

  lead_rows <- which(df$PVALUE < lead_pvalue_threshold)
  if (length(lead_rows) == 0) {
    return(df[integer(0), , drop = FALSE])
  }

  chrom <- df$CHROM
  pos <- df$POS
  pvalue <- df$PVALUE
  row_index <- seq_len(nrow(df))

  is_aberrant <- vapply(
    lead_rows,
    function(r) {
      # Support = any *other* variant on the same chromosome within the window
      # that itself reaches the support threshold.
      supported <- chrom == chrom[r] &
        abs(pos - pos[r]) <= window_bp &
        pvalue <= support_pvalue_threshold &
        row_index != r
      !any(supported, na.rm = TRUE)
    },
    logical(1)
  )

  df[lead_rows[is_aberrant], , drop = FALSE]
}

# Log how many aberrant leads were found before dropping them.
report_aberrant_leads_ = function(
  n_drop,
  lead_pvalue_threshold,
  support_pvalue_threshold,
  window_kb
) {
  if (n_drop == 0) {
    cli::cli_alert_success(
      "No aberrant lead variants found (lead < {lead_pvalue_threshold}; a supporting variant <= {support_pvalue_threshold} exists within {window_kb} kb)."
    )
  } else {
    cli::cli_alert_warning(
      "Dropping {n_drop} aberrant lead variant{?s}: PVALUE < {lead_pvalue_threshold} with no other variant <= {support_pvalue_threshold} within {window_kb} kb."
    )
  }
}

#' @title Exclude aberrant "lonely peak" loci
#' @description
#' Flag and drop genome-wide significant lead variants that have no supporting
#' signal from nearby variants. A genuine association is usually accompanied by
#' correlated (LD) variants with elevated signal, so a very strong lead with no
#' nearby support is suspicious (e.g. a genotyping artifact or a spurious
#' rare-variant singleton). A lead is treated as aberrant when its p-value is
#' below `lead_pvalue_threshold` and no other variant within `window_kb`
#' kilobases on the same chromosome reaches `support_pvalue_threshold`. The
#' number of variants that would be dropped is logged before they are removed.
#' @param x A `GWASFormatter` object, `data.frame`, or `tibble`.
#' @param lead_pvalue_threshold Lead variants with `PVALUE` strictly below this
#'   are examined. Default `5e-10`.
#' @param support_pvalue_threshold A lead is kept (supported) if any other
#'   variant within the window reaches this p-value. Default `5e-5`.
#' @param window_kb Distance in kilobases from the lead, on the same chromosome,
#'   searched for a supporting variant. Default `500`.
#' @param ... Additional arguments passed to methods.
#' @return The input with aberrant lead variants removed. `GWASFormatter`
#'   returns the modified object; eager methods return the filtered data.
#' @export
exclude_aberrant_pvalue_loci = function(x, ...) {
  UseMethod("exclude_aberrant_pvalue_loci")
}

#' @export
exclude_aberrant_pvalue_loci.data.frame = function(
  x,
  lead_pvalue_threshold = 5e-10,
  support_pvalue_threshold = 5e-5,
  window_kb = 500,
  ...
) {
  work <- x
  work$.aberrant_row_id <- seq_len(nrow(work))

  aberrant <- detect_aberrant_leads_(
    work,
    lead_pvalue_threshold,
    support_pvalue_threshold,
    window_kb
  )

  n_drop <- nrow(aberrant)
  report_aberrant_leads_(
    n_drop,
    lead_pvalue_threshold,
    support_pvalue_threshold,
    window_kb
  )

  if (n_drop == 0) {
    return(x)
  }

  x[-aberrant$.aberrant_row_id, , drop = FALSE]
}

#' @export
exclude_aberrant_pvalue_loci.tbl_df = function(x, ...) {
  exclude_aberrant_pvalue_loci.data.frame(x, ...)
}

#' @export
exclude_aberrant_pvalue_loci.GWASFormatter = function(
  x,
  lead_pvalue_threshold = 5e-10,
  support_pvalue_threshold = 5e-5,
  window_kb = 500,
  ...
) {
  start_time <- Sys.time()
  cli::cli_alert_info("Scanning for aberrant 'lonely peak' loci...")

  # Only variants reaching the support threshold can support a strong lead, so
  # collecting this subset is sufficient to decide which leads are aberrant. The
  # leads themselves satisfy the filter (normally lead < support); take the
  # looser of the two thresholds to stay correct even if they are reversed.
  collect_threshold <- max(lead_pvalue_threshold, support_pvalue_threshold)

  candidates <- x$data %>%
    dplyr::filter(PVALUE <= collect_threshold) %>%
    dplyr::select(dplyr::any_of(c("CHROM", "POS", "ID", "PVALUE"))) %>%
    dplyr::collect()

  aberrant <- detect_aberrant_leads_(
    candidates,
    lead_pvalue_threshold,
    support_pvalue_threshold,
    window_kb
  )

  n_drop <- nrow(aberrant)
  report_aberrant_leads_(
    n_drop,
    lead_pvalue_threshold,
    support_pvalue_threshold,
    window_kb
  )

  if (n_drop == 0) {
    return(x)
  }

  if (!"ID" %in% names(aberrant)) {
    cli::cli_abort(
      "Cannot drop aberrant variants: input has no {.field ID} column."
    )
  }

  drop_ids <- aberrant$ID
  materialize_gwas_tbl_(
    x,
    x$data %>% dplyr::filter(!ID %in% drop_ids),
    suffix = "aberrant_excluded"
  )

  elapsed <- round(difftime(Sys.time(), start_time, units = "secs"), 2)
  cli::cli_alert_info("Done in {elapsed} seconds")

  return(x)
}
