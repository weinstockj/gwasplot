# Resolve the output format from an explicit argument or the file extension.
resolve_write_format_ = function(output, format = NULL) {
  if (!is.null(format)) {
    format = match.arg(tolower(format), c("parquet", "csv"))
    return(format)
  }

  ext = tolower(tools::file_ext(output))
  if (ext %in% c("parquet", "pq")) {
    return("parquet")
  }
  if (ext %in% c("csv", "txt", "tsv")) {
    return("csv")
  }

  "parquet"
}

# Abort if the output exists and the caller has not opted into overwriting.
check_overwrite_ = function(output, overwrite) {
  if (!overwrite && file.exists(output)) {
    cli::cli_abort(c(
      "Output file {.path {output}} already exists.",
      "i" = "Set {.code overwrite = TRUE} to replace it."
    ))
  }
}

# Run a DuckDB COPY from a SQL source (a SELECT statement or table name) to disk.
copy_to_file_ = function(con, source_sql, output, format, compression) {
  options_sql = if (format == "parquet") {
    sprintf("(FORMAT PARQUET, COMPRESSION %s)", toupper(compression))
  } else {
    "(FORMAT CSV, HEADER)"
  }

  DBI::dbExecute(
    con,
    sprintf(
      "COPY (%s) TO '%s' %s",
      source_sql,
      gsub("\\\\", "/", output),
      options_sql
    )
  )
}

#' @title Write summary statistics to disk
#' @description Write a (possibly filtered) summary-statistic object to a Parquet
#'   or CSV file. For a `GWASFormatter`, the write streams directly from DuckDB
#'   without collecting into R, so it scales to large data and captures whatever
#'   filtering has been applied (e.g. via [filter_variants()] or
#'   [exclude_difficult_regions()]). `data.frame` and `tibble` inputs are written
#'   through a temporary in-memory DuckDB connection.
#' @param x A `GWASFormatter` object, `data.frame`, or `tibble`.
#' @param output Output file path.
#' @param format Output format, one of `"parquet"` or `"csv"`. Default `NULL`,
#'   which infers the format from the extension of `output` (falling back to
#'   `"parquet"`).
#' @param compression Compression codec for Parquet output (e.g. `"snappy"`,
#'   `"zstd"`, `"gzip"`, `"uncompressed"`). Ignored for CSV. Default `"snappy"`.
#' @param overwrite Logical. If `FALSE` (default) and `output` already exists,
#'   abort rather than replace it.
#' @param ... Additional arguments passed to methods.
#' @return The path to the written file, invisibly.
#' @export
write_summary_statistics = function(
  x,
  output,
  format = NULL,
  compression = "snappy",
  overwrite = FALSE,
  ...
) {
  UseMethod("write_summary_statistics")
}

#' @export
write_summary_statistics.GWASFormatter = function(
  x,
  output,
  format = NULL,
  compression = "snappy",
  overwrite = FALSE,
  ...
) {
  format = resolve_write_format_(output, format)
  check_overwrite_(output, overwrite)

  cli::cli_process_start("Writing summary statistics to {.path {output}}")

  # Render the current lazy query so any applied filters are captured.
  source_sql = dbplyr::sql_render(x$data)
  copy_to_file_(x$con, source_sql, output, format, compression)

  cli::cli_process_done()
  cli::cli_alert_success("Wrote {format} to {.path {output}}")

  invisible(output)
}

#' @export
write_summary_statistics.data.frame = function(
  x,
  output,
  format = NULL,
  compression = "snappy",
  overwrite = FALSE,
  ...
) {
  format = resolve_write_format_(output, format)
  check_overwrite_(output, overwrite)

  cli::cli_process_start("Writing summary statistics to {.path {output}}")

  # Use an isolated in-memory DuckDB connection so we never touch local.duckdb.
  con = DBI::dbConnect(duckdb::duckdb())
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE)

  tmp_name = "gwasplot_write_tmp"
  duckdb::duckdb_register(con, tmp_name, x)
  on.exit(duckdb::duckdb_unregister(con, tmp_name), add = TRUE)

  copy_to_file_(
    con,
    sprintf("SELECT * FROM %s", tmp_name),
    output,
    format,
    compression
  )

  cli::cli_process_done()
  cli::cli_alert_success("Wrote {format} to {.path {output}}")

  invisible(output)
}

#' @export
write_summary_statistics.tbl_df = function(
  x,
  output,
  format = NULL,
  compression = "snappy",
  overwrite = FALSE,
  ...
) {
  write_summary_statistics.data.frame(
    x,
    output,
    format = format,
    compression = compression,
    overwrite = overwrite,
    ...
  )
}
