#' Connect to DuckDB
#' @param read_only Logical. If TRUE, the connection will be read-only.
#' @return A DuckDB connection object.
#' @export 
db_connect = function(read_only = FALSE) {
  con = DBI::dbConnect(duckdb::duckdb(), dbdir = "local.duckdb", read_only = read_only)

  limit = getOption("duckdb_max_memory")

  if (!is.null(limit)) {
    DBI::dbExecute(con, glue("SET memory_limit='{limit}'")) # Set the memory limit to 6GB
  }

  return(con)
}

copy_to_if_missing = function(con, df, name, temporary = FALSE) {
  if (!DBI::dbExistsTable(con, name)) {
    dplyr::copy_to(
      con,
      df,
      name = name,
      temporary = temporary,
      overwrite = FALSE
    )
  }

  dplyr::tbl(con, name)
}

sanitize_table_name_ = function(value, default = "summary_stats") {
  if (is.null(value) || length(value) == 0 || is.na(value) || !nzchar(value)) {
    return(default)
  }

  value = gsub("[^A-Za-z0-9]+", "_", value)
  value = gsub("^_+|_+$", "", value)
  value = tolower(value)

  if (!nzchar(value)) {
    default
  } else {
    value
  }
}

make_unique_table_name_ = function(stem, suffix = NULL, deterministic = FALSE) {
  table_stem = sanitize_table_name_(stem)

  if (!is.null(suffix)) {
    table_stem = sanitize_table_name_(paste(table_stem, suffix, sep = "_"), table_stem)
  }

  if (deterministic) {
    return(table_stem)
  }

  token = paste(sample(c(letters, 0:9), 8, replace = TRUE), collapse = "")
  paste(table_stem, token, sep = "_")
}

make_gwas_table_name_ = function(file_path = NULL, use_cache = FALSE) {
  stem = if (!is.null(file_path)) {
    tools::file_path_sans_ext(basename(file_path))
  } else {
    "summary_stats"
  }

  make_unique_table_name_(paste("summary_stats", stem, sep = "_"), deterministic = use_cache)
}

resolve_gwas_active_tbl_ = function(x, active_table = c("filtered_variants", "summary_stats")) {
  active_table = match.arg(active_table)

  if (active_table == "summary_stats") {
    return(dplyr::tbl(x$con, x$source_table_name))
  }

  x$data
}

materialize_gwas_tbl_ = function(x, tbl, suffix, update_source = FALSE) {
  output_table = make_unique_table_name_(x$table_name, suffix)

  tbl %>%
    dplyr::compute(name = output_table, temporary = FALSE, overwrite = TRUE)

  x$table_name = output_table
  x$data = dplyr::tbl(x$con, output_table)

  if (update_source) {
    x$source_table_name = output_table
  }

  invisible(x)
}
