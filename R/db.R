db_connect = function(read_only = FALSE) {
  con = DBI::dbConnect(duckdb::duckdb(), dbdir = "local.duckdb", read_only = FALSE)

  limit = getOption("duckdb_max_memory")

  DBI::dbExecute(con, glue("SET memory_limit='{limit}'")) # Set the memory limit to 6GB

  return(con)
}
