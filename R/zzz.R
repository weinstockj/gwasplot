.onLoad <- function(libname, pkgname) {
  # Set the default value for the duckdb_max_memory option

  free_RAM = getFreeMemoryGB()

  set_RAM = floor(free_RAM * 0.8) # Set to 80% of free RAM
  #set_RAM = 10
  cli::cli_alert_info(glue::glue("Setting duckdb_max_memory to {set_RAM}GB, using 80% of available system memory."))
  cli::cli_alert_info("Change this with options(duckdb_max_memory = 'XGB')")
  options(duckdb_max_memory = glue("{set_RAM}GB")) # Or any other reasonable default

  # You can also perform other package initialization tasks here
}
