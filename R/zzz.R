.onLoad <- function(libname, pkgname) {
  # Set the default value for the duckdb_max_memory option

  free_RAM = getFreeMemoryGB()

  set_RAM = floor(free_RAM * 0.8) # Set to 80% of free RAM
  set_RAM = 10
  message(glue::glue("Setting duckdb_max_memory to {set_RAM}GB, using 80% of available system memory."))
  options(duckdb_max_memory = glue("{set_RAM}GB")) # Or any other reasonable default

  # You can also perform other package initialization tasks here
}
