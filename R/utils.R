# modified from here: https://stackoverflow.com/questions/27788968/how-would-one-check-the-system-memory-available-using-r-on-a-windows-machine
getFreeMemoryGB <- function() {
  osName <- Sys.info()[["sysname"]]
  if (osName == "Windows") {
    x <- system2("wmic", args =  "OS get FreePhysicalMemory /Value", stdout = TRUE)
    x <- x[grepl("FreePhysicalMemory", x)]
    x <- gsub("FreePhysicalMemory=", "", x, fixed = TRUE)
    x <- gsub("\r", "", x, fixed = TRUE)
    return(as.integer(x) / (1024 ^ 2))
  } else if (osName == 'Linux') {
    x <- system2('free', args='-k', stdout=TRUE)
    x <- strsplit(x[2], " +")[[1]][4]
    return(as.integer(x) / (1024 ^ 2))
  } else {
    stop("Only supported on Windows and Linux")
  }
}


compute_sample = function(g, rows = 100L) {

  df = DBI::dbGetQuery(
    g$con,
    glue("SELECT * FROM summary_stats USING SAMPLE {rows}")
  )

  df
}
