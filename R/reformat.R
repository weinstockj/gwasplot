#' Reformat GWAS summary statistics from regenie or SAIGE
#'
#' This function takes a file path to a parquet or CSV file of GWAS summary statistics
#' output from regenie or SAIGE and reformats it to have a standard set of columns.
#'
#' @param file_path Path to the GWAS summary statistics file (parquet or CSV).
#' @param use_cache Logical. If TRUE, uses the cached table 'summary_stats' instead of reading the file again.
#' @return An R6 object of class `GWASFormatter` containing the reformatted summary statistics.
#' @export
reformat_summary_statistics <- function(file_path, use_cache = FALSE) {
  # Check if the file exists
  if (!file.exists(file_path)) {
    stop(paste("File not found:", file_path))
  }

  gwas = GWASFormatter$new(file_path, use_cache)

  if (!use_cache) {
    gwas$reformat()
  }

  return(gwas)
}

standard_format_names_ = function() {
  list(
    "saige" = c("CHR", "POS", "Allele1", "Allele2", "AC_Allele2", "AF_Allele2",
                "N", "BETA", "SE", "p.value"),
    "regenie" = c("CHROM", "POS", "ALLELE0", "ALLELE1", "A1FREQ", "BETA", "LOG10P")
  )
}



detect_format = function(names) {

  stopifnot(
    is.character(names),
    length(names) > 0,
    length(names) < 100
  )

  standard_names = standard_format_names_()

  all_contained = purrr::map_lgl(
      standard_names,
      ~all(.x %in% names)
    ) %>%
    purrr::set_names(names(standard_names))

  if (sum(all_contained) == 0) {
    stop("No standard format detected. Please provide a file in SAIGE or regenie format.")
  }

  if (sum(all_contained) > 1) {
    stop("Multiple formats detected. Please provide a file in SAIGE or regenie format.")
  }

  if (all_contained["regenie"]) {
    log_info("Detected regenie format")
    return("regenie")
  }

  if (all_contained["saige"]) {
    log_info("Detected saige format")
    return("saige")
  }

  stop("Unknown format. Please provide a file in SAIGE or regenie format.")
}


reformat_lookup = function() {
  list(
    regenie = c(
      "CHROM",
      "POS",
      "REF" = "ALLELE0",
      "ALT" = "ALLELE1",
      "AF_ALT" = "A1FREQ",
      "BETA",
      "LOG10P"
    ),
    saige = c(
      "CHROM" = "CHR",
      "POS",
      "REF" = "Allele1",
      "ALT" = "Allele2",
      "AF_ALT" = "AF_Allele2",
      "BETA",
      "PVALUE" = "p.value"
    )
  )
}

reformat_names = function(ds, format) {

  lookup = reformat_lookup()[[format]]

  ds %>%
    dplyr::select(all_of(lookup))
}

possibly_undo_log10p = function(ds, format) {

  formats_to_undo = c("regenie")

  if (format %in% formats_to_undo) {

    log_info(glue("Converting LOG10PVALUE to PVALUE for {format} format"))

    ds %>%
      dplyr::mutate(
        PVALUE = 10 ^ (-LOG10P)
      )
  } else {
    ds
  }
}

add_ID = function(ds) {
  ds %>%
    dplyr::mutate(
      ID = str_c(CHROM, "_", POS, "_", REF, "_", ALT)
    )
}

possibly_replace_chr23 = function(ds, format) {
  if (format == "regenie") {
    ds %>%
      dplyr::mutate(
        CHROM = dplyr::case_when(
          CHROM == "chr23" ~ "chrX",
          TRUE ~ CHROM
        )
      )
  } else {
    ds
  }
}

GWASFormatter <- R6::R6Class(
  "GWASFormatter",
  public = list(
    con = NULL,
    data = NULL,
    data_names = NULL,
    detected_format = NULL,
    file_path = NULL,
    initialize = function(file_path, use_cache = FALSE) {
      self$file_path = file_path
      self$con = db_connect()

      file_ext <- tolower(tools::file_ext(file_path))

      if (use_cache) {
        cli::cli_alert_info("Using cached table 'summary_stats'")
        self$data = tbl(self$con, "summary_stats")
      } else {
          if (file_ext == "parquet") {
            tryCatch({
              # Read parquet file with optimized column types
              self$data = tbl(self$con, glue("
                SELECT 
                  CAST(CASE 
                    WHEN CHROM LIKE 'chr%' THEN CHROM 
                    ELSE CONCAT('chr', CHROM) 
                  END AS VARCHAR) AS CHROM,
                  CAST(POS AS INTEGER) AS POS,
                  * EXCLUDE(CHROM, POS)
                FROM read_parquet('{file_path}')
              "))
            }, error = function(e) {
              stop(paste("Error reading parquet file:", e$message))
            })
        } else if (file_ext == "csv") {
            tryCatch({
              # Read CSV file with optimized column types
              self$data = tbl(self$con, glue("
                SELECT 
                  CAST(CASE 
                    WHEN CHROM LIKE 'chr%' THEN CHROM 
                    ELSE CONCAT('chr', CHROM) 
                  END AS VARCHAR) AS CHROM,
                  CAST(POS AS INTEGER) AS POS,
                  * EXCLUDE(CHROM, POS)
                FROM read_csv('{file_path}')
              "))
            }, error = function(e) {
              stop(paste("Error reading CSV file:", e$message))
            })
        } else {
            stop("Unsupported file type. Please provide a parquet or CSV file.")
        }
      }

      self$data_names = self$data %>%
        head(n = 2) %>%
        dplyr::collect(.) %>%
        names()

      if(!use_cache) {
        self$detected_format = detect_format(self$data_names)
      } else {
        self$detected_format = "original format unknown"
      }
    },
    reformat = function(use_cache = FALSE) {
      self$data = self$data %>%
        reformat_names(self$detected_format) %>%
        possibly_undo_log10p(self$detected_format) %>%
        possibly_replace_chr23(self$detected_format) %>%
        add_ID() %>%
        dplyr::compute(temporary = FALSE, overwrite = TRUE, name = "summary_stats")
    },
    kill = function() {
      DBI::dbDisconnect(self$con)
    },
    print = function(...) {
      cat("GWAS object\n")
      cat("File path:", self$file_path, "\n")
      cat("Detected format:", self$detected_format, "\n")
      cat("Data names:", paste(self$data_names, collapse = ", "), "\n")
      cat("Data preview:\n")
      print(head(self$data, n = 5) %>% dplyr::collect(.))
    }
  )
)
