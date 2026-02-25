# Helper function for null coalescing operator
`%||%` <- function(x, y) if (is.null(x)) y else x

#' Pull top hits from a GWASFormatter or a data.frame/tibble
#'
#' @param x A GWASFormatter object, data.frame, or tibble.
#' @param threshold The p-value threshold to filter the top hits. Default is 5e-8.
#' @param ... Additional arguments (unused).
#' @return For GWASFormatter, a tibble of filtered hits; for data.frame/tibble, a filtered data.frame/tibble.
#' @export
select_top_hits <- function(x, threshold = 5e-8, ...) {
  UseMethod("select_top_hits")
}

#' @describeIn select_top_hits Method for GWASFormatter objects
#' @export
select_top_hits.GWASFormatter <- function(x, threshold = 5e-8, ...) {
  df = x$data %>%
    dplyr::filter(PVALUE < threshold) %>%
    dplyr::collect(.)

  # ot = purrr::map(df$ID, ~possibly_query(stringr::str_remove(.x, "chr"))) %>%
  #   purrr::compact(.) %>%
  #   dplyr::bind_rows(.) %>%
  #   tibble::as_tibble(.) %>%
  #   dplyr::rename(ID = id) %>%
  #   dplyr::mutate(ID = glue::glue("chr{ID}"))
  #
  # df %>%
  #   dplyr::left_join(ot)
  return(df)
}

#' @describeIn select_top_hits Method for data.frame/tibble objects
#' @export
select_top_hits.data.frame <- function(x, threshold = 5e-8, ...) {
  x %>%
    dplyr::filter(PVALUE < threshold)
}

#' Find the nearest gene for variants
#'
#' @param x A gwas object or tibble containing variant data.
#' @param ... Additional arguments passed to methods.
#' @return The input object with gene annotations added.
#' @export
find_nearest_gene = function(x, threshold = 1e5, pvalue_threshold = NULL, ...) {
  UseMethod("find_nearest_gene")
}

#' @describeIn find_nearest_gene Find the nearest gene for each variant in a gwas object
#' @param threshold The distance threshold to consider a gene as nearest. Default is 1e5.
#' @export
find_nearest_gene.GWASFormatter = function(
  x,
  threshold = 1e5,
  pvalue_threshold = NULL,
  ...
) {
  if (!is.null(pvalue_threshold)) {
    con <- x$con

    # Materialise significant variants inside x$con — no data leaves DuckDB
    sig_query <- x$data %>% dplyr::filter(PVALUE < pvalue_threshold)
    if (!"ID" %in% colnames(sig_query)) {
      sig_query <- dplyr::mutate(sig_query, ID = paste(CHROM, POS, sep = "_"))
    }

    sig_tbl <- sig_query %>%
      dplyr::compute(name = "__sig_variants__", overwrite = TRUE)

    n_sig <- dplyr::count(sig_tbl) %>% dplyr::pull(n)
    if (n_sig == 0L) {
      cli::cli_inform("No variants below pvalue_threshold {pvalue_threshold}.")
      return(x)
    }

    # Load gene reference if not already in this connection
    if (!"human_genes" %in% DBI::dbListTables(con)) {
      dplyr::copy_to(
        con,
        human_genes,
        name = "human_genes",
        temporary = FALSE,
        overwrite = TRUE
      )
    }

    # Gene intervals scoped to this query
    DBI::dbExecute(
      con,
      glue::glue(
        "CREATE OR REPLACE TABLE __gene_intervals__ AS
       SELECT gene_id, gene_name, chrom,
         start - {format(threshold, scientific = FALSE)} AS expanded_start,
         \"end\" + {format(threshold, scientific = FALSE)} AS expanded_end,
         start, \"end\"
       FROM human_genes
       WHERE gene_biotype = 'protein_coding' AND gene_name IS NOT NULL"
      )
    )

    # Range join on the small significant subset — pure DuckDB
    DBI::dbExecute(
      con,
      "CREATE OR REPLACE TABLE __locus_gene_annot__ AS
       WITH ng AS (
         SELECT t.CHROM, t.POS,
           g.gene_id, g.gene_name,
           CASE
             WHEN t.POS >= g.start AND t.POS <= g.\"end\" THEN 0
             WHEN t.POS < g.start THEN g.start - t.POS
             ELSE t.POS - g.\"end\"
           END AS distance
         FROM __sig_variants__ t
         JOIN __gene_intervals__ g
           ON t.CHROM = g.chrom
          AND t.POS >= g.expanded_start
          AND t.POS <= g.expanded_end
       ),
       rk AS (
         SELECT *, ROW_NUMBER() OVER (PARTITION BY ID ORDER BY distance) AS rn
         FROM ng
       )
       SELECT CHROM, POS, gene_id, gene_name,
              CAST(distance AS INTEGER) AS distance
       FROM rk WHERE rn = 1"
    )

    x$data <- x$data %>%
      dplyr::left_join(
        dplyr::tbl(con, "__locus_gene_annot__"),
        by = c("CHROM", "POS")
      )
    return(x)
  }

  # Start timing
  start_time <- Sys.time()
  cli::cli_alert_info("Starting gene annotation...")

  con = db_connect()

  DBI::dbExecute(con, "PRAGMA max_temp_directory_size = '30GB'")

  # Load human genes data
  cli::cli_progress_step("Loading gene reference data")
  dplyr::copy_to(
    con,
    human_genes,
    name = "human_genes",
    temporary = FALSE,
    overwrite = TRUE
  )

  # Create spatial index for genes
  cli::cli_progress_step("Creating optimized gene intervals\n")
  sql_index = glue(
    "
  SELECT 
    gene_id, 
    gene_name,
    gene_biotype,
    chrom,
    start - {format(threshold, scientific = FALSE)} AS expanded_start,
    g.\"end\" + {format(threshold, scientific = FALSE)} AS expanded_end,
    start,
    g.\"end\"
  FROM human_genes g
  WHERE gene_biotype = 'protein_coding' AND gene_name IS NOT NULL
  "
  )

  intervals = dplyr::tbl(con, dplyr::sql(sql_index)) %>%
    dplyr::compute(temporary = FALSE, overwrite = TRUE, name = "gene_intervals")

  DBI::dbExecute(
    con,
    "CREATE INDEX IF NOT EXISTS chrom_start_end ON gene_intervals (chrom, expanded_start, expanded_end)"
  )
  DBI::dbExecute(
    con,
    "CREATE INDEX IF NOT EXISTS chrom_pos ON summary_stats (chrom, POS)"
  )

  cli::cli_progress_step("Finding nearest genes")
  sql = "
  CREATE OR REPLACE TABLE summary_stats_annotated AS
  WITH NearestGenes AS (
    SELECT
      t.*,
      g.gene_id,
      g.gene_name,
      CASE
        -- Variant inside gene
        WHEN t.POS >= g.start AND t.POS <= g.end THEN 0
        -- Variant upstream of gene
        WHEN t.POS < g.start THEN g.start - t.POS
        -- Variant downstream of gene
        ELSE t.POS - g.end
      END AS distance  
    FROM summary_stats t
    -- Efficient range join instead of cross join
    JOIN gene_intervals g ON 
      t.chrom = g.chrom AND 
      t.POS >= g.expanded_start AND 
      t.POS <= g.expanded_end
  ),
  RankedGenes AS (
    SELECT 
      *,
      ROW_NUMBER() OVER (PARTITION BY ID ORDER BY distance) AS rn
    FROM NearestGenes
  )
  SELECT
    CHROM,
    POS,
    ID,
    gene_id,
    gene_name,
    distance
  FROM RankedGenes
  WHERE rn = 1
  "

  cli::cli_progress_step("Updating gwas object with gene annotations\n")
  DBI::dbExecute(con, sql)

  x$data = dplyr::tbl(con, "summary_stats_annotated") %>%
    dplyr::select(ID, gene_id, gene_name, distance) %>%
    dplyr::inner_join(x$data, by = "ID", copy = TRUE)

  # End timing and report
  end_time <- Sys.time()
  elapsed <- round(difftime(end_time, start_time, units = "secs"), 2)
  cli::cli_alert_success("Gene annotation completed in {elapsed} seconds")

  return(x)
}

#' @describeIn find_nearest_gene Find the nearest gene for each variant in a tibble
#' @param chrom_col Name of the column containing chromosome information. Default is "CHROM".
#' @param pos_col Name of the column containing position information. Default is "POS".
#' @param id_col Name of the column containing variant IDs. Default is "ID".
#' @export
find_nearest_gene.tbl_df = function(
  x,
  threshold = 1e5,
  chrom_col = "CHROM",
  pos_col = "POS",
  id_col = "ID",
  pvalue_threshold = NULL,
  ...
) {
  find_nearest_gene.data.frame(
    x,
    threshold,
    chrom_col,
    pos_col,
    id_col,
    pvalue_threshold,
    ...
  )
}

#' @export
find_nearest_gene.data.frame = function(
  x,
  threshold = 1e5,
  chrom_col = "CHROM",
  pos_col = "POS",
  id_col = "ID",
  pvalue_threshold = NULL,
  ...
) {
  if (!is.null(pvalue_threshold)) {
    x <- dplyr::select(x, -dplyr::any_of(c("gene_name", "gene_id", "distance")))
    sig <- dplyr::filter(x, PVALUE < pvalue_threshold)
    if (nrow(sig) == 0) {
      cli::cli_inform("No variants below pvalue_threshold {pvalue_threshold}.")
      return(dplyr::mutate(
        x,
        gene_name = NA_character_,
        gene_id = NA_character_,
        distance = NA_integer_
      ))
    }
    added_id <- !id_col %in% colnames(sig)
    if (added_id) {
      sig <- dplyr::mutate(sig, !!id_col := paste(CHROM, POS, sep = "_"))
    }

    sig_ann <- find_nearest_gene(
      sig,
      threshold = threshold,
      chrom_col = chrom_col,
      pos_col = pos_col,
      id_col = id_col
    )
    gene_cols <- dplyr::select(
      sig_ann,
      CHROM,
      POS,
      gene_name,
      gene_id,
      distance
    )
    return(dplyr::left_join(x, gene_cols, by = c("CHROM", "POS")))
  }

  start_time <- Sys.time()
  cli::cli_alert_info("Starting gene annotation...")

  renamed_df <- x %>%
    dplyr::rename(
      CHROM = !!rlang::sym(chrom_col),
      POS = !!rlang::sym(pos_col),
      ID = !!rlang::sym(id_col)
    )

  # In-memory DuckDB — same fast range-join as GWASFormatter, no file written
  con <- DBI::dbConnect(duckdb::duckdb())
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE)

  cli::cli_progress_step("Loading data into in-memory database")
  DBI::dbWriteTable(con, "summary_stats", renamed_df, overwrite = TRUE)
  dplyr::copy_to(
    con,
    human_genes,
    name = "human_genes",
    temporary = FALSE,
    overwrite = TRUE
  )

  cli::cli_progress_step("Creating optimised gene intervals\n")
  DBI::dbExecute(
    con,
    glue::glue(
      "
    CREATE TABLE gene_intervals AS
    SELECT
      gene_id, gene_name, gene_biotype, chrom,
      start - {format(threshold, scientific = FALSE)} AS expanded_start,
      \"end\" + {format(threshold, scientific = FALSE)}  AS expanded_end,
      start, \"end\"
    FROM human_genes
    WHERE gene_biotype = 'protein_coding' AND gene_name IS NOT NULL
  "
    )
  )
  DBI::dbExecute(
    con,
    "CREATE INDEX idx_gene ON gene_intervals (chrom, expanded_start, expanded_end)"
  )
  DBI::dbExecute(con, "CREATE INDEX idx_snp  ON summary_stats  (CHROM, POS)")

  cli::cli_progress_step("Finding nearest genes")
  result <- DBI::dbGetQuery(
    con,
    "
    WITH NearestGenes AS (
      SELECT
        t.ID,
        g.gene_id,
        g.gene_name,
        CASE
          WHEN t.POS >= g.start AND t.POS <= g.\"end\" THEN 0
          WHEN t.POS < g.start                          THEN g.start - t.POS
          ELSE t.POS - g.\"end\"
        END AS distance
      FROM summary_stats t
      JOIN gene_intervals g
        ON t.CHROM = g.chrom
       AND t.POS  >= g.expanded_start
       AND t.POS  <= g.expanded_end
    ),
    RankedGenes AS (
      SELECT *, ROW_NUMBER() OVER (PARTITION BY ID ORDER BY distance) AS rn
      FROM NearestGenes
    )
    SELECT ID, gene_id, gene_name, distance FROM RankedGenes WHERE rn = 1
  "
  ) %>%
    tibble::as_tibble()

  out <- renamed_df %>%
    dplyr::left_join(result, by = "ID") %>%
    dplyr::mutate(
      distance = dplyr::if_else(
        is.na(distance),
        NA_integer_,
        as.integer(distance)
      )
    )

  end_time <- Sys.time()
  elapsed <- round(difftime(end_time, start_time, units = "secs"), 2)
  cli::cli_alert_success("Gene annotation completed in {elapsed} seconds")
  return(out)
}

#' Annotate data with centromere information
#'
#' @param x A data frame or tibble containing variant data.
#' @param ... Additional arguments passed to methods.
#' @return A data frame with the centromere information.
#' @export
annotate_with_centromere = function(x, ...) {
  UseMethod("annotate_with_centromere")
}

#' @describeIn annotate_with_centromere Annotate a data frame or tibble with centromere information
#' @param chrom_col Name of the column containing chromosome information. Default is "CHROM".
#' @param pos_col Name of the column containing position information. Default is "POS".
#' @export
annotate_with_centromere.data.frame = function(
  x,
  chrom_col = "CHROM",
  pos_col = "POS",
  ...
) {
  # Start timing
  start_time <- Sys.time()
  cli::cli_progress_step("Starting centromere annotation for data frame...")

  # Rename columns if needed
  if (chrom_col != "CHROM" || pos_col != "POS") {
    top_hits <- x %>%
      dplyr::rename(
        CHROM = !!rlang::sym(chrom_col),
        POS = !!rlang::sym(pos_col)
      )
  } else {
    top_hits <- x
  }

  # Prepare ideogram data with centromere information
  ideogram_data <- ideogram %>%
    dplyr::mutate(
      in_centromere = ifelse(stain == "acen", TRUE, FALSE)
    )

  # Perform annotation using dplyr joins
  cli::cli_progress_step("Finding variants in centromere regions")
  result <- top_hits %>%
    # Left join with ideogram to find regions that variants fall into
    dplyr::left_join(
      ideogram_data %>%
        dplyr::select(CHROM = chrom, start, end, name, in_centromere),
      by = dplyr::join_by(CHROM, between(POS, start, end))
    ) %>%
    dplyr::mutate(
      in_centromere = dplyr::if_else(is.na(in_centromere), FALSE, in_centromere)
    )

  # Restore original column names if they were renamed
  if (chrom_col != "CHROM") {
    result <- result %>%
      dplyr::rename(
        !!rlang::sym(chrom_col) := chrom
      )
  }

  # End timing and report
  end_time <- Sys.time()
  elapsed <- round(difftime(end_time, start_time, units = "secs"), 2)
  cli::cli_alert_success("Centromere annotation completed in {elapsed} seconds")

  return(result)
}

#' @describeIn annotate_with_centromere Alias for the data.frame method
#' @export
annotate_with_centromere.tbl_df = function(x, ...) {
  annotate_with_centromere.data.frame(x, ...)
}

#' @describeIn annotate_with_centromere Centromere annotation method for GWASFormatter objects
#' @export
annotate_with_centromere.GWASFormatter = function(x, ...) {
  con = db_connect()

  dplyr::copy_to(
    con,
    ideogram %>%
      dplyr::mutate(
        in_centromere = ifelse(stain == "acen", TRUE, FALSE)
      ),
    name = "ideogram",
    temporary = FALSE,
    overwrite = TRUE
  )

  cen_tbl = dplyr::tbl(con, "ideogram")

  x$data = x$data %>%
    dplyr::left_join(
      cen_tbl %>%
        dplyr::select(chrom, start, end, in_centromere),
      by = dplyr::join_by(chrom, dplyr::between(POS, start, end))
    ) %>%
    dplyr::mutate(
      in_centromere = dplyr::case_when(
        is.na(in_centromere) ~ FALSE,
        TRUE ~ in_centromere
      )
    )

  return(x)
}

#' Annotate data with CHIP gene information
#' @export
annotate_with_chip_genes = function(x, ...) {
  UseMethod("annotate_with_chip_genes")
}

#' @export
annotate_with_chip_genes.data.frame = function(x, ...) {
  if (!"gene_name" %in% names(x)) {
    cli::cli_abort(
      "x must contain a gene_name column; did you run find_nearest_gene?"
    )
  }

  x %>%
    dplyr::mutate(
      is_chip_gene = dplyr::case_when(
        gene_name %in% chip_genes ~ TRUE,
        TRUE ~ FALSE
      )
    )
}

#' Annotate top hits with CHIP gene information
#'
#' @param top_hits A data frame containing the top hits.
#' @return A data frame with the CHIP gene information.
#' @export
annotate_with_chip_genes.GWASFormatter = function(top_hits) {
  if (!"gene_name" %in% names(top_hits)) {
    cli::cli_abort(
      "top_hits must contain a gene_name column; did you run find_nearest_gene?"
    )
  }

  con = db_connect()

  dplyr::copy_to(
    con,
    top_hits,
    name = "top_hits",
    temporary = FALSE,
    overwrite = TRUE
  )

  chip_df = tibble::tibble(
    gene_name = chip_genes
  ) %>%
    dplyr::distinct(.) %>%
    dplyr::mutate(
      is_chip_gene = TRUE
    )

  dplyr::copy_to(
    con,
    chip_df,
    name = "chip_genes",
    temporary = FALSE,
    overwrite = TRUE
  )

  sql = glue(
    "
  SELECT
    t.*,
    c.is_chip_gene
  FROM top_hits t
  LEFT JOIN chip_genes c ON (t.gene_name = c.gene_name)
  "
  )

  DBI::dbGetQuery(con, sql) %>%
    tibble::as_tibble(.) %>%
    dplyr::mutate(
      is_chip_gene = dplyr::case_when(
        is.na(is_chip_gene) ~ FALSE,
        TRUE ~ is_chip_gene
      )
    )
}

#' Annotate data with immunoglobulin gene information
#'
#' @param x A data frame/tibble or GWASFormatter object containing variant data.
#' @param ... Additional arguments passed to methods.
#'
#' @export
annotate_with_immunoglobulin = function(x, ...) {
  UseMethod("annotate_with_immunoglobulin")
}

#' @export
annotate_with_immunoglobulin.data.frame = function(x, ...) {
  result = x %>%
    dplyr::mutate(
      is_IGHV = ifelse(
        CHROM == "chr14" & POS >= 105586437 & POS <= 106879844,
        TRUE,
        FALSE
      ),
      is_IGLV = ifelse(
        CHROM == "chr22" & POS >= 22026076 & POS <= 22922913,
        TRUE,
        FALSE
      ),
    )

  return(result)
}

#' @export
annotate_with_immunoglobulin.tbl_df = function(x, ...) {
  annotate_with_immunoglobulin.data.frame(x, ...)
}

#' @export
annotate_with_immunoglobulin.GWASFormatter = function(x, ...) {
  con = db_connect()

  x$data = x$data %>%
    dplyr::mutate(
      is_IGHV = ifelse(
        CHROM == "chr14" & POS >= 105586437 & POS <= 106879844,
        TRUE,
        FALSE
      ),
      is_IGLV = ifelse(
        CHROM == "chr22" & POS >= 22026076 & POS <= 22922913,
        TRUE,
        FALSE
      )
    ) %>%
    dplyr::compute(temporary = FALSE, overwrite = TRUE, name = "summary_stats")

  return(x)
}

#' Identify independent GWAS loci using a greedy, distance-based algorithm
#'
#' Groups significant variants into independent loci. The lead variant (most
#' significant) for each locus is identified first; all other significant
#' variants on the same chromosome within `window_kb` kilobases are assigned
#' to that locus. No LD information is used.
#'
#' @param x A GWASFormatter object, data.frame, or tibble.
#' @param window_kb Window size in kilobases. Variants within this distance of
#'   the lead variant on the same chromosome are assigned to the same locus.
#'   Default is 500.
#' @param pvalue_threshold P-value significance threshold. Default is 5e-8.
#' @param ... Additional arguments (unused).
#' @return A tibble/data.frame of significant variants with added columns:
#'   \itemize{
#'     \item `locus_id` - Integer locus identifier
#'     \item `is_lead` - Logical, TRUE for the lead (most significant) variant per locus
#'   }
#' @export
identify_loci <- function(x, window_kb = 500, pvalue_threshold = 5e-8, ...) {
  UseMethod("identify_loci")
}

#' @describeIn identify_loci Method for data.frame/tibble objects
#' @export
identify_loci.data.frame <- function(
  x,
  window_kb = 500,
  pvalue_threshold = 5e-8,
  ...
) {
  window_bp <- window_kb * 1000L

  sig <- x %>%
    dplyr::filter(PVALUE <= pvalue_threshold) %>%
    dplyr::arrange(PVALUE)

  if (nrow(sig) == 0L) {
    cli::cli_inform("No variants below threshold {pvalue_threshold}.")
    return(sig %>% dplyr::mutate(locus_id = integer(0), is_lead = logical(0)))
  }

  # Greedy locus assignment
  locus_id <- rep(NA_integer_, nrow(sig))
  is_lead <- rep(FALSE, nrow(sig))
  next_locus <- 1L

  for (i in seq_len(nrow(sig))) {
    if (!is.na(locus_id[i])) {
      next
    }
    # Assign new locus
    locus_id[i] <- next_locus
    is_lead[i] <- TRUE
    # Propagate to nearby variants on the same chromosome
    same_chr <- sig$CHROM == sig$CHROM[i]
    within <- abs(sig$POS - sig$POS[i]) <= window_bp
    locus_id[same_chr & within & is.na(locus_id)] <- next_locus
    next_locus <- next_locus + 1L
  }

  sig %>%
    dplyr::mutate(locus_id = locus_id, is_lead = is_lead)
}

#' @describeIn identify_loci Alias for data.frame method
#' @export
identify_loci.tbl_df <- function(
  x,
  window_kb = 500,
  pvalue_threshold = 5e-8,
  ...
) {
  identify_loci.data.frame(x, window_kb, pvalue_threshold, ...)
}

#' @describeIn identify_loci Method for GWASFormatter objects
#' @export
identify_loci.GWASFormatter <- function(
  x,
  window_kb = 500,
  pvalue_threshold = 5e-8,
  ...
) {
  sig <- x$data %>%
    dplyr::filter(PVALUE <= pvalue_threshold) %>%
    dplyr::collect()

  identify_loci.data.frame(sig, window_kb, pvalue_threshold, ...)
}

#' Annotate top GWAS hits with Open Targets Locus-to-Gene (L2G) predictions
#'
#' Queries the Open Targets Platform to retrieve L2G scores and adds
#' `l2g_gene_id`, `l2g_gene_name`, and `l2g_score` columns.
#'
#' Two methods are available:
#' * `"api"` (default) — queries the GraphQL API once per variant.
#'   Fast for small sets (<100 variants), but slow for larger ones (~30 min
#'   for 1000 variants).
#' * `"bulk"` — uses DuckDB to query the Open Targets Platform parquet files
#'   (FTP). Downloads ~700 MB of L2G predictions on first use, then caches
#'   them permanently. Subsequent queries for any number of variants finish
#'   in seconds.
#'
#' The typical workflow for large sets:
#' ```r
#' hits <- select_top_hits(gwas_obj)
#' hits <- find_nearest_gene(hits)
#' hits <- annotate_with_l2g(hits, method = "bulk")
#' ```
#'
#' @param x A data.frame or tibble of GWAS hits.
#' @param id_col Name of the column containing variant IDs. Default `"ID"`.
#' @param method `"api"` (per-variant GraphQL, default) or `"bulk"` (parquet).
#' @param release Open Targets Platform release string, e.g. `"25.09"`.
#'   `NULL` (default) auto-detects the latest release.
#'   Only used when `method = "bulk"`.
#' @param cache_dir Directory for cached parquet files. Defaults to the
#'   platform user-data directory for gwasplot.
#'   Only used when `method = "bulk"`.
#' @param ask If `TRUE` (default when running interactively), prompt the user
#'   before downloading data. Set `FALSE` to download without prompting (e.g.
#'   in scripts). Only used when `method = "bulk"` and data is not cached.
#' @param ... Additional arguments (unused).
#' @return The input data.frame with three new columns appended:
#'   \describe{
#'     \item{l2g_gene_id}{Ensembl gene ID of the top L2G gene}
#'     \item{l2g_gene_name}{Approved gene symbol}
#'     \item{l2g_score}{L2G score 0–1 (higher = more likely causal)}
#'   }
#'   Variants with no L2G predictions receive `NA` in these columns.
#' @export
annotate_with_l2g <- function(x, id_col = "ID", ...) {
  UseMethod("annotate_with_l2g")
}

#' @export
annotate_with_l2g.data.frame <- function(
  x,
  id_col = "ID",
  method = c("api", "bulk"),
  release = NULL,
  cache_dir = tools::R_user_dir("gwasplot", "data"),
  ask = interactive(),
  ...
) {
  method <- match.arg(method)
  if (method == "api") {
    ids <- x[[id_col]]
    cli::cli_alert_info(
      "Querying Open Targets L2G for {length(ids)} variant(s)..."
    )

    results <- vector("list", length(ids))
    cli::cli_progress_bar("Querying L2G", total = length(ids))
    for (i in seq_along(ids)) {
      results[[i]] <- tryCatch(
        genesForVariant(ids[[i]]),
        error = function(e) {
          tibble::tibble(
            ID = ids[[i]],
            l2g_gene_id = NA_character_,
            l2g_gene_name = NA_character_,
            l2g_score = NA_real_
          )
        }
      )
      cli::cli_progress_update()
    }
    cli::cli_progress_done()

    result_tbl <- dplyr::bind_rows(results)
    dplyr::left_join(x, result_tbl, by = stats::setNames("ID", id_col))
  } else {
    .annotate_l2g_bulk(x, id_col, release, cache_dir, ask)
  }
}

#' @export
annotate_with_l2g.tbl_df <- function(
  x,
  id_col = "ID",
  method = c("api", "bulk"),
  release = NULL,
  cache_dir = tools::R_user_dir("gwasplot", "data"),
  ask = interactive(),
  ...
) {
  annotate_with_l2g.data.frame(x, id_col, method, release, cache_dir, ask, ...)
}

# ── Open Targets Platform bulk helpers (internal) ────────────────────────────

# Detect latest release from the Platform GraphQL meta endpoint.
.ot_release_latest <- function() {
  tryCatch(
    {
      resp <- httr::POST(
        "https://api.platform.opentargets.org/api/v4/graphql",
        httr::content_type_json(),
        body = jsonlite::toJSON(
          list(query = "{ meta { dataVersion { year month } } }"),
          auto_unbox = TRUE
        ),
        httr::timeout(15)
      )
      httr::stop_for_status(resp)
      v <- jsonlite::fromJSON(
        httr::content(resp, "text", encoding = "UTF-8"),
        flatten = TRUE
      )$data$meta$dataVersion
      # year and month are returned as strings (e.g. "25", "12")
      sprintf("%s.%02d", v$year, as.integer(v$month))
    },
    error = function(e) {
      cli::cli_abort(c(
        "!" = "Could not detect the latest Open Targets Platform release.",
        "i" = "Specify manually, e.g. {.code release = \"25.09\"}",
        "x" = conditionMessage(e)
      ))
    }
  )
}

# Enumerate all snappy.parquet URLs in an FTP directory listing.
.ot_ftp_list_parquet <- function(table, release) {
  base_url <- sprintf(
    "https://ftp.ebi.ac.uk/pub/databases/opentargets/platform/%s/output/%s/",
    release,
    table
  )
  tryCatch(
    {
      html <- httr::content(
        httr::GET(base_url, httr::timeout(15)),
        "text",
        encoding = "UTF-8"
      )
      files <- regmatches(
        html,
        gregexpr('href="([^"]+\\.snappy\\.parquet)"', html, perl = TRUE)
      )[[1]]
      names <- sub('href="([^"]+)"', "\\1", files)
      paste0(base_url, names)
    },
    error = function(e) {
      cli::cli_abort(
        "Could not enumerate parquet files for {.val {table}} ({release}): {conditionMessage(e)}"
      )
    }
  )
}

# Build remote parquet glob URL for a given table and release.
# Returns a character vector of explicit URLs (DuckDB httpfs requires this for HTTPS).
.ot_remote_glob <- function(table, release) {
  .ot_ftp_list_parquet(table, release)
}

# Probe FTP to resolve the credibleSet table directory name.
.ot_cs_table_name <- function(release) {
  for (tbl in c("credible_set", "credibleSet")) {
    url <- sprintf(
      "https://ftp.ebi.ac.uk/pub/databases/opentargets/platform/%s/output/%s/",
      release,
      tbl
    )
    resp <- tryCatch(httr::HEAD(url, httr::timeout(10)), error = function(e) {
      NULL
    })
    if (
      !is.null(resp) && httr::status_code(resp) %in% c(200L, 301L, 302L, 303L)
    ) {
      return(tbl)
    }
  }
  cli::cli_abort(
    "Cannot find the credibleSet table on the OT FTP for release {.val {release}}."
  )
}

# Parse an Apache FTP directory listing to sum file sizes and return MB.
.ot_dir_size_mb <- function(table, release) {
  url <- sprintf(
    "https://ftp.ebi.ac.uk/pub/databases/opentargets/platform/%s/output/%s/",
    release,
    table
  )
  tryCatch(
    {
      resp <- httr::GET(url, httr::timeout(15))
      if (httr::http_error(resp)) {
        return(NA_real_)
      }
      html <- httr::content(resp, "text", encoding = "UTF-8")
      # Apache FTP listings show human-readable sizes: "2.4M", "700M", "1.2G", "500K"
      m <- regmatches(
        html,
        gregexpr("([0-9]+\\.?[0-9]*)([KMG])", html, perl = TRUE)
      )[[1]]
      if (length(m) == 0L) {
        return(NA_real_)
      }
      nums <- as.numeric(sub("([0-9.]+)[KMG]", "\\1", m))
      units <- sub("[0-9.]+([KMG])", "\\1", m)
      mb_each <- ifelse(
        units == "K",
        nums / 1024,
        ifelse(units == "M", nums, ifelse(units == "G", nums * 1024, 0))
      )
      round(sum(mb_each, na.rm = TRUE), 0)
    },
    error = function(e) NA_real_
  )
}

# Download parquet files one-by-one via httr (EBI FTP does not support Range
# requests needed by DuckDB httpfs, so we download raw files sequentially).
# Resumes interrupted downloads by skipping already-complete files.
# Writes a '.complete' sentinel on success so callers can distinguish a full
# cache from a partial one.
# If select_sql / where_sql are set, each file is filtered at download time.
.ot_httpfs_copy <- function(
  remote_urls,
  local_dir,
  select_sql = "*",
  where_sql = NULL
) {
  dir.create(local_dir, recursive = TRUE, showWarnings = FALSE)

  filter_on_download <- !identical(select_sql, "*") || !is.null(where_sql)
  sentinel <- file.path(local_dir, ".complete")

  n_todo <- sum(!file.exists(file.path(local_dir, basename(remote_urls))))
  if (n_todo == 0L) {
    writeLines(as.character(Sys.time()), sentinel)
    return(invisible(local_dir))
  }

  cli::cli_progress_bar(
    "Downloading {n_todo} of {length(remote_urls)} file(s)",
    total = n_todo
  )
  for (url in remote_urls) {
    fname <- basename(url)
    dest <- file.path(local_dir, fname)

    if (file.exists(dest)) {
      next
    } # already downloaded

    if (filter_on_download) {
      # Download to a temp file, then filter into the cache dir
      tmp <- tempfile(fileext = ".parquet")
      on.exit(unlink(tmp), add = TRUE)
      resp <- httr::GET(
        url,
        httr::write_disk(tmp, overwrite = TRUE),
        httr::timeout(600)
      )
      httr::stop_for_status(resp)
      con <- DBI::dbConnect(duckdb::duckdb())
      where_clause <- if (!is.null(where_sql)) paste("WHERE", where_sql) else ""
      DBI::dbExecute(
        con,
        sprintf(
          "COPY (SELECT %s FROM read_parquet('%s') %s)
         TO '%s' (FORMAT PARQUET, COMPRESSION SNAPPY)",
          select_sql,
          gsub("\\\\", "/", tmp),
          where_clause,
          dest
        )
      )
      DBI::dbDisconnect(con, shutdown = TRUE)
    } else {
      resp <- httr::GET(
        url,
        httr::write_disk(dest, overwrite = TRUE),
        httr::timeout(600)
      )
      httr::stop_for_status(resp)
    }
    cli::cli_progress_update()
  }
  cli::cli_progress_done()
  writeLines(as.character(Sys.time()), sentinel)
  invisible(local_dir)
}

# Make a DuckDB-compatible parquet expression:
#   - A character vector of HTTPS URLs → DuckDB array literal ['url1', 'url2', ...]
#   - A single local directory path → 'dir/*.parquet'
.parquet_glob <- function(path) {
  if (length(path) > 1 || grepl("^https://", path[1])) {
    # Escape single quotes and wrap each URL
    escaped <- gsub("'", "\\'", path, fixed = TRUE)
    paste0("[", paste0("'", escaped, "'", collapse = ", "), "]")
  } else {
    paste0("'", gsub("\\\\", "/", path), "/*.parquet'")
  }
}

# Core DuckDB query: join credibleSet index + L2G predictions, return the
# top-scoring gene per input variant.
.ot_l2g_query <- function(variant_ids, l2g_path, cs_path) {
  norm_ids <- gsub(":", "_", sub("^chr", "", variant_ids, ignore.case = TRUE))
  use_httpfs <- any(grepl("^https://", c(unlist(l2g_path), unlist(cs_path))))

  con <- DBI::dbConnect(duckdb::duckdb())
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE)
  if (use_httpfs) {
    DBI::dbExecute(con, "INSTALL httpfs; LOAD httpfs;")
  }

  DBI::dbWriteTable(
    con,
    "input_ids",
    data.frame(
      input_id = variant_ids,
      norm_id = norm_ids,
      stringsAsFactors = FALSE
    ),
    overwrite = TRUE,
    temporary = TRUE
  )

  result <- DBI::dbGetQuery(
    con,
    sprintf(
      "WITH cs AS (
       SELECT studyLocusId, variantId
       FROM read_parquet(%s)
       WHERE variantId IN (SELECT norm_id FROM input_ids)
         AND studyType = 'gwas'
     ),
     l2g_agg AS (
       SELECT
         cs.variantId AS norm_id,
         l.geneId     AS l2g_gene_id,
         MAX(l.score) AS l2g_score
       FROM read_parquet(%s) l
       JOIN cs ON l.studyLocusId = cs.studyLocusId
       GROUP BY cs.variantId, l.geneId
     ),
     ranked AS (
       SELECT *,
         ROW_NUMBER() OVER (PARTITION BY norm_id ORDER BY l2g_score DESC) AS rn
       FROM l2g_agg
     )
     SELECT i.input_id AS ID, r.l2g_gene_id, r.l2g_score
     FROM ranked r
     JOIN input_ids i ON r.norm_id = i.norm_id
     WHERE r.rn = 1",
      .parquet_glob(cs_path),
      .parquet_glob(l2g_path)
    )
  ) %>%
    tibble::as_tibble()

  # Attach gene symbols from the built-in reference dataset.
  gene_map <- tibble::tibble(
    l2g_gene_id = human_genes$gene_id,
    l2g_gene_name = human_genes$gene_name
  ) %>%
    dplyr::distinct(l2g_gene_id, .keep_all = TRUE)

  result %>%
    dplyr::left_join(gene_map, by = "l2g_gene_id") %>%
    dplyr::select(ID, l2g_gene_id, l2g_gene_name, l2g_score)
}

# Orchestrates release detection, cache check, user prompt, download, and query.
.annotate_l2g_bulk <- function(x, id_col, release, cache_dir, ask) {
  # 1. Resolve release.
  if (is.null(release)) {
    cli::cli_progress_step("Detecting latest Open Targets Platform release")
    release <- .ot_release_latest()
  }
  cli::cli_alert_info("Open Targets Platform release: {.strong {release}}")

  # 2. Cache paths.
  l2g_cache <- file.path(cache_dir, sprintf("ot_%s_l2g", release))
  cs_cache <- file.path(cache_dir, sprintf("ot_%s_cs_gwas", release))

  is_complete <- function(d) file.exists(file.path(d, ".complete"))
  l2g_cached <- is_complete(l2g_cache)
  cs_cached <- is_complete(cs_cache)

  if (l2g_cached && cs_cached) {
    cli::cli_alert_success("Using cached data from {.path {cache_dir}}")
    l2g_path <- l2g_cache
    cs_path <- cs_cache
  } else {
    # 3. Find the credibleSet table name on FTP.
    cli::cli_progress_step("Locating credibleSet table on Open Targets FTP")
    cs_tbl <- .ot_cs_table_name(release)

    # 4. Estimate download sizes.
    l2g_mb <- .ot_dir_size_mb("l2g_prediction", release)
    cs_mb <- .ot_dir_size_mb(cs_tbl, release)
    l2g_size_str <- if (is.na(l2g_mb)) "~700 MB" else sprintf("%d MB", l2g_mb)
    # credibleSet index: only 3 columns via parquet column pruning (~2% of full table)
    cs_size_str <- if (is.na(cs_mb)) {
      "size unknown (3-col extract)"
    } else {
      sprintf("~%d MB (3-col extract)", max(1L, round(cs_mb / 50L)))
    }

    # 5. Show what will be downloaded and prompt.
    to_download <- c(
      if (!l2g_cached) {
        sprintf("L2G predictions      (%s)  →  %s", l2g_size_str, l2g_cache)
      },
      if (!cs_cached) {
        sprintf("credibleSet index    (%s)  →  %s", cs_size_str, cs_cache)
      }
    )
    cli::cli_inform(c(
      "!" = "Open Targets Platform {.strong {release}} data not found in cache.",
      "i" = "The following will be downloaded:",
      setNames(to_download, rep(" ", length(to_download))),
      "i" = "Files are cached permanently and reused for all future queries."
    ))

    proceed <- if (ask) {
      ans <- readline("Proceed with download? [Y/n]: ")
      nchar(ans) == 0L || tolower(trimws(ans)) == "y"
    } else {
      TRUE
    }

    if (!proceed) {
      cli::cli_abort(c(
        "x" = "Download cancelled.",
        "i" = "Re-run with {.code ask = FALSE} to download automatically, or use",
        "i" = "{.code method = \"api\"} for per-variant GraphQL queries instead."
      ))
    } else {
      if (!l2g_cached) {
        cli::cli_progress_step("Downloading L2G predictions ({l2g_size_str})")
        .ot_httpfs_copy(.ot_remote_glob("l2g_prediction", release), l2g_cache)
      }
      if (!cs_cached) {
        cli::cli_progress_step(
          "Downloading credibleSet GWAS index ({cs_size_str})"
        )
        .ot_httpfs_copy(
          .ot_remote_glob(cs_tbl, release),
          cs_cache,
          select_sql = "studyLocusId, variantId, studyType",
          where_sql = "studyType = 'gwas'"
        )
      }
      l2g_path <- l2g_cache
      cs_path <- cs_cache
    }
  }

  # 6. Query.
  ids <- x[[id_col]]
  cli::cli_progress_step("Querying L2G for {length(ids)} variant(s)")
  result_tbl <- .ot_l2g_query(ids, l2g_path, cs_path)

  n_found <- sum(!is.na(result_tbl$l2g_gene_id))
  cli::cli_alert_success(
    "L2G annotations found for {n_found} of {length(ids)} variant(s)"
  )

  dplyr::left_join(x, result_tbl, by = stats::setNames("ID", id_col))
}

# ── End bulk helpers ──────────────────────────────────────────────────────────

#' Query Open Targets Platform API for Locus-to-Gene (L2G) predictions
#'
#' @param variant_id A string representing the variant ID in format "CHR_POS_REF_ALT"
#'   (e.g., "19_44908822_C_T") or an rsID (e.g., "rs123456").
#' @return A tibble with columns `ID`, `l2g_gene_id`, `l2g_gene_name`, `l2g_score`.
#' @examples
#' \dontrun{
#'   result <- query_ot_api_v2g("19_44908822_C_T")
#'   result <- query_ot_api_v2g("rs123456")
#' }
#' @export
query_ot_api_v2g = function(variant_id = "19_44908822_C_T") {
  genesForVariant(variant_id)
}

genesForVariant <- function(variant_id) {
  na_row <- tibble::tibble(
    ID = variant_id,
    l2g_gene_id = NA_character_,
    l2g_gene_name = NA_character_,
    l2g_score = NA_real_
  )

  # Helper: POST a GraphQL query to the OT Platform API and return parsed data.
  ot_post <- function(query, variables = list()) {
    resp <- httr::POST(
      "https://api.platform.opentargets.org/api/v4/graphql",
      httr::content_type_json(),
      body = jsonlite::toJSON(
        list(query = query, variables = variables),
        auto_unbox = TRUE
      ),
      httr::timeout(30)
    )
    httr::stop_for_status(resp)
    jsonlite::fromJSON(
      httr::content(resp, "text", encoding = "UTF-8"),
      flatten = TRUE
    )$data
  }

  tryCatch(
    {
      # Normalise variant ID: strip leading chr prefix, handle : separator
      normalised_id <- sub("^chr", "", variant_id)
      normalised_id <- gsub(":", "_", normalised_id)

      if (grepl("^rs\\d+$", normalised_id, ignore.case = TRUE)) {
        # Convert rsID to variant ID via Platform API search
        query_searchid <- "query SearchQuery($queryString: String!, $index: Int!, $entityNames: [String!]!) {
        search(queryString: $queryString, entityNames: $entityNames, page: {index: $index, size: 10}) {
          hits {
            object {
              ... on Variant { id __typename }
            }
          }
        }
      }"
        id_result <- ot_post(
          query_searchid,
          list(
            queryString = normalised_id,
            index = 0L,
            entityNames = list("Variant")
          )
        )

        # flatten=TRUE converts hits$object.{id,__typename} into top-level columns
        hits <- id_result$search$hits
        variant_rows <- hits[
          !is.na(hits[["object.__typename"]]) &
            hits[["object.__typename"]] == "Variant",
        ]

        if (nrow(variant_rows) == 0 || is.na(variant_rows[["object.id"]][1])) {
          return(na_row)
        }
        input_variant_id <- variant_rows[["object.id"]][1]
      } else if (grepl("^\\d+_\\d+_[a-zA-Z]+_[a-zA-Z]+$", normalised_id)) {
        input_variant_id <- normalised_id
      } else {
        cli::cli_warn("Unrecognised variant ID format: {variant_id}")
        return(na_row)
      }

      # Simplified query: locus subfield removed (was fetched but never used)
      query <- "query GWASCredibleSetsQuery($variantId: String!, $size: Int!, $index: Int!) {
      variant(variantId: $variantId) {
        credibleSets(studyTypes: [gwas], page: { size: $size, index: $index }) {
          rows {
            l2GPredictions(page: {index: 0, size: 100}) {
              rows {
                score
                target { id approvedSymbol }
              }
            }
          }
        }
      }
    }"

      result <- ot_post(
        query,
        list(variantId = input_variant_id, size = 25L, index = 0L)
      )

      if (is.null(result$variant)) {
        return(na_row)
      }

      cs_rows <- result$variant$credibleSets$rows
      if (
        is.null(cs_rows) ||
          length(cs_rows) == 0 ||
          (is.data.frame(cs_rows) && nrow(cs_rows) == 0)
      ) {
        return(na_row)
      }

      # Aggregate L2G predictions across all credible sets
      all_preds <- list()
      for (i in seq_len(nrow(cs_rows))) {
        l2g_rows <- cs_rows[["l2GPredictions.rows"]][[i]]
        if (
          !is.null(l2g_rows) && is.data.frame(l2g_rows) && nrow(l2g_rows) > 0
        ) {
          all_preds[[length(all_preds) + 1]] <- tibble::tibble(
            gene_id = l2g_rows[["target.id"]],
            gene_name = l2g_rows[["target.approvedSymbol"]],
            l2g_score = l2g_rows[["score"]]
          )
        }
      }

      if (length(all_preds) == 0) {
        return(na_row)
      }

      # Dedup: keep max score per gene (not first-encountered)
      combined <- dplyr::bind_rows(all_preds)
      best <- dplyr::group_by(combined, gene_id)
      best <- dplyr::slice_max(best, l2g_score, n = 1, with_ties = FALSE)
      best <- dplyr::ungroup(best)
      best <- dplyr::arrange(best, dplyr::desc(l2g_score))

      if (nrow(best) == 0) {
        return(na_row)
      }

      tibble::tibble(
        ID = variant_id,
        l2g_gene_id = best$gene_id[[1]],
        l2g_gene_name = best$gene_name[[1]],
        l2g_score = best$l2g_score[[1]]
      )
    },
    error = function(e) {
      cli::cli_warn(
        "genesForVariant failed for {variant_id}: {conditionMessage(e)}"
      )
      na_row
    }
  )
}

#' Query Open Targets Platform API for variant information
#'
#' @param variant_id A string representing the variant ID (e.g., "19_44908822_C_T").
#' @return A data frame containing variant information.
#' @export
query_ot_api_variants <- function(variant_id = "19_44908822_C_T") {
  na_row <- function(id_value = variant_id) {
    data.frame(
      id = id_value %||% NA_character_,
      rsId = NA_character_,
      nearestCodingGene.symbol = NA_character_,
      nearestCodingGene.id = NA_character_,
      nearestCodingGeneDistance = NA_real_,
      nearestGeneDistance = NA_real_,
      mostSevereConsequence = NA_character_,
      caddPhred = NA_real_,
      gnomadNFE = NA_real_,
      gnomadAFR = NA_real_,
      gnomadAMR = NA_real_,
      gnomadEAS = NA_real_,
      stringsAsFactors = FALSE
    )
  }

  # Check if the variant ID argument is empty or null
  if (is.null(variant_id) || variant_id == "") {
    warning("Please provide a value for the variant ID argument.")
    return(na_row(NA_character_))
  }

  # Helper: POST a GraphQL query to the OT Platform API and return parsed data.
  ot_post <- function(query, variables = list()) {
    resp <- httr::POST(
      "https://api.platform.opentargets.org/api/v4/graphql",
      httr::content_type_json(),
      body = jsonlite::toJSON(
        list(query = query, variables = variables),
        auto_unbox = TRUE
      ),
      httr::timeout(30)
    )
    httr::stop_for_status(resp)
    jsonlite::fromJSON(
      httr::content(resp, "text", encoding = "UTF-8"),
      flatten = TRUE
    )$data
  }

  # Try-catch block for handling connection timeout
  tryCatch(
    {
      cli::cli_progress_step(
        "Connecting to the Open Targets Platform GraphQL API...",
        spinner = TRUE
      )

      # Normalize variant IDs like "chr1:123_A_T" -> "1_123_A_T"
      normalized_variant_id <- gsub(
        ":",
        "_",
        sub("^chr", "", variant_id, ignore.case = TRUE)
      )

      # Check variant id format
      if (
        grepl(pattern = "^rs\\d+$", normalized_variant_id, ignore.case = TRUE)
      ) {
        # Convert rs id to variant id using new Platform API search
        query_searchid <- "query SearchQuery($queryString: String!, $index: Int!, $entityNames: [String!]!) {
        search(
          queryString: $queryString
          entityNames: $entityNames
          page: {index: $index, size: 10}
        ) {
          total
          hits {
            id
            object {
              ... on Variant {
                id
                variantDescription
                referenceAllele
                alternateAllele
                rsIds
                __typename
              }
            }
          }
        }
      }"

        id_result <- ot_post(
          query_searchid,
          list(
            queryString = normalized_variant_id,
            index = 0L,
            entityNames = list("Variant")
          )
        )

        variant_hits <- id_result$search$hits
        if (!is.data.frame(variant_hits) || nrow(variant_hits) == 0) {
          stop(paste("No variant found for rsID:", normalized_variant_id))
        }

        # flatten=TRUE usually produces object.__typename / object.id columns.
        if (
          "object.__typename" %in%
            names(variant_hits) &&
            "object.id" %in% names(variant_hits)
        ) {
          variant_objects <- variant_hits[
            !is.na(variant_hits[["object.__typename"]]) &
              variant_hits[["object.__typename"]] == "Variant",
            ,
            drop = FALSE
          ]
          if (
            nrow(variant_objects) == 0 ||
              is.na(variant_objects[["object.id"]][1])
          ) {
            stop(paste("No variant found for rsID:", normalized_variant_id))
          }
          input_variant_id <- variant_objects[["object.id"]][1]
        } else if ("id" %in% names(variant_hits)) {
          # Fallback for alternate flattening behavior.
          candidate <- variant_hits[["id"]][1]
          if (
            is.na(candidate) ||
              !grepl("^\\d+_\\d+_[A-Za-z]+_[A-Za-z]+$", candidate)
          ) {
            stop(paste("No variant found for rsID:", normalized_variant_id))
          }
          input_variant_id <- candidate
        } else {
          stop(paste("No variant found for rsID:", normalized_variant_id))
        }
      } else if (
        grepl(
          pattern = "^\\d+_\\d+_[a-zA-Z]+_[a-zA-Z]+$",
          normalized_variant_id
        )
      ) {
        input_variant_id <- normalized_variant_id
      } else {
        stop("\nPlease provide a variant ID")
      }

      # Check if the input_variant_id is null or empty
      if (is.null(input_variant_id) || input_variant_id == "") {
        stop(
          "There is no variant ID defined for this rsID by Open Targets Platform"
        )
      }

      # Define the query for new Platform API
      query <- "query VariantInfoQuery($variantId: String!) {
      variant(variantId: $variantId) {
        id
        rsIds
        chromosome
        position
        referenceAllele
        alternateAllele
        variantDescription
        mostSevereConsequence {
          id
          label
        }
        alleleFrequencies {
          populationName
          alleleFrequency
        }
        transcriptConsequences {
          target {
            id
            approvedSymbol
          }
          distanceFromTss
          impact
        }
      }
    }"

      cli::cli_progress_step("Downloading variant data...", spinner = TRUE)
      var_info <- ot_post(query, list(variantId = input_variant_id))

      # Process and return data frame compatible with old format
      if (!is.null(var_info$variant)) {
        variant_data <- var_info$variant

        # Find nearest gene from transcript consequences
        nearest_gene_symbol <- NA_character_
        nearest_gene_id <- NA_character_
        nearest_gene_distance <- NA_real_

        tc <- variant_data$transcriptConsequences
        if (!is.null(tc) && length(tc) > 0) {
          # Support both flattened data.frame and nested list responses.
          tc_df <- if (is.data.frame(tc)) {
            tc
          } else if (is.list(tc)) {
            tryCatch(dplyr::bind_rows(tc), error = function(e) data.frame())
          } else {
            data.frame()
          }

          if (is.data.frame(tc_df) && nrow(tc_df) > 0) {
            if (!"distanceFromTss" %in% names(tc_df)) {
              tc_df$distanceFromTss <- NA_real_
            }
            if (!"target.approvedSymbol" %in% names(tc_df)) {
              tc_df$target.approvedSymbol <- NA_character_
            }
            if (!"target.id" %in% names(tc_df)) {
              tc_df$target.id <- NA_character_
            }

            dist <- suppressWarnings(abs(as.numeric(tc_df$distanceFromTss)))
            dist[is.na(dist)] <- Inf
            min_dist_idx <- which.min(dist)

            if (length(min_dist_idx) == 1 && is.finite(dist[min_dist_idx])) {
              nearest_gene_symbol <- tc_df$target.approvedSymbol[[min_dist_idx]]
              nearest_gene_id <- tc_df$target.id[[min_dist_idx]]
              nearest_gene_distance <- as.numeric(tc_df$distanceFromTss[[
                min_dist_idx
              ]])
            }
          }
        }

        # Extract allele frequencies
        gnomad_nfe <- NA_real_
        gnomad_afr <- NA_real_
        gnomad_amr <- NA_real_
        gnomad_eas <- NA_real_

        af <- variant_data$alleleFrequencies
        if (!is.null(af) && length(af) > 0) {
          af_df <- if (is.data.frame(af)) {
            af
          } else if (is.list(af)) {
            tryCatch(dplyr::bind_rows(af), error = function(e) data.frame())
          } else {
            data.frame()
          }

          if (is.data.frame(af_df) && nrow(af_df) > 0) {
            if (!"populationName" %in% names(af_df)) {
              af_df$populationName <- NA_character_
            }
            if (!"alleleFrequency" %in% names(af_df)) {
              af_df$alleleFrequency <- NA_real_
            }

            if (any(grepl("NFE", af_df$populationName, ignore.case = TRUE))) {
              gnomad_nfe <- af_df$alleleFrequency[grepl(
                "NFE",
                af_df$populationName,
                ignore.case = TRUE
              )][1]
            }
            if (any(grepl("AFR", af_df$populationName, ignore.case = TRUE))) {
              gnomad_afr <- af_df$alleleFrequency[grepl(
                "AFR",
                af_df$populationName,
                ignore.case = TRUE
              )][1]
            }
            if (any(grepl("AMR", af_df$populationName, ignore.case = TRUE))) {
              gnomad_amr <- af_df$alleleFrequency[grepl(
                "AMR",
                af_df$populationName,
                ignore.case = TRUE
              )][1]
            }
            if (any(grepl("EAS", af_df$populationName, ignore.case = TRUE))) {
              gnomad_eas <- af_df$alleleFrequency[grepl(
                "EAS",
                af_df$populationName,
                ignore.case = TRUE
              )][1]
            }
          }
        }

        severe_label <- NA_character_
        msc <- variant_data$mostSevereConsequence
        if (is.data.frame(msc) && "label" %in% names(msc) && nrow(msc) > 0) {
          severe_label <- msc$label[[1]]
        } else if (is.list(msc) && !is.null(msc$label)) {
          severe_label <- msc$label
        }

        # Create data frame compatible with old format
        df_var_info <- data.frame(
          id = variant_data$id %||% NA_character_,
          rsId = paste(variant_data$rsIds, collapse = ",") %||% NA_character_,
          nearestCodingGene.symbol = nearest_gene_symbol,
          nearestCodingGene.id = nearest_gene_id,
          nearestCodingGeneDistance = nearest_gene_distance,
          nearestGeneDistance = nearest_gene_distance, # Same as coding gene distance for compatibility
          mostSevereConsequence = severe_label,
          caddPhred = NA_real_, # CADD scores not available in Platform API
          gnomadNFE = gnomad_nfe,
          gnomadAFR = gnomad_afr,
          gnomadAMR = gnomad_amr,
          gnomadEAS = gnomad_eas,
          stringsAsFactors = FALSE
        )

        return(df_var_info)
      } else {
        warning("No variant information found for ID: ", input_variant_id)
        return(na_row(input_variant_id))
      }
    },
    error = function(e) {
      # Handling connection timeout
      if (grepl("Timeout was reached", e$message)) {
        warning(
          "Connection timeout reached while connecting to the Open Targets Platform GraphQL API."
        )
      } else {
        warning("Error querying Open Targets Platform API: ", e$message)
      }
      return(na_row(variant_id))
    }
  )
}

#' Annotate variants with functional consequences using Ensembl VEP API
#'
#' This function queries the Ensembl Variant Effect Predictor (VEP) REST API
#' to annotate a list of variants with their functional consequences. It processes
#' variants in batches to respect API rate limits and returns detailed consequence
#' information including transcript effects, protein changes, and regulatory impacts.
#'
#' @param variants A character vector containing variant information in the format
#'   "chr_pos_ref_alt" (e.g., "chr21_26960070_G_A" or "21_26960070_G_A").
#' @param batch_size Integer specifying the number of variants to process per API call.
#'   Default is 200 (max 1000). Smaller batches are more reliable for large datasets.
#' @param flag_pick Logical. If TRUE, only report the transcript with the PICK flag.
#'   Default is TRUE to get the single best transcript per variant.
#' @param include_hgvs Logical. If TRUE, include HGVS nomenclature in results. Default is TRUE.
#' @param include_domains Logical. If TRUE, include protein domain information. Default is FALSE.
#' @param include_regulatory Logical. If TRUE, include regulatory feature consequences. Default is FALSE.
#' @param sleep_time Numeric. Time in seconds to wait between API calls to respect rate limits.
#'   Default is 1 second. Increase if you encounter rate limiting.
#' @param verbose Logical. If TRUE, print progress messages. Default is TRUE.
#'
#' @return A tibble containing variant consequences with the following columns:
#'   \itemize{
#'     \item ID - Original variant string
#'     \item CHROM - Chromosome (extracted from variant)
#'     \item POS - Position (extracted from variant)
#'     \item REF - Reference allele (extracted from variant)
#'     \item ALT - Alternate allele (extracted from variant)
#'     \item most_severe_consequence - Most severe consequence for the variant
#'     \item gene_id - Ensembl gene ID
#'     \item gene_symbol - Gene symbol
#'     \item transcript_id - Ensembl transcript ID
#'     \item consequence_terms - All consequence terms (comma-separated)
#'     \item impact - Impact level (HIGH, MODERATE, LOW, MODIFIER)
#'     \item protein_position - Position in protein sequence
#'     \item amino_acids - Reference and alternate amino acids
#'     \item codons - Reference and alternate codons
#'     \item existing_variation - Known variant IDs (e.g., rsIDs)
#'     \item hgvsc - HGVS coding sequence nomenclature (if include_hgvs=TRUE)
#'     \item hgvsp - HGVS protein nomenclature (if include_hgvs=TRUE)
#'     \item domains - Protein domains affected (if include_domains=TRUE)
#'   }
#'
#' @details
#' The function handles variant input in the format "chr_pos_ref_alt":
#' \itemize{
#'   \item Input format: "chr21_26960070_G_A" or "21_26960070_G_A"
#'   \item Automatic rate limiting to respect Ensembl's 15 requests/second limit
#'   \item Batch processing for efficient handling of large variant lists
#'   \item Error handling and retry logic for failed requests
#' }
#'
#' @examples
#' \dontrun{
#' # Example 1: Single variant
#' variants <- "chr21_26960070_G_A"
#' consequences <- annotate_variants_ensembl(variants)
#'
#' # Example 2: Multiple variants
#' variants <- c(
#'   "chr21_26960070_G_A",
#'   "21_26965148_G_A",
#'   "chrX_155066068_C_T"
#' )
#' consequences <- annotate_variants_ensembl(variants)
#'
#' # Example 3: With additional options
#' consequences <- annotate_variants_ensembl(
#'   variants,
#'   flag_pick = TRUE,
#'   include_domains = TRUE,
#'   batch_size = 100
#' )
#' }
#'
#' @note
#' \itemize{
#'   \item Respects Ensembl API rate limits (15 requests/second)
#'   \item Large variant lists are automatically batched
#'   \item Requires internet connection to Ensembl REST API
#'   \item For very large datasets (>10,000 variants), consider using local VEP installation
#'   \item Only supports human variants (homo_sapiens)
#' }
#'
#' @references
#' McLaren et al. (2016). The Ensembl Variant Effect Predictor.
#' Genome Biology 17, 122. doi:10.1186/s13059-016-0974-4
#'
#' @seealso
#' \url{https://rest.ensembl.org/documentation/info/vep_region_post}
#' \url{https://github.com/Ensembl/ensembl-vep}
#'
#' @export
annotate_variants_ensembl <- function(
  variants,
  batch_size = 200,
  flag_pick = TRUE, # Use flag_pick instead of canonical_only
  include_hgvs = TRUE,
  include_domains = FALSE,
  include_regulatory = FALSE,
  sleep_time = 1,
  verbose = FALSE
) {
  # Load required packages
  # Validate inputs
  if (length(variants) == 0) {
    stop("No variants provided")
  }

  if (batch_size > 1000) {
    warning(
      "Batch size > 1000 may cause API errors. Consider using smaller batches."
    )
    batch_size <- 1000
  }

  # Convert input to standard format
  variant_strings <- format_variants_for_vep(variants)

  if (verbose) {
    message(sprintf(
      "Annotating %d variants using Ensembl VEP API",
      length(variant_strings)
    ))
    message(sprintf("Processing in batches of %d variants", batch_size))
  }

  # Split variants into batches
  n_variants <- length(variant_strings)
  n_batches <- ceiling(n_variants / batch_size)

  all_results <- list()

  # Process each batch
  for (i in 1:n_batches) {
    if (verbose) {
      message(sprintf("Processing batch %d of %d...", i, n_batches))
    }

    # Calculate batch indices
    start_idx <- (i - 1) * batch_size + 1
    end_idx <- min(i * batch_size, n_variants)
    batch_variants <- variant_strings[start_idx:end_idx]

    # Query VEP API for this batch
    batch_result <- query_vep_batch(
      batch_variants,
      flag_pick = flag_pick,
      include_hgvs = include_hgvs,
      include_domains = include_domains,
      include_regulatory = include_regulatory,
      verbose = verbose
    )

    if (!is.null(batch_result)) {
      all_results[[i]] <- batch_result
    }

    # Rate limiting: sleep between requests (except for the last batch)
    if (i < n_batches && sleep_time > 0) {
      Sys.sleep(sleep_time)
    }
  }

  # Combine all results
  if (length(all_results) == 0) {
    warning("No results returned from VEP API")
    return(tibble::tibble())
  }

  final_results <- dplyr::bind_rows(all_results)

  # Parse original variants to extract CHROM, POS, REF, ALT
  # Convert VCF format back to original variant format for joining
  final_results$original_variant <- sapply(
    final_results$input,
    function(vcf_string) {
      # VCF format: "chr pos id ref alt qual filter info"
      parts <- strsplit(vcf_string, " ")[[1]]
      if (length(parts) >= 5) {
        chr <- parts[1]
        pos <- parts[2]
        ref <- parts[4]
        alt <- parts[5]

        # Add chr prefix if not present
        if (!grepl("^chr", chr)) {
          chr <- paste0("chr", chr)
        }

        return(paste(chr, pos, ref, alt, sep = "_"))
      }
      return(NA_character_)
    }
  )

  # Parse original variants for CHROM, POS, REF, ALT
  parsed_variants <- parse_variants_for_output(variants)

  # Join with parsed variant information
  final_results <- final_results %>%
    dplyr::left_join(
      parsed_variants,
      by = c("original_variant" = "original")
    ) %>%
    dplyr::select(
      ID = original_variant,
      CHROM,
      POS,
      REF,
      ALT,
      most_severe_consequence,
      gene_id,
      gene_symbol,
      transcript_id,
      consequence_terms,
      impact,
      protein_position,
      amino_acids,
      codons,
      existing_variation,
      hgvsc,
      hgvsp,
      domains
    ) %>%
    tibble::as_tibble()

  if (verbose) {
    message(sprintf("Successfully annotated %d variants", nrow(final_results)))
  }

  return(final_results)
}

#' Format variants for VEP API input
#'
#' @param variants Input variants in format "chr_pos_ref_alt"
#' @return Character vector of VCF-formatted variant strings
#' @keywords internal
format_variants_for_vep <- function(variants) {
  if (!is.character(variants)) {
    stop("Variants must be a character vector in format 'chr_pos_ref_alt'")
  }

  # Parse variants in format "chr_pos_ref_alt"
  vcf_strings <- sapply(
    variants,
    function(variant) {
      parts <- strsplit(variant, "_")[[1]]

      if (length(parts) != 4) {
        stop(sprintf(
          "Invalid variant format: %s. Expected format: 'chr_pos_ref_alt'",
          variant
        ))
      }

      chr <- parts[1]
      pos <- parts[2]
      ref <- parts[3]
      alt <- parts[4]

      # Remove 'chr' prefix if present
      chr <- gsub("^chr", "", chr)

      # Create VCF-like string: chr pos id ref alt qual filter info
      paste(chr, pos, ".", ref, alt, ".", ".", ".")
    },
    USE.NAMES = FALSE
  )

  return(vcf_strings)
}

#' Parse variants to extract CHROM, POS, REF, ALT for output
#'
#' @param variants Input variants in format "chr_pos_ref_alt"
#' @return data.frame with original variant and parsed components
#' @keywords internal
parse_variants_for_output <- function(variants) {
  if (!is.character(variants)) {
    stop("Variants must be a character vector in format 'chr_pos_ref_alt'")
  }

  # Parse variants in format "chr_pos_ref_alt"
  parsed <- data.frame(
    original = variants,
    CHROM = character(length(variants)),
    POS = integer(length(variants)),
    REF = character(length(variants)),
    ALT = character(length(variants)),
    stringsAsFactors = FALSE
  )

  for (i in seq_along(variants)) {
    variant <- variants[i]
    parts <- strsplit(variant, "_")[[1]]

    if (length(parts) != 4) {
      stop(sprintf(
        "Invalid variant format: %s. Expected format: 'chr_pos_ref_alt'",
        variant
      ))
    }

    parsed$CHROM[i] <- parts[1]
    parsed$POS[i] <- as.integer(parts[2])
    parsed$REF[i] <- parts[3]
    parsed$ALT[i] <- parts[4]
  }

  return(parsed)
}

#' Query VEP API for a batch of variants
#'
#' @param variant_batch Character vector of VCF-formatted variants
#' @param flag_pick Include flag_pick to mark the selected transcript
#' @param include_hgvs Include HGVS nomenclature
#' @param include_domains Include protein domains
#' @param include_regulatory Include regulatory consequences
#' @return data.frame with VEP results
#' @keywords internal
query_vep_batch <- function(
  variant_batch,
  flag_pick = FALSE,
  include_hgvs = TRUE,
  include_domains = FALSE,
  include_regulatory = FALSE,
  verbose = TRUE
) {
  # Construct API URL (fixed to homo_sapiens)
  base_url <- "https://rest.ensembl.org"
  endpoint <- "/vep/homo_sapiens/region"
  url <- paste0(base_url, endpoint)

  # Construct request body
  request_body <- list(variants = variant_batch)

  # Add VEP options as URL parameters
  query_params <- list()

  if (flag_pick) {
    query_params[["flag_pick"]] <- "1"
  }

  if (include_hgvs) {
    query_params[["hgvs"]] <- "1"
  }

  if (include_domains) {
    query_params[["domains"]] <- "1"
  }

  if (include_regulatory) {
    query_params[["regulatory"]] <- "1"
  }

  # Set headers
  headers <- c(
    "Content-Type" = "application/json",
    "Accept" = "application/json"
  )

  # Make API request with error handling
  tryCatch(
    {
      response <- httr::POST(
        url = url,
        query = query_params,
        body = jsonlite::toJSON(request_body, auto_unbox = FALSE),
        httr::add_headers(.headers = headers),
        httr::timeout(120) # 2 minute timeout
      )

      # Check for rate limiting
      if (httr::status_code(response) == 429) {
        retry_after <- httr::headers(response)[["retry-after"]]
        if (!is.null(retry_after)) {
          wait_time <- as.numeric(retry_after)
          message(sprintf("Rate limited. Waiting %d seconds...", wait_time))
          Sys.sleep(wait_time)

          # Retry the request
          response <- httr::POST(
            url = url,
            query = query_params,
            body = jsonlite::toJSON(request_body, auto_unbox = FALSE),
            httr::add_headers(.headers = headers),
            httr::timeout(120)
          )
        }
      }

      # Check response status
      if (httr::status_code(response) != 200) {
        # Get response content for debugging
        error_content <- httr::content(response, "text", encoding = "UTF-8")
        warning(sprintf(
          "API request failed with status %d. URL: %s\nRequest body: %s\nError response: %s",
          httr::status_code(response),
          url,
          jsonlite::toJSON(request_body, auto_unbox = FALSE),
          substr(error_content, 1, 500)
        )) # Truncate long error messages
        return(NULL)
      }

      # Parse response
      content <- httr::content(response, "text", encoding = "UTF-8")

      # Debug: show raw content structure
      if (verbose) {
        message("Raw API response (first 500 chars):")
        message(substr(content, 1, 500))
      }

      vep_results <- jsonlite::fromJSON(content, flatten = TRUE) # Try flattened first

      # Debug: print structure to understand the response
      if (verbose) {
        message("API response structure:")
        message(paste("Length:", length(vep_results)))
        message(paste("Class:", class(vep_results)))
        if (is.data.frame(vep_results)) {
          message(paste(
            "Column names:",
            paste(names(vep_results), collapse = ", ")
          ))
          message(paste("Number of rows:", nrow(vep_results)))
          # Show first few column names for debugging
          if (ncol(vep_results) > 10) {
            message(paste(
              "First 10 columns:",
              paste(names(vep_results)[1:10], collapse = ", ")
            ))
          }
        } else if (length(vep_results) > 0) {
          message(paste(
            "First element names:",
            paste(names(vep_results[[1]]), collapse = ", ")
          ))
        }
      }

      # Process and format results
      formatted_results <- format_vep_results(vep_results, flag_pick)

      return(formatted_results)
    },
    error = function(e) {
      warning(sprintf("Error querying VEP API: %s", e$message))
      return(NULL)
    }
  )
}

#' Format VEP API results into a tidy data.frame
#'
#' @param vep_results Raw results from VEP API (can be data.frame or list)
#' @param flag_pick Logical. If TRUE, filter for transcripts with PICK flag
#' @return Formatted data.frame
#' @keywords internal
format_vep_results <- function(vep_results, flag_pick = FALSE) {
  if (is.null(vep_results) || length(vep_results) == 0) {
    return(data.frame())
  }

  # Handle case where VEP returns a data.frame directly (flattened results)
  if (is.data.frame(vep_results)) {
    # Check if transcript_consequences is a column with nested data
    if ("transcript_consequences" %in% names(vep_results)) {
      # transcript_consequences contains nested list data
      result_list <- list()

      for (i in 1:nrow(vep_results)) {
        row_data <- vep_results[i, ]
        input_variant <- row_data$input %||% NA_character_
        most_severe <- row_data$most_severe_consequence %||% NA_character_

        # Extract colocated variants
        existing_variation <- ""
        if (
          !is.null(row_data$colocated_variants) &&
            length(row_data$colocated_variants[[1]]) > 0
        ) {
          colocated_data <- row_data$colocated_variants[[1]]
          if (
            is.data.frame(colocated_data) && "id" %in% names(colocated_data)
          ) {
            existing_ids <- colocated_data$id[!is.na(colocated_data$id)]
            existing_variation <- paste(existing_ids, collapse = ",")
          }
        }

        # Extract transcript consequences
        if (
          !is.null(row_data$transcript_consequences) &&
            length(row_data$transcript_consequences[[1]]) > 0
        ) {
          transcript_data <- row_data$transcript_consequences[[1]]

          # If transcript_data is a data.frame with multiple rows (multiple transcripts)
          if (is.data.frame(transcript_data)) {
            # Filter for transcripts with PICK flag if requested
            if (flag_pick && "PICK" %in% names(transcript_data)) {
              pick_transcripts <- transcript_data[
                transcript_data$PICK == 1,
                ,
                drop = FALSE
              ]
              if (nrow(pick_transcripts) > 0) {
                transcript_data <- pick_transcripts
              } else {
                # If no PICK transcripts found, take the first one
                transcript_data <- transcript_data[1, , drop = FALSE]
              }
            } else if (flag_pick) {
              # If flag_pick is TRUE but no PICK column, just take the first transcript
              transcript_data <- transcript_data[1, , drop = FALSE]
            }

            for (j in 1:nrow(transcript_data)) {
              tc <- transcript_data[j, ]

              # Extract consequence terms
              consequence_terms <- ""
              if (
                !is.null(tc$consequence_terms) &&
                  length(tc$consequence_terms[[1]]) > 0
              ) {
                consequence_terms <- paste(
                  tc$consequence_terms[[1]],
                  collapse = ","
                )
              }

              result_list[[length(result_list) + 1]] <- data.frame(
                input = input_variant,
                most_severe_consequence = most_severe,
                gene_id = tc$gene_id %||% NA_character_,
                gene_symbol = tc$gene_symbol %||% NA_character_,
                transcript_id = tc$transcript_id %||% NA_character_,
                consequence_terms = consequence_terms,
                impact = tc$impact %||% NA_character_,
                protein_position = as.character(tc$protein_start %||% NA),
                amino_acids = tc$amino_acids %||% NA_character_,
                codons = tc$codons %||% NA_character_,
                existing_variation = existing_variation,
                hgvsc = tc$hgvsc %||% NA_character_,
                hgvsp = tc$hgvsp %||% NA_character_,
                domains = "", # Will be processed separately if needed
                stringsAsFactors = FALSE
              )
            }
          }
        } else {
          # No transcript consequences
          result_list[[length(result_list) + 1]] <- data.frame(
            input = input_variant,
            most_severe_consequence = most_severe,
            gene_id = NA_character_,
            gene_symbol = NA_character_,
            transcript_id = NA_character_,
            consequence_terms = NA_character_,
            impact = NA_character_,
            protein_position = NA_character_,
            amino_acids = NA_character_,
            codons = NA_character_,
            existing_variation = existing_variation,
            hgvsc = NA_character_,
            hgvsp = NA_character_,
            domains = NA_character_,
            stringsAsFactors = FALSE
          )
        }
      }

      return(dplyr::bind_rows(result_list))
    }

    # Fallback: Extract relevant columns and create standardized output
    result_df <- data.frame(
      input = vep_results$input %||% NA_character_,
      most_severe_consequence = vep_results$most_severe_consequence %||%
        NA_character_,
      gene_id = vep_results$transcript_consequences.gene_id %||%
        vep_results$gene_id %||%
        NA_character_,
      gene_symbol = vep_results$transcript_consequences.gene_symbol %||%
        vep_results$gene_symbol %||%
        NA_character_,
      transcript_id = vep_results$transcript_consequences.transcript_id %||%
        vep_results$transcript_id %||%
        NA_character_,
      consequence_terms = vep_results$transcript_consequences.consequence_terms %||%
        vep_results$consequence_terms %||%
        NA_character_,
      impact = vep_results$transcript_consequences.impact %||%
        vep_results$impact %||%
        NA_character_,
      protein_position = as.character(
        vep_results$transcript_consequences.protein_start %||%
          vep_results$protein_start %||%
          NA
      ),
      amino_acids = vep_results$transcript_consequences.amino_acids %||%
        vep_results$amino_acids %||%
        NA_character_,
      codons = vep_results$transcript_consequences.codons %||%
        vep_results$codons %||%
        NA_character_,
      existing_variation = "", # Will be filled below
      hgvsc = vep_results$transcript_consequences.hgvsc %||%
        vep_results$hgvsc %||%
        NA_character_,
      hgvsp = vep_results$transcript_consequences.hgvsp %||%
        vep_results$hgvsp %||%
        NA_character_,
      domains = "", # Will be filled below
      stringsAsFactors = FALSE
    )

    # Handle existing variation (colocated variants)
    colocated_cols <- grep(
      "colocated_variants",
      names(vep_results),
      value = TRUE
    )
    if (length(colocated_cols) > 0) {
      # Try to extract IDs from colocated variants columns
      id_cols <- grep(
        "colocated_variants.*\\.id",
        names(vep_results),
        value = TRUE
      )
      if (length(id_cols) > 0) {
        existing_ids <- unlist(vep_results[id_cols])
        result_df$existing_variation <- paste(
          existing_ids[!is.na(existing_ids)],
          collapse = ","
        )
      }
    }

    return(result_df)
  }

  # Handle case where VEP returns a list (original structure)
  # Initialize result list
  result_list <- list()

  # Process each variant result
  for (i in seq_along(vep_results)) {
    variant_result <- vep_results[[i]]

    # Extract basic variant information
    input_variant <- variant_result$input %||% NA_character_
    most_severe <- variant_result$most_severe_consequence %||% NA_character_

    # Extract colocated variants (existing variation)
    existing_variation <- ""
    if (
      !is.null(variant_result$colocated_variants) &&
        length(variant_result$colocated_variants) > 0
    ) {
      # Handle case where colocated_variants is a data.frame or list
      if (is.data.frame(variant_result$colocated_variants)) {
        existing_ids <- variant_result$colocated_variants$id
      } else if (is.list(variant_result$colocated_variants)) {
        existing_ids <- sapply(variant_result$colocated_variants, function(x) {
          x$id %||% NA
        })
      } else {
        existing_ids <- character(0)
      }
      existing_variation <- paste(
        existing_ids[!is.na(existing_ids)],
        collapse = ","
      )
    }

    # Check if there are transcript consequences
    if (
      is.null(variant_result$transcript_consequences) ||
        length(variant_result$transcript_consequences) == 0
    ) {
      # No transcript consequences, create minimal record
      result_list[[i]] <- data.frame(
        input = input_variant,
        most_severe_consequence = most_severe,
        gene_id = NA_character_,
        gene_symbol = NA_character_,
        transcript_id = NA_character_,
        consequence_terms = NA_character_,
        impact = NA_character_,
        protein_position = NA_character_,
        amino_acids = NA_character_,
        codons = NA_character_,
        existing_variation = existing_variation,
        hgvsc = NA_character_,
        hgvsp = NA_character_,
        domains = NA_character_,
        stringsAsFactors = FALSE
      )
    } else {
      # Process transcript consequences
      transcript_data <- variant_result$transcript_consequences

      # Handle case where transcript_consequences is a list of lists
      if (is.list(transcript_data) && !is.data.frame(transcript_data)) {
        # Convert list of transcript consequences to data.frame
        transcript_rows <- list()
        for (j in seq_along(transcript_data)) {
          tc <- transcript_data[[j]]

          # Extract consequence terms
          consequence_terms <- ""
          if (!is.null(tc$consequence_terms)) {
            if (is.character(tc$consequence_terms)) {
              consequence_terms <- paste(tc$consequence_terms, collapse = ",")
            } else if (is.list(tc$consequence_terms)) {
              consequence_terms <- paste(
                unlist(tc$consequence_terms),
                collapse = ","
              )
            }
          }

          # Extract domains
          domains <- ""
          if (!is.null(tc$domains) && length(tc$domains) > 0) {
            if (is.list(tc$domains)) {
              domain_names <- sapply(tc$domains, function(d) d$db %||% "")
              domains <- paste(domain_names[domain_names != ""], collapse = ",")
            }
          }

          transcript_rows[[j]] <- data.frame(
            input = input_variant,
            most_severe_consequence = most_severe,
            gene_id = tc$gene_id %||% NA_character_,
            gene_symbol = tc$gene_symbol %||% NA_character_,
            transcript_id = tc$transcript_id %||% NA_character_,
            consequence_terms = consequence_terms,
            impact = tc$impact %||% NA_character_,
            protein_position = as.character(tc$protein_start %||% NA),
            amino_acids = tc$amino_acids %||% NA_character_,
            codons = tc$codons %||% NA_character_,
            existing_variation = existing_variation,
            hgvsc = tc$hgvsc %||% NA_character_,
            hgvsp = tc$hgvsp %||% NA_character_,
            domains = domains,
            stringsAsFactors = FALSE
          )
        }
        result_list[[i]] <- dplyr::bind_rows(transcript_rows)
      } else {
        # transcript_data is already a data.frame, process directly
        result_list[[i]] <- data.frame(
          input = input_variant,
          most_severe_consequence = most_severe,
          gene_id = transcript_data$gene_id %||% NA_character_,
          gene_symbol = transcript_data$gene_symbol %||% NA_character_,
          transcript_id = transcript_data$transcript_id %||% NA_character_,
          consequence_terms = sapply(
            transcript_data$consequence_terms %||% list(NA),
            function(x) paste(x, collapse = ",")
          ),
          impact = transcript_data$impact %||% NA_character_,
          protein_position = as.character(
            transcript_data$protein_start %||% NA
          ),
          amino_acids = transcript_data$amino_acids %||% NA_character_,
          codons = transcript_data$codons %||% NA_character_,
          existing_variation = existing_variation,
          hgvsc = transcript_data$hgvsc %||% NA_character_,
          hgvsp = transcript_data$hgvsp %||% NA_character_,
          domains = sapply(transcript_data$domains %||% list(NA), function(x) {
            if (is.list(x)) {
              paste(sapply(x, function(d) d$db %||% ""), collapse = ",")
            } else {
              NA_character_
            }
          }),
          stringsAsFactors = FALSE
        )
      }
    }
  }

  # Combine all results
  if (length(result_list) > 0) {
    final_result <- dplyr::bind_rows(result_list)
  } else {
    final_result <- data.frame()
  }

  return(final_result)
}
