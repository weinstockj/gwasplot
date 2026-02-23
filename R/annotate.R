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
find_nearest_gene = function(x, ...) {
  UseMethod("find_nearest_gene")
}

#' @describeIn find_nearest_gene Find the nearest gene for each variant in a gwas object
#' @param threshold The distance threshold to consider a gene as nearest. Default is 1e5.
#' @export
find_nearest_gene.GWASFormatter = function(x, threshold = 1e5, ...) {
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
  ...
) {
  find_nearest_gene.data.frame(x, threshold, chrom_col, pos_col, id_col, ...)
}

#' @export
find_nearest_gene.data.frame = function(
  x,
  threshold = 1e5,
  chrom_col = "CHROM",
  pos_col = "POS",
  id_col = "ID",
  ...
) {
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

#' Query Open Targets Platform API for Locus-to-Gene (L2G) predictions
#'
#' This function retrieves L2G predictions from credible sets containing the specified variant.
#' It replaces the deprecated V2G functionality from the old Open Targets Genetics API.
#'
#' @param variant_id A string representing the variant ID in format "CHR_POS_REF_ALT"
#'   (e.g., "19_44908822_C_T") or an rsID (e.g., "rs123456").
#' @param pageindex Pagination index (currently unused in new API structure).
#' @param pagesize Pagination size (currently unused in new API structure).
#' @return A list containing L2G predictions and associated data from credible sets.
#' @note This function has been migrated from the deprecated Open Targets Genetics API
#'   to the new Open Targets Platform API. The functionality has changed from direct
#'   variant-to-gene mapping to locus-to-gene predictions via credible sets.
#' @examples
#' \dontrun{
#'   # Query by variant ID
#'   result <- query_ot_api_v2g("19_44908822_C_T")
#'
#'   # Query by rsID
#'   result <- query_ot_api_v2g("rs123456")
#' }
#' @export
query_ot_api_v2g = function(
  variant_id = "19_44908822_C_T",
  pageindex = 0,
  pagesize = 20
) {
  # modified frOm: https://github.com/amirfeizi/otargen/blob/main/R/indexVariantsAndStudiesForTagVariant.R
  genesForVariant(variant_id)
}

genesForVariant <- function(variant_id) {
  ## Set up to query Open Targets Platform API (migrated from deprecated Genetics API)
  tryCatch(
    {
      cli::cli_progress_step(
        "Connecting to the Open Targets Platform GraphQL API...",
        spinner = TRUE
      )
      otg_cli <- ghql::GraphqlClient$new(
        url = "https://api.platform.opentargets.org/api/v4/graphql"
      )
      otg_qry <- ghql::Query$new()

      # Check variant id format
      if (grepl(pattern = "rs\\d+", variant_id)) {
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

        variables <- list(
          queryString = variant_id,
          index = 0L,
          entityNames = list("Variant")
        )

        otg_qry$query(name = "convertid", x = query_searchid)
        id_result <- jsonlite::fromJSON(
          otg_cli$exec(otg_qry$queries$convertid, variables),
          flatten = TRUE
        )$data

        # Filter for variant objects and get the first variant ID
        variant_hits <- id_result$search$hits
        variant_objects <- variant_hits[
          sapply(variant_hits$object, function(x) {
            !is.null(x$`__typename`) && x$`__typename` == "Variant"
          }),
        ]

        if (
          length(variant_objects) == 0 ||
            is.null(variant_objects$object[[1]]$id)
        ) {
          stop(paste("No variant found for rsID:", variant_id))
        }

        input_variant_id <- variant_objects$object[[1]]$id
      } else if (grepl(pattern = "\\d+_\\d+_[a-zA-Z]+_[a-zA-Z]+", variant_id)) {
        input_variant_id <- variant_id
      } else {
        stop("\nPlease provide a variant ID.")
      }

      # New Platform API query for L2G predictions via credible sets
      # This replaces the old genesForVariant query which is no longer available
      query <- "query GWASCredibleSetsQuery($variantId: String!, $size: Int!, $index: Int!) {
      variant(variantId: $variantId) {
        id
        referenceAllele
        alternateAllele
        chromosome
        position
        credibleSets(studyTypes: [gwas], page: { size: $size, index: $index }) {
          count
          rows {
            studyLocusId
            pValueMantissa
            pValueExponent
            beta
            finemappingMethod
            confidence
            study {
              traitFromSource
              id
              diseases {
                name
                id
              }
            }
            l2GPredictions(page: {index: 0, size: 100}) {
              rows {
                score
                target {
                  id
                  approvedSymbol
                }
                features {
                  name
                  value
                  shapValue
                }
              }
            }
            locus(variantIds: [$variantId]) {
              rows {
                posteriorProbability
                pValueExponent
                pValueMantissa
                beta
                is95CredibleSet
                is99CredibleSet
                variant {
                  id
                  rsIds
                }
              }
            }
          }
        }
      }
    }"

      ## Execute the query

      result_pkg <- list()

      variables <- list(
        variantId = input_variant_id,
        size = 500L,
        index = 0L
      )

      otg_qry$query(name = "l2g_query", x = query)
      cli::cli_progress_step(
        paste0("Downloading L2G data for ", variant_id, " ..."),
        spinner = TRUE
      )

      result <- jsonlite::fromJSON(
        otg_cli$exec(otg_qry$queries$l2g_query, variables),
        flatten = TRUE
      )$data

      # Extract L2G predictions from credible sets
      if (
        !is.null(result$variant$credibleSets$rows) &&
          length(result$variant$credibleSets$rows) > 0
      ) {
        # Aggregate L2G predictions across all credible sets
        all_l2g_predictions <- list()

        for (i in seq_along(result$variant$credibleSets$rows)) {
          credible_set <- result$variant$credibleSets$rows[[i]]
          if (
            !is.null(credible_set$l2GPredictions$rows) &&
              length(credible_set$l2GPredictions$rows) > 0
          ) {
            l2g_data <- credible_set$l2GPredictions$rows

            # Add credible set context
            for (j in seq_along(l2g_data)) {
              l2g_data[[j]]$studyLocusId <- credible_set$studyLocusId
              l2g_data[[
                j
              ]]$traitFromSource <- credible_set$study$traitFromSource
              l2g_data[[j]]$variant <- input_variant_id
            }

            all_l2g_predictions <- c(all_l2g_predictions, l2g_data)
          }
        }

        if (length(all_l2g_predictions) > 0) {
          # Convert to data frame format similar to old genesForVariant output
          result_df <- data.frame(
            gene.symbol = sapply(all_l2g_predictions, function(x) {
              x$target$approvedSymbol %||% NA_character_
            }),
            gene.id = sapply(all_l2g_predictions, function(x) {
              x$target$id %||% NA_character_
            }),
            variant = sapply(all_l2g_predictions, function(x) {
              x$variant %||% NA_character_
            }),
            overallScore = sapply(all_l2g_predictions, function(x) {
              x$score %||% NA_real_
            }),
            studyLocusId = sapply(all_l2g_predictions, function(x) {
              x$studyLocusId %||% NA_character_
            }),
            traitFromSource = sapply(all_l2g_predictions, function(x) {
              x$traitFromSource %||% NA_character_
            }),
            stringsAsFactors = FALSE
          )

          # Remove duplicates and sort by score
          result_df <- result_df[
            !duplicated(paste(result_df$gene.symbol, result_df$variant)),
          ]
          result_df <- result_df[
            order(result_df$overallScore, decreasing = TRUE),
          ]

          # Create simplified output similar to old format
          result_core <- result_df %>%
            dplyr::select(gene.symbol, variant, overallScore, gene.id) %>%
            dplyr::arrange(desc(overallScore))

          # Note: QTL data structure has changed significantly in Platform API
          # For backward compatibility, create empty QTL columns
          result_qtl <- result_df %>%
            dplyr::select(gene.symbol, variant) %>%
            dplyr::mutate(qtls = NA_character_)

          # Note: Intervals data structure has changed significantly in Platform API
          # For backward compatibility, create empty intervals columns
          result_intervals <- result_df %>%
            dplyr::select(gene.symbol, variant) %>%
            dplyr::mutate(intervals = NA_character_)

          # Note: Functional predictions and distances are not directly available in new Platform API
          # For backward compatibility, create empty columns
          result_functionalPredictions <- result_df %>%
            dplyr::select(gene.symbol, variant) %>%
            dplyr::mutate(functionalPredictions = NA_character_)

          result_distances <- result_df %>%
            dplyr::select(gene.symbol, variant) %>%
            dplyr::mutate(distances = NA_character_)

          result_pkg$core <- result_core
          result_pkg$qtls <- result_qtl
          result_pkg$intervals <- result_intervals
          result_pkg$functionalPredictions <- result_functionalPredictions
          result_pkg$distances <- result_distances

          return(result_pkg)
        } else {
          warning("No L2G predictions found for variant: ", variant_id)
          return(NULL)
        }
      } else {
        warning("No credible sets found for variant: ", variant_id)
        return(NULL)
      }
    },
    error = function(e) {
      if (grepl("timeout", e$message, ignore.case = TRUE)) {
        stop(
          "Connection timeout reached while connecting to the Open Targets Platform GraphQL API."
        )
      } else {
        stop("Error querying Open Targets Platform API: ", e$message)
      }
    }
  )
}

#' Query Open Targets Platform API for variant information
#'
#' @param variant_id A string representing the variant ID (e.g., "19_44908822_C_T").
#' @return A data frame containing variant information.
#' @export
query_ot_api_variants <- function(variant_id = "19_44908822_C_T") {
  # Check if the variant ID argument is empty or null
  if (is.null(variant_id) || variant_id == "") {
    message("Please provide a value for the variant ID argument.")
    return(NULL)
  }
  print(variant_id)

  # Try-catch block for handling connection timeout
  tryCatch(
    {
      # Set up to query Open Targets Platform API (migrated from deprecated Genetics API)
      cli::cli_progress_step(
        "Connecting to the Open Targets Platform GraphQL API...",
        spinner = TRUE
      )
      otg_cli <- ghql::GraphqlClient$new(
        url = "https://api.platform.opentargets.org/api/v4/graphql"
      )
      otg_qry <- ghql::Query$new()

      # Check variant id format
      if (grepl(pattern = "rs\\d+", variant_id)) {
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

        variables <- list(
          queryString = variant_id,
          index = 0L,
          entityNames = list("Variant")
        )

        otg_qry$query(name = "rsi2vid", x = query_searchid)
        id_result <- jsonlite::fromJSON(
          otg_cli$exec(otg_qry$queries$rsi2vid, variables),
          flatten = TRUE
        )$data

        # Filter for variant objects and get the first variant ID
        variant_hits <- id_result$search$hits
        variant_objects <- variant_hits[
          sapply(variant_hits$object, function(x) {
            !is.null(x$`__typename`) && x$`__typename` == "Variant"
          }),
        ]

        if (
          length(variant_objects) == 0 ||
            is.null(variant_objects$object[[1]]$id)
        ) {
          stop(paste("No variant found for rsID:", variant_id))
        }

        input_variant_id <- variant_objects$object[[1]]$id
      } else if (grepl(pattern = "\\d+_\\d+_[a-zA-Z]+_[a-zA-Z]+", variant_id)) {
        input_variant_id <- variant_id
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

      # Execute the query
      variables <- list(variantId = input_variant_id)
      otg_qry$query(name = "variantInfoquery", x = query)
      cli::cli_progress_step("Downloading variant data...", spinner = TRUE)
      var_info <- jsonlite::fromJSON(
        otg_cli$exec(otg_qry$queries$variantInfoquery, variables),
        flatten = TRUE
      )$data

      # Process and return data frame compatible with old format
      if (!is.null(var_info$variant)) {
        variant_data <- var_info$variant

        # Find nearest gene from transcript consequences
        nearest_gene_symbol <- NA_character_
        nearest_gene_id <- NA_character_
        nearest_gene_distance <- NA_real_

        if (
          !is.null(variant_data$transcriptConsequences) &&
            length(variant_data$transcriptConsequences) > 0
        ) {
          # Find the transcript consequence with minimum distance
          min_dist_idx <- which.min(sapply(
            variant_data$transcriptConsequences,
            function(x) abs(x$distanceFromTss %||% Inf)
          ))
          if (length(min_dist_idx) > 0) {
            nearest_tc <- variant_data$transcriptConsequences[[min_dist_idx]]
            nearest_gene_symbol <- nearest_tc$target$approvedSymbol
            nearest_gene_id <- nearest_tc$target$id
            nearest_gene_distance <- nearest_tc$distanceFromTss
          }
        }

        # Extract allele frequencies
        gnomad_nfe <- NA_real_
        gnomad_afr <- NA_real_
        gnomad_amr <- NA_real_
        gnomad_eas <- NA_real_

        if (!is.null(variant_data$alleleFrequencies)) {
          for (af in variant_data$alleleFrequencies) {
            if (grepl("NFE", af$populationName, ignore.case = TRUE)) {
              gnomad_nfe <- af$alleleFrequency
            }
            if (grepl("AFR", af$populationName, ignore.case = TRUE)) {
              gnomad_afr <- af$alleleFrequency
            }
            if (grepl("AMR", af$populationName, ignore.case = TRUE)) {
              gnomad_amr <- af$alleleFrequency
            }
            if (grepl("EAS", af$populationName, ignore.case = TRUE)) {
              gnomad_eas <- af$alleleFrequency
            }
          }
        }

        # Create data frame compatible with old format
        df_var_info <- data.frame(
          id = variant_data$id %||% NA_character_,
          rsId = paste(variant_data$rsIds, collapse = ",") %||% NA_character_,
          nearestCodingGene.symbol = nearest_gene_symbol,
          nearestCodingGene.id = nearest_gene_id,
          nearestCodingGeneDistance = nearest_gene_distance,
          nearestGeneDistance = nearest_gene_distance, # Same as coding gene distance for compatibility
          mostSevereConsequence = variant_data$mostSevereConsequence$label %||%
            NA_character_,
          caddPhred = NA_real_, # CADD scores not available in Platform API
          gnomadNFE = gnomad_nfe,
          gnomadAFR = gnomad_afr,
          gnomadAMR = gnomad_amr,
          gnomadEAS = gnomad_eas,
          stringsAsFactors = FALSE
        )

        return(df_var_info)
      } else {
        stop("No variant information found for ID: ", input_variant_id)
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
      return(NULL)
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
