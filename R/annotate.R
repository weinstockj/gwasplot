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
  sql_index = glue("
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
  WHERE gene_biotype = 'protein_coding'
  ")

  intervals = dplyr::tbl(con, dplyr::sql(sql_index)) %>%
    dplyr::compute(temporary = FALSE, overwrite = TRUE, name = "gene_intervals")

  DBI::dbExecute(
    con,
    "CREATE INDEX IF NOT EXISTS chrom_start_end ON gene_intervals (chrom, expanded_start, expanded_end)"
  )
  DBI::dbExecute(con, "CREATE INDEX IF NOT EXISTS chrom_pos ON summary_stats (chrom, POS)")
  
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
find_nearest_gene.tbl_df = function(x, threshold = 1e5, chrom_col = "CHROM", pos_col = "POS", id_col = "ID", ...) {
  # Start timing
  start_time <- Sys.time()
  cli::cli_alert_info("Starting gene annotation for tibble...")
  
  # Rename columns for consistency
  renamed_df <- x %>%
    dplyr::rename(
      CHROM = !!rlang::sym(chrom_col),
      POS = !!rlang::sym(pos_col),
      ID = !!rlang::sym(id_col)
    )

  other_cols = setdiff(names(renamed_df), c("CHROM", "POS", "ID", "gene_id" , "gene_name", "distance"))
  
  # Load and filter human genes data
  genes_df <- human_genes %>%
    dplyr::filter(gene_biotype == "protein_coding") %>%
    dplyr::mutate(
      expanded_start = start - threshold,
      expanded_end = end + threshold
    )
  
  nearest_genes <- renamed_df %>%
    dplyr::nest_by(CHROM, POS, ID, .keep = TRUE) %>%
    dplyr::reframe({
      nearby_genes <- genes_df %>%
        dplyr::filter(
          chrom == CHROM,
          expanded_start <= POS,
          expanded_end >= POS
        ) %>%
        dplyr::mutate(
          distance = dplyr::case_when(
            POS >= start & POS <= end ~ 0L,
            POS < start             ~ as.integer(start - POS),
            TRUE                     ~ as.integer(POS - end)
          )
        ) %>%
        dplyr::arrange(distance)
      
      if (nrow(nearby_genes) == 0L) {
        tibble::tibble(
          gene_id   = NA_character_,
          gene_name = NA_character_,
          distance  = NA_integer_
        )
      } else {
        nearby_genes %>%
          dplyr::slice(1) %>%
          dplyr::select(gene_id, gene_name, distance)
      }
    })
  
  renamed_df = renamed_df %>%
    dplyr::left_join(
      nearest_genes,
      by = c("CHROM", "POS", "ID")
    ) %>%
    dplyr::mutate(
      distance = dplyr::if_else(is.na(distance), NA_integer_, as.integer(distance))
    )
  # End timing and report
  end_time <- Sys.time()
  elapsed <- round(difftime(end_time, start_time, units = "secs"), 2)
  cli::cli_alert_success("Gene annotation completed in {elapsed} seconds")
  
  return(renamed_df)
}

#' @export
find_nearest_gene.data.frame = function(x, threshold = 1e5, chrom_col = "CHROM", pos_col = "POS", id_col = "ID", ...) {
  find_nearest_gene.tbl_df(x, threshold, chrom_col, pos_col, id_col, ...)
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
annotate_with_centromere.data.frame = function(x, chrom_col = "CHROM", pos_col = "POS", ...) {
  # Start timing
  start_time <- Sys.time()
  cli::cli_progress_step("Starting centromere annotation for data frame...")
  
  # Rename columns if needed
  if (chrom_col != "chrom" || pos_col != "POS") {
    top_hits <- x %>%
      dplyr::rename(
        chrom = !!rlang::sym(chrom_col),
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
        dplyr::select(chrom, start, end, name, in_centromere),
      by = dplyr::join_by(chrom, between(POS, start, end))
    ) %>%
    dplyr::mutate(
      in_centromere = dplyr::if_else(is.na(in_centromere), FALSE, in_centromere)
    )
  
  # Restore original column names if they were renamed
  if (chrom_col != "chrom") {
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

#' Annotate top hits with CHIP gene information
#' 
#' @param top_hits A data frame containing the top hits.
#' @return A data frame with the CHIP gene information.
#' @export
annotate_with_chip_genes = function(top_hits) {

  if (!"gene_name" %in% names(top_hits)) {
    cli::cli_abort("top_hits must contain a gene_name column; did you run find_nearest_gene?")
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

  sql = glue("
  SELECT
    t.*,
    c.is_chip_gene
  FROM top_hits t
  LEFT JOIN chip_genes c ON (t.gene_name = c.gene_name)
  ")

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
      is_IGHV = ifelse(chrom == "chr14" & POS >= 105586437 & POS <= 106879844, TRUE, FALSE),
      is_IGLV = ifelse(chrom == "chr22" & POS >= 22026076 & POS <= 22922913, TRUE, FALSE),
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
      is_IGHV = ifelse(chrom == "chr14" & POS >= 105586437 & POS <= 106879844, TRUE, FALSE),
      is_IGLV = ifelse(chrom == "chr22" & POS >= 22026076 & POS <= 22922913, TRUE, FALSE)
    ) %>%
    dplyr::compute(temporary = FALSE, overwrite = TRUE, name = "summary_stats")


  return(x)
}

query_ot_api_v2g = function(variant_id = "19_44908822_C_T", pageindex = 0, pagesize = 20) {


  # modified frOm: https://github.com/amirfeizi/otargen/blob/main/R/indexVariantsAndStudiesForTagVariant.R
  genesForVariant(variant_id)
}

genesForVariant <- function(variant_id) {
  ## Set up to query Open Targets Genetics API
  tryCatch({
    cli::cli_progress_step("Connecting to the Open Targets Genetics GraphQL API...", spinner = TRUE)
    otg_cli <- ghql::GraphqlClient$new(url = "https://api.genetics.opentargets.org/graphql")
    otg_qry <- ghql::Query$new()

    # Check variant id format
    if (grepl(pattern = "rs\\d+", variant_id)) {
      # Convert rs id to variant id
      query_searchid <- "query ConvertRSIDtoVID($queryString:String!) {
        search(queryString:$queryString){
          totalVariants
          variants{
            id
          }
        }
      }"

      variables <- list(queryString = variant_id)

      otg_qry$query(name = "convertid", x = query_searchid)
      id_result <- jsonlite::fromJSON(otg_cli$exec(otg_qry$queries$convertid, variables), flatten = TRUE)$data
      input_variant_id <- id_result$search$variants$id
    } else if (grepl(pattern = "\\d+_\\d+_[a-zA-Z]+_[a-zA-Z]+", variant_id)) {
      input_variant_id <- variant_id
    } else {
      stop("\nPlease provide a variant ID.")
    }

    query <- "query v2gquery($variantId: String!){
  genesForVariant(variantId: $variantId) {
    gene{
      id
      symbol
    }
    variant
    overallScore
    qtls{
      typeId
      aggregatedScore
      tissues{
        tissue{
          id
          name
        }
        quantile
        beta
        pval
      }
    }
    intervals{
      typeId
      sourceId
      aggregatedScore
      tissues{
        tissue{
          id
          name
        }
        quantile
        score
      }
    }
    functionalPredictions{
      typeId
      sourceId
      aggregatedScore
      tissues{
        tissue{
          id
          name
        }
        maxEffectLabel
        maxEffectScore
      }
    }
    distances{
      typeId
      sourceId
      aggregatedScore
      tissues{
        tissue{
          id
          name
        }
        distance
        score
        quantile
      }
    }
  }
}"

    ## Execute the query

    result_pkg <- list()

    variables <- list(variantId = input_variant_id)

    otg_qry$query(name = "v2g_query", x = query)
    cli::cli_progress_step(paste0("Downloading data for ", variant_id, " ..."), spinner = TRUE)

    result <- jsonlite::fromJSON(otg_cli$exec(otg_qry$queries$v2g_query, variables), flatten = TRUE)$data
    result_df <- as.data.frame(result$genesForVariant)

    if (nrow(result_df) != 0) {
      # parsing the nested JSON output in tidy data table format
      result_core <- result_df %>%
        dplyr::select(gene.symbol, variant, overallScore, gene.id) %>%
        dplyr::arrange(desc(overallScore))

      # qtl
      qtls_is_empty <-  all(sapply(result_df$qtls, function(x) length(x) == 0))

      if (qtls_is_empty) {
        result_qtl <- result_df %>%
          dplyr::select(gene.symbol, variant, qtls)
        result_qtl <- result_qtl %>%
          mutate(qtls = ifelse(sapply(qtls, length) == 0, NA_character_, toString(qtls)))
      } else {
        result_qtl <- result_df %>%
          dplyr::select(gene.symbol, variant, qtls) %>%
          tidyr::unnest(qtls, names_sep = '.', keep_empty = TRUE) %>%
          dplyr::rename("typeId" = "qtls.typeId",
                        "aggregatedScore" = "qtls.aggregatedScore")

        if ("qtls.tissues" %in% colnames(result_qtl)) {
      result_qtl <- result_qtl %>%
            tidyr::unnest(qtls.tissues, names_sep = '_', keep_empty = TRUE ) %>%
            dplyr::rename("tissues_id" = "qtls.tissues_tissue.id",
                          "tissues_name" = "qtls.tissues_tissue.name")
          base::colnames(result_qtl) <- stringr::str_replace_all(colnames(result_qtl), "qtls.", "")
        }
      }

      # intervals
      ints_is_empty <-  all(sapply(result_df$intervals, function(x) length(x) == 0))

      if (ints_is_empty) {
        result_intervals <- result_df %>%
          dplyr::select(gene.symbol, variant, intervals)
        result_intervals <- result_intervals %>%
          mutate(intervals = ifelse(sapply(intervals, length) == 0, NA_character_, toString(intervals)))
      } else {
        result_intervals <- result_df %>%
          dplyr::select(gene.symbol, variant, intervals) %>%
          tidyr::unnest(intervals, names_sep = '.', keep_empty = TRUE) %>%
          dplyr::rename("typeId" = "intervals.typeId",
                        "aggregatedScore" = "intervals.aggregatedScore")

        if ("intervals.tissues" %in% colnames(result_intervals)) {
        result_intervals <- result_intervals %>%
            tidyr::unnest(intervals.tissues, names_sep = '_', keep_empty = TRUE) %>%
            dplyr::rename("tissues_id" = "intervals.tissues_tissue.id",
                          "tissues_name" = "intervals.tissues_tissue.name")
          base::colnames(result_intervals) <- stringr::str_replace_all(colnames(result_intervals), "intervals.", "")
        }
      }

      # distances
      dists_is_empty <- all(sapply(result_df$distances, function(x) length(x) == 0))

      if (dists_is_empty) {
        result_distances <- result_df %>%
          dplyr::select(gene.symbol, variant, distances)
        result_distances <- result_distances %>%
          mutate(distances = ifelse(sapply(distances, length) == 0, NA_character_, toString(distances)))
      } else {
        result_distances <- result_df %>%
          dplyr::select(gene.symbol, variant, distances) %>%
          tidyr::unnest(distances, names_sep = '.', keep_empty = TRUE) %>%
          dplyr::rename("typeId" = "distances.typeId",
                        "aggregatedScore" = "distances.aggregatedScore")

        if ("distances.tissues" %in% colnames(result_distances)) {
          result_distances <- result_distances %>%
            tidyr::unnest(distances.tissues, names_sep = '_', keep_empty = TRUE) %>%
            dplyr::rename("tissues_id" = "distances.tissues_tissue.id",
                          "tissues_name" = "distances.tissues_tissue.name")
          base::colnames(result_distances) <- stringr::str_replace_all(colnames(result_distances), "distances.", "")
        }
      }

      # result_functionalPredictions
      funcPreds_is_empty <- all(sapply(result_df$functionalPredictions, function(x) length(x) == 0))

      if (funcPreds_is_empty) {
        result_functionalPredictions <- result_df %>%
          dplyr::select(gene.symbol, variant, functionalPredictions)
        result_functionalPredictions <- result_functionalPredictions %>%
          mutate(functionalPredictions = ifelse(sapply(functionalPredictions, length) == 0, NA_character_, toString(functionalPredictions)))
      } else {
        result_functionalPredictions <- result_df %>%
          dplyr::select(gene.symbol, variant, functionalPredictions) %>%
          tidyr::unnest(functionalPredictions, names_sep = '.', keep_empty = TRUE) %>%
          dplyr::rename("typeId" = "functionalPredictions.typeId",
                        "aggregatedScore" = "functionalPredictions.aggregatedScore")

        if ("functionalPredictions.tissues" %in% colnames(result_functionalPredictions)) {
          result_functionalPredictions <- result_functionalPredictions %>%
            tidyr::unnest(functionalPredictions.tissues, names_sep = '_', keep_empty = TRUE) %>%
            dplyr::rename("tissues_id" = "functionalPredictions.tissues_tissue.id",
                          "tissues_name" = "functionalPredictions.tissues_tissue.name")
          base::colnames(result_functionalPredictions) <- stringr::str_replace_all(colnames(result_functionalPredictions), "functionalPredictions.", "")
        }
      }

      result_pkg <- list(v2g = result_core, tssd = result_distances,
                         qtls = result_qtl, chromatin = result_intervals,
                         functionalpred = result_functionalPredictions)
    }
    cli::cli_progress_update()
    return(result_pkg)

  }, error = function(e) {
    # Handling connection timeout
    if(grepl("Timeout was reached", e$message)) {
      stop("Connection timeout reached while connecting to the Open Targets Genetics GraphQL API.")
    } else {
      stop(e) # Handle other types of errors
    }
  })
}

#' Query Open Targets Genetics API for variant information
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
  tryCatch({
    # Set up to query Open Targets Genetics API
    cli::cli_progress_step("Connecting to the Open Targets Genetics GraphQL API...", spinner = TRUE)
    otg_cli <- ghql::GraphqlClient$new(url = "https://api.genetics.opentargets.org/graphql")
    otg_qry <- ghql::Query$new()

    # Check variant id format
    if (grepl(pattern = "rs\\d+", variant_id)) {
      # Convert rs id to variant id
      query_searchid <- "query rsi2vid($queryString:String!) {
      search(queryString:$queryString){
        totalVariants
        variants{
          id
          }
        }
      }"

      variables <- list(queryString = variant_id)
      otg_qry$query(name = "rsi2vid", x = query_searchid)
      id_result <- jsonlite::fromJSON(otg_cli$exec(otg_qry$queries$rsi2vid, variables), flatten = TRUE)$data
      input_variant_id <- id_result$search$variants$id
    } else if (grepl(pattern = "\\d+_\\d+_[a-zA-Z]+_[a-zA-Z]+", variant_id)) {
      input_variant_id <- variant_id
    } else {
      stop("\nPlease provide a variant ID")
    }

    # Check if the input_variant_id is null or empty
    if (is.null(input_variant_id) || input_variant_id == "") {
      stop("There is no variant ID defined for this rsID by Open Target Genetics")
    }

    # Define the query
    query <- "query variantInfoquery($variantId: String!){
    variantInfo(variantId: $variantId){
      rsId
      id
      nearestGeneDistance
      nearestCodingGene{
        id
        symbol
      }
      nearestCodingGeneDistance
      mostSevereConsequence
      caddPhred
      gnomadNFE
      gnomadAFR
      gnomadAMR
      gnomadEAS
    }
  }"

    # Execute the query
    variables <- list(variantId = input_variant_id)
    otg_qry$query(name = "variantInfoquery", x = query)
    cli::cli_progress_step("Downloading data...", spinner = TRUE)
    var_info <- jsonlite::fromJSON(otg_cli$exec(otg_qry$queries$variantInfoquery, variables), flatten = TRUE)$data

    # Flatten and return data frame
    flat_var_info <- unlist(var_info$variantInfo)
    df_var_info <- as.data.frame(t(flat_var_info))
    names(df_var_info) <- names(flat_var_info)
    return(df_var_info)

  }, error = function(e) {
    # Handling connection timeout
    if(grepl("Timeout was reached", e$message)) {
      warning("Connection timeout reached while connecting to the Open Targets Genetics GraphQL API.")
    } else {
      warning(e) # Handle other types of errors
    }
    return(NULL)
  })
}
