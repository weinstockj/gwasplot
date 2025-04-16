#' Pull top hits from a gwas object
#' 
#' @param gwas A gwas object containing the data to plot.
#' @param threshold The p-value threshold to filter the top hits. Default is 5e-8.
#' @return A data frame with the top hits.
#' @export
annotate_top_hits = function(gwas, threshold = 5e-8) {
  df = gwas$data %>%
    dplyr::filter(PVALUE < threshold) %>%
    dplyr::collect(.)

  # possibly_query = purrr::possibly(
  #   query_ot_api_variants,
  #   otherwise = NULL,
  #   quiet = TRUE
  # )
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

#' Find the nearest gene for each top hit
#' 
#' @param top_hits A data frame containing the top hits.
#' @param threshold The distance threshold to consider a gene as nearest. Default is 1e5.
#' @return A data frame with the nearest gene information.
#' @export
find_nearest_gene = function(top_hits, threshold = 1e5) {

  con = db_connect()


  dplyr::copy_to(
    con,
    top_hits,
    name = "top_hits",
    temporary = FALSE,
    overwrite = TRUE
  )

  dplyr::copy_to(
    con,
    human_genes,
    name = "human_genes",
    temporary = FALSE,
    overwrite = TRUE
  )

  sql = glue("
  WITH VariantCross AS (
      SELECT
       t.*,
       LEAST(ABS(t.pos - g.start), ABS(t.pos - g.end)) AS distance,
       g.gene_name,
       g.gene_id
      FROM top_hits t
      CROSS JOIN human_genes g
      WHERE t.chrom = g.chrom
      AND g.gene_biotype = 'protein_coding'
  ),
  VariantCrossRanked AS (
      SELECT *,
      ROW_NUMBER() OVER (PARTITION BY ID ORDER BY distance) AS rn
      FROM VariantCross
  )
  SELECT
    *
  FROM VariantCrossRanked
  WHERE rn = 1
  ")

  DBI::dbGetQuery(con, sql) %>%
    tibble::as_tibble(.) %>%
    dplyr::mutate(
      gene_id = dplyr::case_when(
        distance > threshold ~ NA_character_,
        TRUE ~ gene_id
      ),
      gene_name = dplyr::case_when(
        distance > threshold ~ NA_character_,
        TRUE ~ gene_name
      ),
      distance = dplyr::case_when(
        distance > threshold ~ NA_real_,
        TRUE ~ distance
      )
    ) %>%
    dplyr::select(-rn)
}

#' Annotate top hits with centromere information
#' 
#' @param top_hits A data frame containing the top hits.
#' @return A data frame with the centromere information.
#' @export
annotate_with_centromere = function(top_hits) {
  con = db_connect()

  dplyr::copy_to(
    con,
    top_hits,
    name = "top_hits",
    temporary = FALSE,
    overwrite = TRUE
  )

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

  sql = glue("
  SELECT
    t.*,
    c.name as ideogram_name,
    c.in_centromere
  FROM top_hits t
  LEFT JOIN ideogram c ON
  (t.chrom = c.chrom) AND (t.POS BETWEEN c.start AND c.end)
  ")

  DBI::dbGetQuery(con, sql) %>%
    tibble::as_tibble(.) %>%
    dplyr::mutate(
      in_centromere = dplyr::case_when(
        is.na(in_centromere) ~ FALSE,
        TRUE ~ in_centromere
      )
    )
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
    dplyr::distinct(.)

  dplyr::copy_to(
    con,
    chip_genes %>%
      dplyr::mutate(
        is_chip_gene = TRUE
      ),
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
