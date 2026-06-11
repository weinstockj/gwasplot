#' Fixed-effects meta-analysis of two GWAS datasets
#'
#' @param x A `GWASFormatter` object, data.frame, or tibble containing standardized GWAS results.
#' @param y A second `GWASFormatter` object, data.frame, or tibble containing standardized GWAS results.
#' @param keep Output row-set policy. `"overlap"` (default) keeps only variants
#'   present in both studies (intersection). `"union"` keeps variants present in
#'   either study; a variant found in only one study passes through with that
#'   study's estimate and `N_studies = 1` (and `matched_by` of `"x_only"` or
#'   `"y_only"`).
#' @param allow_swaps Logical. If `TRUE`, simple `REF`/`ALT` swaps are harmonized by flipping the second study's `BETA`.
#'   `GWASFormatter` inputs require the DuckDB `stochastic` community extension to be installed and loadable.
#' @param ... Additional arguments passed to methods.
#' @return A GWAS-shaped meta-analysis result. `GWASFormatter` methods return a modified `GWASFormatter`; eager methods return a tibble.
#' @export
meta_analyze_fe <- function(x, y, ...) {
  UseMethod("meta_analyze_fe")
}

meta_required_columns_ <- c("CHROM", "POS", "REF", "ALT", "BETA", "SE")
meta_optional_columns_ <- c("ID", "PVALUE", "AF_ALT", "phenotype", "filename")

duckdb_stochastic_failure_message_ <- function(con, error) {
  duckdb_version <- tryCatch(
    DBI::dbGetQuery(con, "SELECT version() AS version")$version[[1]],
    error = function(version_error) NA_character_
  )

  version_note <- if (is.na(duckdb_version)) {
    ""
  } else {
    glue(" Current DuckDB version: {duckdb_version}.")
  }

  c(
    "DuckDB meta-analysis requires the `stochastic` community extension.",
    "Run `INSTALL stochastic FROM community; LOAD stochastic;` in DuckDB before calling `meta_analyze_fe()` on `GWASFormatter` inputs.",
    glue(
      "DuckDB could not install or load `stochastic` in this session.{version_note}"
    ),
    glue("Original DuckDB error: {conditionMessage(error)}")
  )
}

duckdb_check_stochastic_ <- function(con) {
  cached <- attr(con, "gwasplot_stochastic_available", exact = TRUE)

  if (!is.null(cached)) {
    return(cached)
  }

  status <- tryCatch(
    {
      DBI::dbExecute(con, "INSTALL stochastic FROM community")
      DBI::dbExecute(con, "LOAD stochastic")
      DBI::dbGetQuery(
        con,
        "SELECT dist_normal_cdf_complement(0.0, 1.0, 0.0) AS ok"
      )
      TRUE
    },
    error = function(e) {
      duckdb_stochastic_failure_message_(con, e)
    }
  )

  attr(con, "gwasplot_stochastic_available") <- status
  status
}

duckdb_require_stochastic_ <- function(con) {
  status <- duckdb_check_stochastic_(con)

  if (isTRUE(status)) {
    return(invisible(TRUE))
  }

  cli::cli_abort(status)
}

duckdb_check_stochastic_onload_ <- function() {
  if (
    !requireNamespace("DBI", quietly = TRUE) ||
      !requireNamespace("duckdb", quietly = TRUE)
  ) {
    return(invisible(FALSE))
  }

  tryCatch(
    {
      con <- DBI::dbConnect(duckdb::duckdb(), dbdir = ":memory:")
      on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE)

      status <- duckdb_check_stochastic_(con)
      .gwasplot_state$stochastic_status <- status

      if (!isTRUE(status)) {
        cli::cli_alert_warning(status[[1]])
      }

      invisible(isTRUE(status))
    },
    error = function(e) {
      .gwasplot_state$stochastic_status <- conditionMessage(e)
      cli::cli_alert_warning(
        "DuckDB stochastic preflight failed during package load."
      )
      invisible(FALSE)
    }
  )
}

meta_two_sided_pvalue_sql_ <- function(z_expr) {
  glue("2 * dist_normal_cdf_complement(0.0, 1.0, {z_expr})")
}

# DuckDB CTEs that normalize each study's alleles for joining.
meta_input_ctes_sql_ <- function(x_input_table, y_input_table) {
  glue(
    "
  WITH x_input AS (
    SELECT
      CHROM, POS, REF, ALT, ID, BETA, SE, AF_ALT,
      LEAST(REF, ALT) AS allele_a,
      GREATEST(REF, ALT) AS allele_b
    FROM {x_input_table}
  ),
  y_input AS (
    SELECT
      CHROM, POS, REF, ALT, ID, BETA, SE, AF_ALT,
      LEAST(REF, ALT) AS allele_a,
      GREATEST(REF, ALT) AS allele_b
    FROM {y_input_table}
  )"
  )
}

# Build the meta-analysis SQL. The overlap path inner-joins (both studies
# required); the union path full-outer-joins and lets the inverse-variance
# weighting drop a missing/invalid study to zero weight so single-study
# variants pass through with N_studies = 1.
build_meta_sql_ <- function(
  keep,
  result_table,
  x_input_table,
  y_input_table,
  swap_condition,
  pvalue_sql
) {
  ctes = meta_input_ctes_sql_(x_input_table, y_input_table)

  if (keep == "overlap") {
    return(glue(
      "
  CREATE OR REPLACE TABLE {result_table} AS
  {ctes},
  matched AS (
    SELECT
      x.CHROM,
      x.POS,
      x.REF AS REF,
      x.ALT AS ALT,
      x.ID AS ID,
      x.AF_ALT AS AF_ALT,
      x.BETA AS BETA_1,
      x.SE AS SE_1,
      y.BETA AS BETA_2,
      y.SE AS SE_2,
      CASE
        WHEN x.REF = y.REF AND x.ALT = y.ALT THEN 'exact'
        WHEN {swap_condition} THEN 'swap'
        ELSE 'mismatch'
      END AS matched_by,
      CASE WHEN {swap_condition} THEN TRUE ELSE FALSE END AS flip_2,
      CASE WHEN {swap_condition} THEN -y.BETA ELSE y.BETA END AS BETA_2_harmonized
    FROM x_input x
    INNER JOIN y_input y
      ON x.CHROM = y.CHROM
     AND x.POS = y.POS
     AND x.allele_a = y.allele_a
     AND x.allele_b = y.allele_b
  ),
  meta AS (
    SELECT
      CHROM, POS, REF, ALT, ID, AF_ALT,
      ((1 / (SE_1 * SE_1)) * BETA_1 + (1 / (SE_2 * SE_2)) * BETA_2_harmonized) /
        ((1 / (SE_1 * SE_1)) + (1 / (SE_2 * SE_2))) AS BETA,
      SQRT(1 / ((1 / (SE_1 * SE_1)) + (1 / (SE_2 * SE_2)))) AS SE,
      matched_by, flip_2, BETA_1, SE_1, BETA_2, SE_2,
      abs(((1 / (SE_1 * SE_1)) * BETA_1 + (1 / (SE_2 * SE_2)) * BETA_2_harmonized) /
        ((1 / (SE_1 * SE_1)) + (1 / (SE_2 * SE_2))) /
        SQRT(1 / ((1 / (SE_1 * SE_1)) + (1 / (SE_2 * SE_2))))) AS z_abs
    FROM matched
    WHERE matched_by != 'mismatch'
      AND SE_1 IS NOT NULL AND SE_2 IS NOT NULL
      AND BETA_1 IS NOT NULL AND BETA_2_harmonized IS NOT NULL
      AND SE_1 > 0 AND SE_2 > 0
  )
  SELECT
    CHROM, POS, REF, ALT, ID, AF_ALT, BETA, SE,
    {pvalue_sql} AS PVALUE,
    2 AS N_studies,
    matched_by, flip_2, BETA_1, SE_1, BETA_2, SE_2
  FROM meta
  "
    ))
  }

  glue(
    "
  CREATE OR REPLACE TABLE {result_table} AS
  {ctes},
  matched AS (
    SELECT
      COALESCE(x.CHROM, y.CHROM) AS CHROM,
      COALESCE(x.POS, y.POS) AS POS,
      CASE WHEN x.REF IS NOT NULL THEN x.REF ELSE y.REF END AS REF,
      CASE WHEN x.REF IS NOT NULL THEN x.ALT ELSE y.ALT END AS ALT,
      CASE WHEN x.REF IS NOT NULL THEN x.ID ELSE y.ID END AS ID,
      CASE WHEN x.REF IS NOT NULL THEN x.AF_ALT ELSE y.AF_ALT END AS AF_ALT,
      x.BETA AS BETA_1,
      x.SE AS SE_1,
      y.BETA AS BETA_2,
      y.SE AS SE_2,
      (x.REF IS NOT NULL) AS present_1,
      (y.REF IS NOT NULL) AS present_2,
      CASE
        WHEN x.REF IS NOT NULL AND y.REF IS NOT NULL AND x.REF = y.REF AND x.ALT = y.ALT THEN 'exact'
        WHEN x.REF IS NOT NULL AND y.REF IS NOT NULL AND {swap_condition} THEN 'swap'
        WHEN x.REF IS NOT NULL AND y.REF IS NOT NULL THEN 'mismatch'
        WHEN x.REF IS NOT NULL THEN 'x_only'
        ELSE 'y_only'
      END AS matched_by,
      CASE
        WHEN x.REF IS NOT NULL AND y.REF IS NOT NULL AND {swap_condition} THEN TRUE
        ELSE FALSE
      END AS flip_2,
      CASE
        WHEN y.REF IS NULL THEN NULL
        WHEN x.REF IS NOT NULL AND {swap_condition} THEN -y.BETA
        ELSE y.BETA
      END AS BETA_2_harmonized
    FROM x_input x
    FULL OUTER JOIN y_input y
      ON x.CHROM = y.CHROM
     AND x.POS = y.POS
     AND x.allele_a = y.allele_a
     AND x.allele_b = y.allele_b
  ),
  weighted AS (
    SELECT
      *,
      CASE
        WHEN present_1 AND SE_1 IS NOT NULL AND BETA_1 IS NOT NULL AND SE_1 > 0
        THEN 1 / (SE_1 * SE_1) ELSE 0
      END AS w1,
      CASE
        WHEN present_2 AND matched_by != 'mismatch'
         AND SE_2 IS NOT NULL AND BETA_2_harmonized IS NOT NULL AND SE_2 > 0
        THEN 1 / (SE_2 * SE_2) ELSE 0
      END AS w2
    FROM matched
  ),
  meta AS (
    SELECT
      CHROM, POS, REF, ALT, ID, AF_ALT,
      ((CASE WHEN w1 > 0 THEN w1 * BETA_1 ELSE 0 END) +
       (CASE WHEN w2 > 0 THEN w2 * BETA_2_harmonized ELSE 0 END)) / (w1 + w2) AS BETA,
      SQRT(1 / (w1 + w2)) AS SE,
      matched_by, flip_2, BETA_1, SE_1, BETA_2, SE_2,
      ((CASE WHEN w1 > 0 THEN 1 ELSE 0 END) + (CASE WHEN w2 > 0 THEN 1 ELSE 0 END)) AS N_studies,
      abs(((CASE WHEN w1 > 0 THEN w1 * BETA_1 ELSE 0 END) +
           (CASE WHEN w2 > 0 THEN w2 * BETA_2_harmonized ELSE 0 END)) / (w1 + w2) /
          SQRT(1 / (w1 + w2))) AS z_abs
    FROM weighted
    WHERE (w1 + w2) > 0
  )
  SELECT
    CHROM, POS, REF, ALT, ID, AF_ALT, BETA, SE,
    {pvalue_sql} AS PVALUE,
    N_studies,
    matched_by, flip_2, BETA_1, SE_1, BETA_2, SE_2
  FROM meta
  "
  )
}

prepare_meta_input_ <- function(x, label) {
  input <- tibble::as_tibble(x)
  missing_cols <- setdiff(meta_required_columns_, names(input))

  if (length(missing_cols) > 0) {
    cli::cli_abort(
      "{label} is missing required columns: {toString(missing_cols)}"
    )
  }

  keep_cols <- unique(c(meta_required_columns_, meta_optional_columns_))
  input <- input[, intersect(keep_cols, names(input)), drop = FALSE]

  input$CHROM <- as.character(input$CHROM)
  input$POS <- as.integer(input$POS)
  input$REF <- toupper(as.character(input$REF))
  input$ALT <- toupper(as.character(input$ALT))

  if (!"ID" %in% names(input)) {
    input <- add_ID(input)
  }

  input$.meta_ref <- pmin(input$REF, input$ALT)
  input$.meta_alt <- pmax(input$REF, input$ALT)

  duplicate_keys <- duplicated(input[c(
    "CHROM",
    "POS",
    ".meta_ref",
    ".meta_alt"
  )])
  if (any(duplicate_keys)) {
    cli::cli_abort(
      "{label} contains duplicated variants after allele harmonization."
    )
  }

  input
}

harmonize_meta_inputs_ <- function(x, y, allow_swaps = TRUE, keep = "overlap") {
  left <- prepare_meta_input_(x, "x")
  right <- prepare_meta_input_(y, "y")

  join_fn <- if (keep == "union") dplyr::full_join else dplyr::inner_join
  joined <- join_fn(
    left,
    right,
    by = c("CHROM", "POS", ".meta_ref", ".meta_alt"),
    suffix = c("_1", "_2")
  )

  if (nrow(joined) == 0) {
    cli::cli_abort("No variants were found between x and y.")
  }

  present_1 <- !is.na(joined$REF_1)
  present_2 <- !is.na(joined$REF_2)
  both <- present_1 & present_2

  exact_match <- both &
    joined$REF_1 == joined$REF_2 &
    joined$ALT_1 == joined$ALT_2
  swapped_match <- both &
    joined$REF_1 == joined$ALT_2 &
    joined$ALT_1 == joined$REF_2
  exact_match <- !is.na(exact_match) & exact_match
  swapped_match <- !is.na(swapped_match) & swapped_match

  if (!allow_swaps) {
    swapped_match[] <- FALSE
  }

  joined$matched_by <- dplyr::case_when(
    exact_match ~ "exact",
    swapped_match ~ "swap",
    both ~ "mismatch",
    present_1 ~ "x_only",
    TRUE ~ "y_only"
  )

  # For overlap, drop variants whose alleles cannot be reconciled. For union,
  # keep them (their study-2 contribution is dropped downstream via weighting).
  if (keep == "overlap") {
    joined <- joined[joined$matched_by != "mismatch", , drop = FALSE]

    if (nrow(joined) == 0) {
      cli::cli_abort("No harmonized variants remained after allele matching.")
    }
  }

  joined$present_1 <- !is.na(joined$REF_1)
  joined$present_2 <- !is.na(joined$REF_2)
  joined$flip_2 <- joined$matched_by == "swap"

  # Report each variant in study-1 orientation when present, else study-2's.
  joined$REF <- ifelse(joined$present_1, joined$REF_1, joined$REF_2)
  joined$ALT <- ifelse(joined$present_1, joined$ALT_1, joined$ALT_2)
  joined$ID <- ifelse(joined$present_1, joined$ID_1, joined$ID_2)

  joined$BETA_2_harmonized <- ifelse(
    joined$flip_2,
    -joined$BETA_2,
    joined$BETA_2
  )

  if ("AF_ALT_1" %in% names(joined) || "AF_ALT_2" %in% names(joined)) {
    af_1 <- if ("AF_ALT_1" %in% names(joined)) joined$AF_ALT_1 else NA_real_
    af_2 <- if ("AF_ALT_2" %in% names(joined)) joined$AF_ALT_2 else NA_real_
    joined$AF_ALT <- ifelse(joined$present_1, af_1, af_2)
  }

  joined
}

compute_fixed_effects_meta_ <- function(harmonized, keep = "overlap") {
  valid_1 <-
    harmonized$present_1 &
    is.finite(harmonized$BETA_1) &
    is.finite(harmonized$SE_1) &
    harmonized$SE_1 > 0
  valid_2 <-
    harmonized$present_2 &
    is.finite(harmonized$BETA_2_harmonized) &
    is.finite(harmonized$SE_2) &
    harmonized$SE_2 > 0
  valid_1 <- !is.na(valid_1) & valid_1
  valid_2 <- !is.na(valid_2) & valid_2

  keep_rows <- if (keep == "union") valid_1 | valid_2 else valid_1 & valid_2

  harmonized <- harmonized[keep_rows, , drop = FALSE]
  valid_1 <- valid_1[keep_rows]
  valid_2 <- valid_2[keep_rows]

  if (nrow(harmonized) == 0) {
    cli::cli_abort(
      "No variants with finite BETA and SE values were available for meta-analysis."
    )
  }

  # Inverse-variance weights; a study with weight 0 drops out of the average,
  # so a variant present in only one study returns that study's estimate.
  weight_1 <- ifelse(valid_1, 1 / (harmonized$SE_1^2), 0)
  weight_2 <- ifelse(valid_2, 1 / (harmonized$SE_2^2), 0)
  contrib_1 <- ifelse(valid_1, weight_1 * harmonized$BETA_1, 0)
  contrib_2 <- ifelse(valid_2, weight_2 * harmonized$BETA_2_harmonized, 0)

  beta <- (contrib_1 + contrib_2) / (weight_1 + weight_2)
  se <- sqrt(1 / (weight_1 + weight_2))
  z_score <- beta / se
  pvalue <- 2 * stats::pnorm(-abs(z_score))
  n_studies <- as.integer(valid_1) + as.integer(valid_2)

  result <- tibble::tibble(
    CHROM = harmonized$CHROM,
    POS = as.integer(harmonized$POS),
    REF = harmonized$REF,
    ALT = harmonized$ALT,
    ID = harmonized$ID,
    BETA = beta,
    SE = se,
    PVALUE = pvalue,
    N_studies = n_studies,
    matched_by = harmonized$matched_by,
    flip_2 = harmonized$flip_2,
    BETA_1 = harmonized$BETA_1,
    SE_1 = harmonized$SE_1,
    BETA_2 = harmonized$BETA_2,
    SE_2 = harmonized$SE_2
  )

  if ("AF_ALT" %in% names(harmonized)) {
    result$AF_ALT <- harmonized$AF_ALT
    result <- result[, c(
      "CHROM",
      "POS",
      "REF",
      "ALT",
      "ID",
      "AF_ALT",
      "BETA",
      "SE",
      "PVALUE",
      "N_studies",
      "matched_by",
      "flip_2",
      "BETA_1",
      "SE_1",
      "BETA_2",
      "SE_2"
    )]
  }

  result
}

materialize_meta_input_table_ <- function(gwas, con, table_name, select_cols) {
  if (DBI::dbExistsTable(con, gwas$table_name)) {
    dplyr::tbl(con, gwas$table_name) %>%
      dplyr::select(dplyr::any_of(select_cols)) %>%
      dplyr::compute(name = table_name, temporary = FALSE, overwrite = TRUE)
  } else {
    dplyr::copy_to(
      con,
      dplyr::collect(dplyr::select(gwas$data, dplyr::any_of(select_cols))),
      name = table_name,
      temporary = FALSE,
      overwrite = TRUE
    )
  }

  index_name = sanitize_table_name_(
    paste(table_name, "idx", sep = "_"),
    paste0(table_name, "_idx")
  )
  DBI::dbExecute(
    con,
    glue(
      "CREATE INDEX IF NOT EXISTS {index_name} ON {table_name} (CHROM, POS, REF, ALT)"
    )
  )

  invisible(table_name)
}

#' @export
meta_analyze_fe.data.frame <- function(
  x,
  y,
  keep = c("overlap", "union"),
  allow_swaps = TRUE,
  ...
) {
  if (!inherits(y, c("data.frame", "tbl_df"))) {
    cli::cli_abort(
      "For data.frame inputs, y must also be a data.frame or tibble."
    )
  }

  keep <- match.arg(keep)

  harmonized <- harmonize_meta_inputs_(
    x,
    y,
    allow_swaps = allow_swaps,
    keep = keep
  )
  compute_fixed_effects_meta_(harmonized, keep = keep)
}

#' @export
meta_analyze_fe.tbl_df <- function(
  x,
  y,
  keep = c("overlap", "union"),
  allow_swaps = TRUE,
  ...
) {
  meta_analyze_fe.data.frame(x, y, keep = keep, allow_swaps = allow_swaps, ...)
}

#' @export
meta_analyze_fe.GWASFormatter <- function(
  x,
  y,
  keep = c("overlap", "union"),
  allow_swaps = TRUE,
  ...
) {
  keep <- match.arg(keep)
  select_cols <- c(
    "CHROM",
    "POS",
    "REF",
    "ALT",
    "ID",
    "BETA",
    "SE",
    "PVALUE",
    "AF_ALT"
  )

  if (!inherits(y, c("GWASFormatter", "data.frame", "tbl_df"))) {
    cli::cli_abort(
      "For GWASFormatter inputs, y must be a GWASFormatter, data.frame, or tibble."
    )
  }

  if (!inherits(y, "GWASFormatter")) {
    result = meta_analyze_fe.data.frame(
      dplyr::collect(dplyr::select(x$data, dplyr::any_of(select_cols))),
      tibble::as_tibble(y),
      keep = keep,
      allow_swaps = allow_swaps,
      ...
    )

    result_table = make_unique_table_name_(x$table_name, "meta")
    dplyr::copy_to(
      x$con,
      result,
      name = result_table,
      temporary = FALSE,
      overwrite = TRUE
    )
    x$table_name = result_table
    x$source_table_name = result_table
    x$data = dplyr::tbl(x$con, result_table)
    return(x)
  }

  x_input_table = make_unique_table_name_(x$table_name, "meta_input")
  y_input_table = make_unique_table_name_(y$table_name, "meta_input")
  result_table = make_unique_table_name_(x$table_name, "meta")
  duckdb_require_stochastic_(x$con)

  materialize_meta_input_table_(x, x$con, x_input_table, select_cols)
  materialize_meta_input_table_(y, x$con, y_input_table, select_cols)

  swap_condition = if (allow_swaps) {
    "x.REF = y.ALT AND x.ALT = y.REF"
  } else {
    "FALSE"
  }

  pvalue_sql = meta_two_sided_pvalue_sql_("z_abs")

  sql = build_meta_sql_(
    keep = keep,
    result_table = result_table,
    x_input_table = x_input_table,
    y_input_table = y_input_table,
    swap_condition = swap_condition,
    pvalue_sql = pvalue_sql
  )

  DBI::dbExecute(x$con, sql)

  x$table_name = result_table
  x$source_table_name = result_table
  x$data = dplyr::tbl(x$con, result_table)
  x
}
