test_that("meta_analyze_fe harmonizes swapped alleles and computes fixed-effects results", {
  study1 <- tibble::tibble(
    CHROM = c("chr1", "chr1"),
    POS = c(100L, 200L),
    REF = c("A", "A"),
    ALT = c("C", "G"),
    BETA = c(0.1, 0.2),
    SE = c(0.1, 0.1),
    PVALUE = c(0.01, 0.01)
  ) %>%
    add_ID()

  study2 <- tibble::tibble(
    CHROM = c("chr1", "chr1"),
    POS = c(100L, 200L),
    REF = c("A", "G"),
    ALT = c("C", "A"),
    BETA = c(0.2, -0.1),
    SE = c(0.2, 0.2),
    PVALUE = c(0.02, 0.02)
  ) %>%
    add_ID()

  result <- meta_analyze_fe(study1, study2)

  expect_equal(nrow(result), 2L)
  expect_equal(result$matched_by, c("exact", "swap"))
  expect_equal(result$flip_2, c(FALSE, TRUE))
  expect_equal(result$BETA, c(0.12, 0.18), tolerance = 1e-8)
  expect_equal(result$SE, c(sqrt(1 / 125), sqrt(1 / 125)), tolerance = 1e-8)
  expect_true(all(
    c("CHROM", "POS", "REF", "ALT", "ID", "BETA", "SE", "PVALUE") %in%
      names(result)
  ))
})

test_that("meta_analyze_fe keep = 'union' passes through single-study variants", {
  study1 <- tibble::tibble(
    CHROM = c("chr1", "chr1"),
    POS = c(100L, 200L),
    REF = c("A", "A"),
    ALT = c("C", "G"),
    BETA = c(0.1, 0.2),
    SE = c(0.1, 0.1),
    PVALUE = c(0.01, 0.01)
  ) %>%
    add_ID()

  # study2 shares variant at POS 100 (swapped) and adds a unique one at POS 300.
  study2 <- tibble::tibble(
    CHROM = c("chr1", "chr1"),
    POS = c(100L, 300L),
    REF = c("A", "T"),
    ALT = c("C", "G"),
    BETA = c(0.2, 0.4),
    SE = c(0.2, 0.1),
    PVALUE = c(0.02, 0.02)
  ) %>%
    add_ID()

  overlap <- meta_analyze_fe(study1, study2, keep = "overlap")
  expect_equal(nrow(overlap), 1L)
  expect_equal(overlap$POS, 100L)

  union <- meta_analyze_fe(study1, study2, keep = "union")
  union <- union[order(union$POS), ]

  # Union = study1's 2 variants + study2's unique POS 300 = 3 rows.
  expect_equal(union$POS, c(100L, 200L, 300L))
  expect_equal(union$matched_by, c("exact", "x_only", "y_only"))
  expect_equal(union$N_studies, c(2L, 1L, 1L))

  # POS 100: both studies, inverse-variance weighted (w1 = 100, w2 = 25).
  expect_equal(union$BETA[[1]], 0.12, tolerance = 1e-8)
  expect_equal(union$SE[[1]], sqrt(1 / 125), tolerance = 1e-8)

  # POS 200: study1 only -> study1's estimate verbatim.
  expect_equal(union$BETA[[2]], 0.2, tolerance = 1e-8)
  expect_equal(union$SE[[2]], 0.1, tolerance = 1e-8)

  # POS 300: study2 only -> study2's estimate verbatim.
  expect_equal(union$BETA[[3]], 0.4, tolerance = 1e-8)
  expect_equal(union$SE[[3]], 0.1, tolerance = 1e-8)
})

test_that("meta_analyze_fe.GWASFormatter returns a GWAS-compatible result", {
  con <- DBI::dbConnect(duckdb::duckdb(), dbdir = ":memory:")
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE)

  tryCatch(
    duckdb_require_stochastic_(con),
    error = function(e) {
      testthat::skip(conditionMessage(e))
    }
  )

  workdir <- tempfile("meta-gwas-")
  dir.create(workdir)

  old_wd <- setwd(workdir)
  on.exit(setwd(old_wd), add = TRUE)

  file1 <- tempfile(fileext = ".csv")
  file2 <- tempfile(fileext = ".csv")

  utils::write.csv(
    data.frame(
      CHROM = c("chr1", "chr1"),
      POS = c(100L, 200L),
      ALLELE0 = c("A", "A"),
      ALLELE1 = c("C", "G"),
      A1FREQ = c(0.1, 0.2),
      BETA = c(0.1, 0.2),
      SE = c(0.1, 0.1),
      LOG10P = c(2, 2)
    ),
    file1,
    row.names = FALSE,
    quote = TRUE
  )

  utils::write.csv(
    data.frame(
      CHROM = c("chr1", "chr1"),
      POS = c(100L, 200L),
      ALLELE0 = c("A", "G"),
      ALLELE1 = c("C", "A"),
      A1FREQ = c(0.1, 0.8),
      BETA = c(0.2, -0.1),
      SE = c(0.2, 0.2),
      LOG10P = c(2, 2)
    ),
    file2,
    row.names = FALSE,
    quote = TRUE
  )

  gwas1 <- reformat_summary_statistics(file1, use_cache = FALSE)
  gwas2 <- reformat_summary_statistics(file2, use_cache = FALSE)

  expect_false(identical(gwas1$table_name, gwas2$table_name))

  on.exit(
    {
      if (DBI::dbIsValid(gwas1$con)) {
        DBI::dbDisconnect(gwas1$con, shutdown = TRUE)
      }

      if (DBI::dbIsValid(gwas2$con)) {
        DBI::dbDisconnect(gwas2$con, shutdown = TRUE)
      }
    },
    add = TRUE
  )

  meta <- meta_analyze_fe(gwas1, gwas2)
  result <- meta$data %>%
    dplyr::select(CHROM, POS, BETA, SE, PVALUE, matched_by, flip_2) %>%
    dplyr::collect()

  expect_equal(nrow(result), 2L)
  expect_equal(result$matched_by, c("exact", "swap"))
  expect_equal(result$flip_2, c(FALSE, TRUE))
  expect_equal(result$BETA, c(0.12, 0.18), tolerance = 1e-8)

  selected <- select_top_hits(result, threshold = 1)
  expect_equal(nrow(selected), 2L)
})

test_that("meta_analyze_fe.GWASFormatter supports keep = 'union'", {
  con <- DBI::dbConnect(duckdb::duckdb(), dbdir = ":memory:")
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE)

  tryCatch(
    duckdb_require_stochastic_(con),
    error = function(e) {
      testthat::skip(conditionMessage(e))
    }
  )

  workdir <- tempfile("meta-gwas-union-")
  dir.create(workdir)

  old_wd <- setwd(workdir)
  on.exit(setwd(old_wd), add = TRUE)

  file1 <- tempfile(fileext = ".csv")
  file2 <- tempfile(fileext = ".csv")

  # Shared variant at POS 100; study1-only at POS 200; study2-only at POS 300.
  utils::write.csv(
    data.frame(
      CHROM = c("chr1", "chr1"),
      POS = c(100L, 200L),
      ALLELE0 = c("A", "A"),
      ALLELE1 = c("C", "G"),
      A1FREQ = c(0.1, 0.2),
      BETA = c(0.1, 0.2),
      SE = c(0.1, 0.1),
      LOG10P = c(2, 2)
    ),
    file1,
    row.names = FALSE,
    quote = TRUE
  )

  utils::write.csv(
    data.frame(
      CHROM = c("chr1", "chr1"),
      POS = c(100L, 300L),
      ALLELE0 = c("A", "T"),
      ALLELE1 = c("C", "G"),
      A1FREQ = c(0.1, 0.3),
      BETA = c(0.2, 0.4),
      SE = c(0.2, 0.1),
      LOG10P = c(2, 2)
    ),
    file2,
    row.names = FALSE,
    quote = TRUE
  )

  gwas1 <- reformat_summary_statistics(file1, use_cache = FALSE)
  gwas2 <- reformat_summary_statistics(file2, use_cache = FALSE)

  on.exit(
    {
      if (DBI::dbIsValid(gwas1$con)) {
        DBI::dbDisconnect(gwas1$con, shutdown = TRUE)
      }

      if (DBI::dbIsValid(gwas2$con)) {
        DBI::dbDisconnect(gwas2$con, shutdown = TRUE)
      }
    },
    add = TRUE
  )

  meta <- meta_analyze_fe(gwas1, gwas2, keep = "union")
  result <- meta$data %>%
    dplyr::select(POS, BETA, SE, N_studies, matched_by) %>%
    dplyr::arrange(POS) %>%
    dplyr::collect()

  expect_equal(result$POS, c(100L, 200L, 300L))
  expect_equal(result$matched_by, c("exact", "x_only", "y_only"))
  expect_equal(result$N_studies, c(2L, 1L, 1L))
  expect_equal(result$BETA, c(0.12, 0.2, 0.4), tolerance = 1e-8)
  expect_equal(result$SE, c(sqrt(1 / 125), 0.1, 0.1), tolerance = 1e-8)
})

test_that("meta_analyze_fe.GWASFormatter keeps extreme-tail p-values finite", {
  con <- DBI::dbConnect(duckdb::duckdb(), dbdir = ":memory:")
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE)

  tryCatch(
    duckdb_require_stochastic_(con),
    error = function(e) {
      testthat::skip(conditionMessage(e))
    }
  )

  workdir <- tempfile("meta-gwas-tail-")
  dir.create(workdir)

  old_wd <- setwd(workdir)
  on.exit(setwd(old_wd), add = TRUE)

  file1 <- tempfile(fileext = ".csv")
  file2 <- tempfile(fileext = ".csv")

  utils::write.csv(
    data.frame(
      CHROM = "chr1",
      POS = 100L,
      ALLELE0 = "A",
      ALLELE1 = "C",
      A1FREQ = 0.1,
      BETA = 1.0,
      SE = 0.1,
      LOG10P = 2
    ),
    file1,
    row.names = FALSE,
    quote = TRUE
  )

  utils::write.csv(
    data.frame(
      CHROM = "chr1",
      POS = 100L,
      ALLELE0 = "A",
      ALLELE1 = "C",
      A1FREQ = 0.1,
      BETA = 1.0,
      SE = 0.1,
      LOG10P = 2
    ),
    file2,
    row.names = FALSE,
    quote = TRUE
  )

  gwas1 <- reformat_summary_statistics(file1, use_cache = FALSE)
  gwas2 <- reformat_summary_statistics(file2, use_cache = FALSE)

  on.exit(
    {
      if (DBI::dbIsValid(gwas1$con)) {
        DBI::dbDisconnect(gwas1$con, shutdown = TRUE)
      }

      if (DBI::dbIsValid(gwas2$con)) {
        DBI::dbDisconnect(gwas2$con, shutdown = TRUE)
      }
    },
    add = TRUE
  )

  meta <- meta_analyze_fe(gwas1, gwas2)
  result <- meta$data %>%
    dplyr::select(BETA, SE, PVALUE) %>%
    dplyr::collect()

  expected_p <- 2 * stats::pnorm(-abs(result$BETA / result$SE))

  expect_gt(result$PVALUE[[1]], 0)
  expect_equal(result$PVALUE, expected_p, tolerance = 1e-6)
})
