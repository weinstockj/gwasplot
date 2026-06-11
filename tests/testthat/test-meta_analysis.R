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
  expect_true(all(c("CHROM", "POS", "REF", "ALT", "ID", "BETA", "SE", "PVALUE") %in% names(result)))
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