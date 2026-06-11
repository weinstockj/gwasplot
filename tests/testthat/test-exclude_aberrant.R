test_that("exclude_aberrant_pvalue_loci drops lonely peaks but keeps supported leads", {
  df <- tibble::tibble(
    CHROM = c("chr1", "chr1", "chr2", "chr2", "chr3"),
    POS = c(1000000L, 2000000L, 5000000L, 5050000L, 9000000L),
    ID = c("v1", "v2", "v3", "v4", "v5"),
    PVALUE = c(1e-12, 0.4, 1e-12, 1e-6, 1e-9)
  )

  result <- exclude_aberrant_pvalue_loci(df)

  # v1: strong lead on chr1 with nearest neighbour 1 Mb away (> 500 kb) -> dropped.
  # v3: strong lead on chr2 with a supporter (v4, 1e-6) 50 kb away         -> kept.
  # v5: 1e-9 is above the lead threshold (5e-10), so never examined        -> kept.
  expect_equal(sort(result$ID), c("v2", "v3", "v4", "v5"))
  expect_false("v1" %in% result$ID)
})

test_that("exclude_aberrant_pvalue_loci respects custom thresholds and window", {
  df <- tibble::tibble(
    CHROM = c("chr1", "chr1"),
    POS = c(1000000L, 1100000L),
    ID = c("lead", "neighbour"),
    PVALUE = c(1e-12, 1e-4)
  )

  # Default support threshold 5e-5: neighbour (1e-4) does NOT support -> dropped.
  expect_equal(nrow(exclude_aberrant_pvalue_loci(df)), 1L)
  expect_false("lead" %in% exclude_aberrant_pvalue_loci(df)$ID)

  # Looser support threshold 1e-3: neighbour (1e-4) now supports -> kept.
  kept <- exclude_aberrant_pvalue_loci(df, support_pvalue_threshold = 1e-3)
  expect_equal(nrow(kept), 2L)

  # Tight window 50 kb: neighbour is 100 kb away, out of window -> dropped.
  tight <- exclude_aberrant_pvalue_loci(
    df,
    support_pvalue_threshold = 1e-3,
    window_kb = 50
  )
  expect_false("lead" %in% tight$ID)
})

test_that("exclude_aberrant_pvalue_loci is a no-op when no strong leads exist", {
  df <- tibble::tibble(
    CHROM = c("chr1", "chr1"),
    POS = c(1000000L, 2000000L),
    ID = c("a", "b"),
    PVALUE = c(1e-7, 0.5)
  )

  expect_equal(exclude_aberrant_pvalue_loci(df), df)
})

test_that("exclude_aberrant_pvalue_loci.GWASFormatter drops aberrant leads from the table", {
  con <- DBI::dbConnect(duckdb::duckdb(), dbdir = ":memory:")
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE)

  workdir <- tempfile("aberrant-")
  dir.create(workdir)
  old_wd <- setwd(workdir)
  on.exit(setwd(old_wd), add = TRUE)

  file1 <- tempfile(fileext = ".csv")

  # LOG10P = -log10(PVALUE). Lonely peak at chr1:1e6 (1e-12, no support within
  # 500 kb); supported lead at chr2:5e6 (1e-12) with neighbour at chr2:5.05e6.
  utils::write.csv(
    data.frame(
      CHROM = c("chr1", "chr1", "chr2", "chr2"),
      POS = c(1000000L, 2000000L, 5000000L, 5050000L),
      ALLELE0 = c("A", "A", "A", "A"),
      ALLELE1 = c("C", "C", "C", "C"),
      A1FREQ = c(0.1, 0.2, 0.3, 0.4),
      BETA = c(0.5, 0.01, 0.5, 0.2),
      SE = c(0.05, 0.05, 0.05, 0.05),
      LOG10P = c(12, 0.4, 12, 6)
    ),
    file1,
    row.names = FALSE,
    quote = TRUE
  )

  gwas <- reformat_summary_statistics(file1, use_cache = FALSE)
  on.exit(
    if (DBI::dbIsValid(gwas$con)) DBI::dbDisconnect(gwas$con, shutdown = TRUE),
    add = TRUE
  )

  exclude_aberrant_pvalue_loci(gwas)

  kept <- gwas$data %>%
    dplyr::select(CHROM, POS) %>%
    dplyr::arrange(CHROM, POS) %>%
    dplyr::collect()

  # The lonely peak at chr1:1,000,000 is removed; everything else remains.
  expect_equal(nrow(kept), 3L)
  expect_false(any(kept$CHROM == "chr1" & kept$POS == 1000000L))
  expect_true(any(kept$CHROM == "chr2" & kept$POS == 5000000L))
})
