test_that("find_nearest_gene.tbl_df handles custom columns and duplicate rows", {
  gene <- human_genes %>%
    dplyr::filter(gene_biotype == "protein_coding", !is.na(gene_name)) %>%
    dplyr::slice(1)

  variants <- tibble::tibble(
    chr = c(gene$chrom, gene$chrom, gene$chrom, gene$chrom),
    pos = c(gene$start, gene$start - 10, gene$end + 10, gene$start),
    variant_id = c("inside", "upstream", "downstream", "inside")
  )

  result <- find_nearest_gene(
    variants,
    threshold = 100000,
    chrom_col = "chr",
    pos_col = "pos",
    id_col = "variant_id"
  )

  expect_equal(nrow(result), 4L)
  expect_true(all(c("CHROM", "POS", "ID", "gene_id", "gene_name", "distance") %in% names(result)))
  expect_equal(result$distance[result$ID == "inside"], c(0L, 0L))
  expect_equal(result$distance[result$ID == "upstream"], 10L)
  expect_equal(result$distance[result$ID == "downstream"], 10L)
})

test_that("find_nearest_gene.GWASFormatter preserves unmatched variants", {
  gene <- human_genes %>%
    dplyr::filter(gene_biotype == "protein_coding", !is.na(gene_name)) %>%
    dplyr::slice(1)

  workdir <- tempfile("nearest-gene-")
  dir.create(workdir)

  old_wd <- setwd(workdir)
  on.exit(setwd(old_wd), add = TRUE)

  input_file <- tempfile(fileext = ".csv")
  utils::write.csv(
    data.frame(
      CHROM = c(gene$chrom, "chrUn"),
      POS = c(as.integer(gene$start), 1L),
      ALLELE0 = c("A", "A"),
      ALLELE1 = c("C", "C"),
      A1FREQ = c(0.1, 0.2),
      BETA = c(0.01, 0.02),
      SE = c(0.001, 0.002),
      LOG10P = c(8, 3)
    ),
    input_file,
    row.names = FALSE,
    quote = TRUE
  )

  gwas <- reformat_summary_statistics(input_file, use_cache = FALSE)
  on.exit(
    {
      if (DBI::dbIsValid(gwas$con)) {
        DBI::dbDisconnect(gwas$con, shutdown = TRUE)
      }
    },
    add = TRUE
  )

  annotated <- find_nearest_gene(gwas, threshold = 100000)
  result <- annotated$data %>%
    dplyr::select(CHROM, POS, ID, gene_name, distance) %>%
    dplyr::collect()

  expect_equal(nrow(result), 2L)
  expect_true(is.na(result$gene_name[result$CHROM == "chrUn"]))
  expect_true(is.na(result$distance[result$CHROM == "chrUn"]))
})