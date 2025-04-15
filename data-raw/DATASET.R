## code to prepare `DATASET` dataset goes here

download.file("https://ftp.ensembl.org/pub/release-113/gtf/homo_sapiens/Homo_sapiens.GRCh38.113.gtf.gz", "homo_sapiens.GRCh38.113.gtf.gz")
annot = rtracklayer::import("homo_sapiens.GRCh38.113.gtf.gz")

human_genes = as.data.frame(annot) %>%
  dplyr::filter(type == "gene") %>%
  dplyr::select(
    chrom = seqnames,
    start,
    end,
    gene_name,
    gene_id,
    gene_biotype
  ) %>%
  dplyr::mutate(chrom = glue("chr{chrom}")) %>%
  dplyr::distinct(.)

usethis::use_data(human_genes, overwrite = TRUE)

# download centromere data from UCSC genome browser
ideogram = arrow::read_tsv_arrow("centromere.txt") %>%
  dplyr::select(
    chrom = `#chrom`,
    start = chromStart,
    end = chromEnd,
    name,
    stain = gieStain
  ) %>%
  dplyr::distinct(.)

usethis::use_data(ideogram, overwrite = TRUE)
