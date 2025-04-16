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

chip_genes = c("ASXL1", "ASXL2", "BCOR", "BCORL1", "BRAF", "BRCC3", "CBL", "CBLB", 
                 "CEBPA", "CREBBP", "CSF1R", "CSF3R", "CTCF", "CUX1", "DNMT3A", "EED", 
                 "EP300", "ETNK1", "ETV6", "EZH2", "FLT3", "GATA1", "GATA2", "GATA3", 
                 "GNA13", "GNAS", "GNB1", "IDH1", "IDH2", "IKZF1", "IKZF2", "IKZF3", 
                 "JAK1", "JAK2", "JAK3", "KDM6A", "KIT", "KRAS", "LUC7L2", "KMT2A", 
                 "KMT2D", "MPL", "NF1", "NPM1", "NRAS", "PDS5B", "PDSS2", "PHF6", 
                 "PHIP", "PPM1D", "PRPF40B", "PRPF8", "PTEN", "PTPN11", "RAD21", 
                 "RUNX1", "SETBP1", "SETD2", "SETDB1", "SF1", "SF3A1", "SF3B1", 
                 "SRSF2", "SMC1A", "SMC3", "STAG1", "STAG2", "SUZ12", "TET2", "TP53", 
                 "U2AF1", "U2AF2", "WT1", "ZRSR2", "ZBTB33", "YLPM1", "SRCAP", "ZNF318")

usethis::use_data(chip_genes, overwrite = TRUE)