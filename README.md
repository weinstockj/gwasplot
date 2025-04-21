# Overview

**gwasplot** provides fast and efficient tools for visualizing and annotating GWAS 
summary statistics. It offers functions to reformat data, create plots, and annotate 
top hits with genomic features. 

It is designed to work with output from GWAS performed on WGS data with many rare variants, using [duckdb](https://duckdb.org/) to manipulate data. 

For more information, visit the [GitHub repository](https://github.com/weinstockj/gwasplot) or the [project website](https://weinstockj.github.io/gwasplot/).

# Installation

Install the development version of **gwasplot** from GitHub:

```r
# Install remotes if not already installed
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}

# Install gwasplot from GitHub
remotes::install_github("weinstockj/gwasplot")
```


## Input formats

`gwasplot` expects output from either regenie or saige in csv or parquet format (parquet is generally recommended). 

Regenie output looks like this:

```r
dplyr::glimpse(df)
Rows: 5
Columns: 12
$ CHROM     <chr> "chr3", "chr3", "chr3", "chr3", "chr3"
$ POS       <int> 90000011, 90000261, 90000274, 90000324, 90000366
$ ID        <chr> "chr3-90000011-T-C", "chr3-90000261-GATT-G", "chr3-900002…"
$ ALLELE0   <chr> "T", "GATT", "T", "C", "T"
$ ALLELE1   <chr> "C", "G", "A", "A", "G"
$ A1FREQ    <dbl> 1.37775e-04, 2.91090e-04, 2.38259e-04, 1.62637e-04, 4.66158e…
$ N         <int> 482670, 482669, 482669, 482669, 482669
$ BETA      <dbl> -0.00157451, -0.06463220, -0.06452010, 0.02796820, -0.021063…
$ SE        <dbl> 0.0730887, 0.0479369, 0.0540183, 0.0646778, 0.1213880
$ CHISQ     <dbl> 0.000464077, 1.817850000, 1.426620000, 0.186990000, 0.030108…
$ LOG10P    <dbl> 0.00752913, 0.75063100, 0.63391900, 0.17689500, 0.06437000
```

gwastools will read in data from either regenie or saige and standarize the column names,
and then create a table called `summary_stats` in a `duckdb` database. 

# Usage
Below are some sample usage examples:

Reformat GWAS Summary Statistics
Use the reformat_summary_statistics function to read and reformat GWAS summary statistics from a parquet or CSV file:

```r
library(gwasplot)

# Reformat summary statistics from a file
gwas_stats <- reformat_summary_statistics("path/to/your/file.parquet")
print(gwas_stats )

GWAS object
File path: ../concatenated_results.parquet
Detected format: regenie
Data names: CHROM, POS, ID, ALLELE0, ALLELE1, A1FREQ, N, BETA, SE, CHISQ, LOG10P, phenotype
Data preview:
# A tibble: 5 × 9
  CHROM      POS REF   ALT      AF_ALT     BETA  LOG10P PVALUE ID
  <chr>    <dbl> <chr> <chr>     <dbl>    <dbl>   <dbl>  <dbl> <chr>
1 chr3  90000011 T     C     0.000138  -0.00157 0.00753  0.983 chr3_90000011_T_C
2 chr3  90000261 GATT  G     0.000291  -0.0646  0.751    0.178 chr3_90000261_GA…
3 chr3  90000274 T     A     0.000238  -0.0645  0.634    0.232 chr3_90000274_T_A
4 chr3  90000324 C     A     0.000163   0.0280  0.177    0.665 chr3_90000324_C_A
5 chr3  90000366 T     G     0.0000466 -0.0211  0.0644   0.862 chr3_90000366_T_G
```

Then create a manhattan plot:

```r
manhattan(gwas = gwas_stats, output_prefix = "my_manhattan")
```