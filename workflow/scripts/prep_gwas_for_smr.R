#--------------------------------------------------------------------------------------
#
#    Prep GWAS for SMR
#
#--------------------------------------------------------------------------------------

## Set up logging for smk  ------------------------------------------------------------
if (exists("snakemake")) {
  log_smk <- function() {
    if (exists("snakemake") & length(snakemake@log) != 0) {
      log <- file(snakemake@log[1][[1]], open = "wt")
      sink(log, append = TRUE)
      sink(log, append = TRUE, type = "message")
    }
  }
}
log_smk()
message('\n\nMunge eQTL, snp and gene level data for SMR input ...')

# Load packages
library(tidyverse)

# Input and output paths
gwas_in <- snakemake@input[['gwas']]
frq_in <- snakemake@input[['frq']]
gwas_out <- snakemake@output[['query']]

# Read in data
gwas_tbl <- read_tsv(gwas_in) |>
  rename(b = BETA, se = SE, p = PVAL) |>
  select(SNP, CHR, BP, A1, A2, b, se, p, N)

frq_tbl <- read_tsv(frq_in)

# Join GWAS with reference frequencies by CHR and SNP
mrg_tbl <- gwas_tbl %>%
  left_join(frq_tbl, by = c("CHR", "SNP"), suffix = c("", ".ref"))

# Assign freq based on allele alignment
mrg_tbl <- mrg_tbl %>%
  mutate(
    freq = case_when(
      A1 == A1.ref & A2 == A2.ref ~ MAF,
      A1 == A2.ref & A2 == A1.ref ~ 1 - MAF,
      TRUE ~ NA_real_
    )
  )

# Check for mismatches or missing frequencies
mismatches <- mrg_tbl %>%
  filter(is.na(freq) & !is.na(A1.ref))
if (nrow(mismatches) > 0) {
  warning(str_glue("{nrow(mismatches)} SNPs have mismatched alleles or missing frequencies and will be excluded"))
}

# Filter out SNPs with missing freq
result <- mrg_tbl %>%
  filter(!is.na(freq)) %>%
  select(SNP, A1, A2, freq, b, se, p, N)

# Write output
write_tsv(result, output_file, na = "NA")

message(str_glue("Formatted GWAS data written to {output_file}"))

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------