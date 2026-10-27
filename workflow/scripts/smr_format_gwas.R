#--------------------------------------------------------------------------------------
#
#    Prepare GWAS for SMR
#
#--------------------------------------------------------------------------------------
#
# Pipeline:  10SMR | Rule: format_gwas (one job per GWAS)
#            Upstream:   07PREP-GWAS ({gwas}_hg38.tsv); cat_refs (1000G EUR hg38 .frq)
#            Downstream: smr
#
# Purpose:   Convert one GWAS to SMR .ma format (SNP A1 A2 freq b se p N):
#            - Take b from BETA, or log(OR) when only OR is given
#            - Drop non-rsIDs, duplicate rsIDs and alt contigs (clean_gwas())
#            - Set freq to the 1000G frequency of the GWAS A1; alleles and b are
#              unchanged (align_freq_to_ref())
#            - Drop SNPs with no matching 1000G frequency
#
# Inputs:    gwas         {gwas}_hg38.tsv (SNP, CHR, BP, A1, A2, BETA or OR, SE, PVAL, N)
#            frq          1000G.EUR.hg38.frq
#
# Outputs:   ma           {gwas}_hg38.ma
#
# Notes:     SMR checks alleles across the GWAS, eQTL and LD reference itself, so the
#            GWAS allele order is kept. 
#
#--------------------------------------------------------------------------------------

## Set up logging for Snakemake
if (exists("snakemake")) {
  log_smk <- function() {
    if (exists("snakemake") && length(snakemake@log) != 0) {
      log <- file(snakemake@log[[1]], open = "wt")
      sink(log, append = TRUE)
      sink(log, append = TRUE, type = "message")
    }
  }
  log_smk()
}

message("\n\nMunge GWAS for SMR input ...")

##  Load packages, functions and variables  -------------------------------------------
library(tidyverse)

# Load functions
source(file.path(snakemake@scriptdir, "gwas_functions.R"))
source(file.path(snakemake@scriptdir, "smr_functions.R"))

# Input and output paths
gwas_in <- snakemake@input[["gwas"]]
frq_in <- snakemake@input[["frq"]]
gwas_out <- snakemake@output[[1]]

# Print what each variable is set to
message("\nVariables")
message("============================")
tibble(
  variable = c("gwas_in", "frq_in", "gwas_out"),
  value    = c(gwas_in, frq_in, gwas_out)) |>
  knitr::kable(format = "simple") |>
  cat(sep = "\n")
message("============================\n")

## Load data  -------------------------------------------------------------------------
message("Loading data ...\n")

gwas_tbl <- read_tsv(gwas_in, guess_max = 100000, col_types = cols(.default = "c"))

# Determine effect column
if ("BETA" %in% names(gwas_tbl)) {
  gwas_tbl <- gwas_tbl %>%
    mutate(
      BP   = as.integer(BP),
      PVAL = as.numeric(PVAL),
      SE   = as.numeric(SE),
      N    = as.integer(N),
      BETA = as.numeric(BETA)
    )
} else if ("OR" %in% names(gwas_tbl)) {
  gwas_tbl <- gwas_tbl %>%
    mutate(
      BP   = as.integer(BP),
      PVAL = as.numeric(PVAL),
      SE   = as.numeric(SE),
      N    = as.integer(N),
      OR   = as.numeric(OR),
      BETA = log(OR)
    )
} else {
  stop("Neither BETA nor OR column found in GWAS file: ", gwas_in)
}

gwas_tbl <- gwas_tbl %>%
  rename(b = BETA, se = SE, p = PVAL) %>%
  select(SNP, CHR, BP, A1, A2, b, se, p, N)

frq_tbl <- read_table(frq_in, col_types = cols(.default = "c", MAF = "d", NCHROBS = "i")) 

## GWAS checks  -----------------------------------------------------------------------
gwas_tbl <- clean_gwas(gwas_tbl)

# Join GWAS with reference frequencies by CHR and SNP
# New strategy: 
# - no allele filpping - 80% are flipped
# - prioritise GWAS allele order over reference despite 80% mismatch
# - But flip MAF if necessary
message("\nAppending freq info to GWAS ...\n")
mrg_tbl <- align_freq_to_ref(gwas_tbl, frq_tbl, by = c("CHR", "SNP"), freq_col = "freq",
                             label = "GWAS")

## Write outputs  ---------------------------------------------------------------------
# Filter out SNPs with missing freq
message("Filtering SNPs with missing or mismatched frequencies ...\n")
result <- mrg_tbl %>%
  filter(!is.na(freq)) %>%
  select(SNP, A1, A2, freq, b, se, p, N)

# Write output
message(str_glue("Writing GWAS in .ma format to {gwas_out} ...\n"))
write_tsv(result, gwas_out, na = "NA")

message(str_glue("Formatted GWAS data written to {gwas_out}"))
message("Done.\n")

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
