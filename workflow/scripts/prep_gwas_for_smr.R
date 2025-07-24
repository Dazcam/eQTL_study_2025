#--------------------------------------------------------------------------------------
#
#    Prep GWAS for SMR
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

# Load packages
library(tidyverse)

# Input and output paths
gwas_in <- snakemake@input[["gwas"]]
frq_in <- snakemake@input[["frq"]]
gwas_out <- snakemake@output[[1]]

# Read in data
message("\nLoading data ...\n")
gwas_tbl <- read_tsv("../results/05SLDSR/gwas/scz_hg38.tsv", 
                     col_types = cols(.default = "c", BP = "i", PVAL = "d", 
                                      BETA = "d", SE = "d", Z = "d", N = "i")) |>
  rename(b = BETA, se = SE, p = PVAL) |>
  select(SNP, CHR, BP, A1, A2, b, se, p, N)

frq_tbl <- read_table('../results/07SMR/smr_input/allele_frq_hg38_ref.txt', 
                      col_types = cols(.default = "c", MAF = "d", NCHROBS = "i")) 

# Join GWAS with reference frequencies by CHR and SNP
message("Appending freq info to GWAS ...\n")
mrg_tbl <- gwas_tbl %>%
  left_join(frq_tbl, by = c("CHR", "SNP"), suffix = c("", ".ref"))

# Check initial mismatches to determine strand flip handling
initial_mismatch_tbl <- mrg_tbl %>%
  filter(!(A1 == A1.ref & A2 == A2.ref | A1 == A2.ref & A2 == A1.ref) & !is.na(A1.ref))
use_strand_flip <- nrow(initial_mismatch_tbl) > 100
message(str_glue("Initial mismatches: {nrow(initial_mismatch_tbl)}. Strand flip handling: {ifelse(use_strand_flip, 'enabled', 'disabled')}.\n"))

# Assign freq and status based on allele alignment
message("Assigning frequencies and checking for allele mismatches ...\n")
mrg_tbl <- mrg_tbl %>%
  mutate(
    freq = case_when(
      A1 == A1.ref & A2 == A2.ref ~ MAF,           # Alleles align directly
      A1 == A2.ref & A2 == A1.ref ~ 1 - MAF,       # Alleles are flipped
      # Strand flips for A/T or C/G SNPs (only if mismatches > 100)
      use_strand_flip & (
        (A1 == "A" & A2 == "T" & A1.ref == "T" & A2.ref == "A") |
          (A1 == "T" & A2 == "A" & A1.ref == "A" & A2.ref == "T") |
          (A1 == "C" & A2 == "G" & A1.ref == "G" & A2.ref == "C") |
          (A1 == "G" & A2 == "C" & A1.ref == "C" & A2.ref == "G")
      ) ~ 1 - MAF,                                # Treat as flipped
      TRUE ~ NA_real_                              # Mismatch or missing
    ),
    allele_status = case_when(
      A1 == A1.ref & A2 == A2.ref ~ "aligned",
      A1 == A2.ref & A2 == A1.ref ~ "flipped",
      use_strand_flip & (
        (A1 == "A" & A2 == "T" & A1.ref == "T" & A2.ref == "A") |
          (A1 == "T" & A2 == "A" & A1.ref == "A" & A2.ref == "T") |
          (A1 == "C" & A2 == "G" & A1.ref == "G" & A2.ref == "C") |
          (A1 == "G" & A2 == "C" & A1.ref == "C" & A2.ref == "G")
      ) ~ "strand_flipped",
      is.na(A1.ref) ~ "missing_reference",
      TRUE ~ "mismatched"
    )
  )

# Log allele alignment statistics
allele_counts <- mrg_tbl %>%
  count(allele_status) %>%
  mutate(percentage = n / sum(n) * 100)
message("Allele alignment summary:\n")
walk2(allele_counts$allele_status, allele_counts$percentage, 
      ~ message(str_glue("{.x}: {.y}% (n = {allele_counts$n[allele_counts$allele_status == .x]})")))

# Check for true allele mismatches (excluding missing reference)
message("Checking for true allele mismatches ...\n")
mismatch_tbl <- mrg_tbl %>%
  filter(is.na(freq) & !is.na(A1.ref))
if (nrow(mismatch_tbl) > 0) {
  warning(str_glue("\n{nrow(mismatch_tbl)} SNPs have mismatched alleles (e.g., non-strand-flip incompatibilities) and will be excluded"))
}

# Check for missing reference SNPs
missing_ref_tbl <- mrg_tbl %>%
  filter(is.na(freq) & is.na(A1.ref))
if (nrow(missing_ref_tbl) > 0) {
  warning(str_glue("\n{nrow(missing_ref_tbl)} SNPs are missing from the reference and will be excluded"))
}

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