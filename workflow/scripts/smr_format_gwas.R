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
message("GWAS loaded from: ", gwas_in)
message("Allele freq file loaded from: ", frq_in)

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
# non-rsIDs
non_rsids <- gwas_tbl %>%
  filter(!str_starts(SNP, "rs")) %>%
  distinct(SNP) %>%
  nrow()
if (non_rsids > 0) {
  message(str_glue("Found {non_rsids} non-rsID SNPs. Removing them ..."))
  gwas_tbl <- gwas_tbl %>% filter(str_starts(SNP, "rs"))
} else {
  message("No non-rsID SNPs found.")
}

# duplicate rsIDs
duplicate_rsids <- gwas_tbl %>%
  group_by(SNP) %>%
  filter(n() > 1) %>%
  distinct(SNP) %>%
  nrow()
if (duplicate_rsids > 0) {
  message(str_glue("Found {duplicate_rsids} duplicate rsIDs. Keeping first occurrence ..."))
  gwas_tbl <- gwas_tbl %>% distinct(SNP, .keep_all = TRUE)
} else {
  message("No duplicate rsIDs found.")
}

# Alt Chrs 
message("Checking GWAS chromosomes ...\n")
gwas_tbl |> 
  group_by(CHR) |> 
  count() |> 
  print(n = Inf)

alt_chroms <- gwas_tbl %>%
  filter(str_ends(CHR, "alt")) %>%
  distinct(CHR) %>%
  pull(CHR)

if (length(alt_chroms) > 0) {
  message("\nChromosomes ending with 'alt' found: ", paste(alt_chroms, collapse = ", "))
  message("\nCount of SNPs with 'alt' chromosomes: ", 
          nrow(gwas_tbl %>% filter(str_ends(CHR, "alt"))))
  message("Removing SNPs with 'alt' chromosomes ... ")
  gwas_tbl <- gwas_tbl %>% filter(!str_ends(CHR, "alt"))
} else {
  message("No chromosomes ending with 'alt' found.")
}

# Join GWAS with reference frequencies by CHR and SNP
message("\nAppending freq info to GWAS ...\n")
mrg_tbl <- gwas_tbl %>%
  left_join(frq_tbl, by = c("CHR", "SNP"), suffix = c("", ".ref"))

# New strategy: 
# - no allele filpping - 80% are flipped
# - prioritise GWAS allele order over reference despite 80% mismatch
# - But flip MAF if necessary
mrg_tbl <- mrg_tbl %>%
  mutate(
    allele_status = case_when(
      A1 == A1.ref & A2 == A2.ref ~ "aligned",
      A1 == A2.ref & A2 == A1.ref ~ "flipped",
      TRUE ~ "mismatched"
    ),
    freq = case_when(
      allele_status == "aligned" ~ MAF,
      allele_status == "flipped" ~ 1 - MAF,
      TRUE ~ NA_real_
    )
    # b remains unchanged
  )

# Report allele flips
message("\n=== GWAS frequency alignment summary ===")
mrg_tbl %>%
  count(allele_status) %>%
  mutate(percent = round(100 * n / sum(n), 1)) %>%
  arrange(desc(n)) %>%
  print(n = Inf)
message("========================================\n")

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
