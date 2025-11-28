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

# # Assign freq and status based on allele alignment
# message("Assigning frequencies and adjusting alleles ...\n")
# mrg_tbl <- mrg_tbl %>%
#   mutate(
#     A1_temp = case_when(
#       A1 == A1.ref & A2 == A2.ref ~ A1,           # Alleles align directly
#       A1 == A2.ref & A2 == A1.ref ~ A1.ref,       # Flipped, use reference alleles
#       TRUE ~ NA_character_                        # Mismatch or missing
#     ),
#     A2_temp = case_when(
#       A1 == A1.ref & A2 == A2.ref ~ A2,           # Alleles align directly
#       A1 == A2.ref & A2 == A1.ref ~ A2.ref,       # Flipped, use reference alleles
#       TRUE ~ NA_character_                        # Mismatch or missing
#     ),
#     freq = case_when(
#       !is.na(A1.ref) & (A1 == A1.ref & A2 == A2.ref | A1 == A2.ref & A2 == A1.ref) ~ MAF,  # Use MAF for aligned or flipped
#       TRUE ~ NA_real_                             # Mismatch or missing
#     ),
#     b = case_when(
#       A1 == A1.ref & A2 == A2.ref ~ b,            # Same alleles
#       A1 == A2.ref & A2 == A1.ref ~ -b,           # Flipped, reverse effect
#       TRUE ~ NA_real_                             # Mismatch or missing
#     ),
#     allele_status = case_when(
#       A1 == A1.ref & A2 == A2.ref ~ "aligned",
#       A1 == A2.ref & A2 == A1.ref ~ "flipped",
#       is.na(A1.ref) ~ "missing_reference",
#       TRUE ~ "mismatched"
#     ),
#     A1 = A1_temp,                                 # Rename to maintain naming convention
#     A2 = A2_temp
#   ) %>%
#   select(-A1_temp, -A2_temp)
# 
# # Log allele alignment statistics
# message("Allele alignment summary:\n")
# allele_counts <- mrg_tbl %>%
#   count(allele_status) %>%
#   mutate(percentage = n / nrow(gwas_tbl) * 100)
# walk2(allele_counts$allele_status, allele_counts$percentage,
#       ~ message(str_glue("{.x}: {.y}% (n = {allele_counts$n[allele_counts$allele_status == .x]})")))
# 
# # Check for true allele mismatches (excluding missing reference)
# message("\nChecking for true allele mismatches ...\n")
# mismatch_tbl <- mrg_tbl %>%
#   filter(is.na(freq) & !is.na(A1.ref))
# if (nrow(mismatch_tbl) > 0) {
#   warning(str_glue("\n{nrow(mismatch_tbl)} SNPs have mismatched alleles (e.g., non-strand-flip incompatibilities) and will be excluded"))
# }
# 
# # Check for missing reference SNPs
# missing_ref_tbl <- mrg_tbl %>%
#   filter(is.na(freq) & is.na(A1.ref))
# if (nrow(missing_ref_tbl) > 0) {
#   warning(str_glue("\n{nrow(missing_ref_tbl)} SNPs are missing from the reference and will be excluded"))
# }

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
