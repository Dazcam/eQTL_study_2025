#--------------------------------------------------------------------------------------
#
#    GWAS helper functions (shared across pipelines)
#
#--------------------------------------------------------------------------------------
#
# Functions: clean_gwas()  Drop non-rsIDs, duplicate rsIDs (keeps first) and SNPs on
#                          'alt' contigs
#
# Packages:  dplyr and stringr must be loaded by the calling script
#
# Used in:   smr_format_gwas.R (10SMR), ctwas_run.R (12CTWAS),
#            check_gwas_mismatched_snps.R (one-off check)
#
#--------------------------------------------------------------------------------------

## GWAS cleaning  -----------------------------------------------------------------------
# Removes non-rsID SNPs, duplicate rsIDs (keeps first) and SNPs on 'alt' contigs.
# Expects SNP and CHR columns. Returns the cleaned GWAS table.
clean_gwas <- function(gwas_tbl) {

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

  gwas_tbl
}

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
