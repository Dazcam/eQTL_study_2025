#--------------------------------------------------------------------------------------
#
#    Create SNP lookup table for snATAC-seq enrichment
#
#--------------------------------------------------------------------------------------

# 1. Load sig. eQTL an restrict to FDR < 0.05
# 2. Use BiomaRt to get hg38 SNP coordinates for 1

# Notes: Needs to run on head node BiomaRt requires net access

## Info  ------------------------------------------------------------------------------

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

message("\n\nObtain SNP lookup file for cell-specific sig. eQTL ...")

# -------------------------------------------------------------------------------------

library(tidyverse)
library(biomaRt)
library(GenomicRanges)

# Input and output paths
qtl_indep <- snakemake@input[["qtl_indep"]] 
out_file <- as.character(snakemake@output)
cell_type <- snakemake@wildcards[["cell_type"]]

# Check variable assignment
message("\nVariables")
cat("============================")
tibble(
  variable = c("qtl_indep", "out_file"),
  value    = c(qtl_indep, out_file)) |> 
  knitr::kable(format = "simple", align = "l") |>
  print()
message("\n============================\n")

message("\n--- Generating SNP lookup file ---")
message("\nProcessing: ", qtl_indep)

# Load significant eQTL SNPs
my_eqtl <- read_tsv(qtl_indep) %>%
  filter(rank == 2) %>%
  dplyr::select(variant_id)

if (nrow(my_eqtl) < 20) {
  message("Skipping ", cell_type, ": only ", nrow(my_eqtl), " secondary independent eQTLs (need ≥20)")
  
  # Create tbl for stub file
  empty_lookup <- tibble(
    refsnp_id = character(),
    chr_name = character(),
    chrom_start = integer(),
    chrom_end = integer()
  )
  
  message("Writing empty SNP lookup stub ...")
  write_rds(empty_lookup, out_file)
  message("Done.")

} else {
  
  message("Proceeding with ", nrow(my_eqtl), " secondary independent eQTLs")

  # Initialize Ensembl once
  temp_cache <- tempdir()
  Sys.setenv(BIOMART_CACHE = temp_cache)
  biomartCacheClear() # As running on head node need to use temp cache
  ensembl <- useEnsembl("ENSEMBL_MART_SNP", dataset = "hsapiens_snp")
  
  # SNP lookup
  message("Pulling out SNPs with BiomaRt ...")
  snp_lookup_hg38 <- getBM(
    attributes = c("refsnp_id", "chr_name", "chrom_start", "chrom_end"),
    filters = "snp_filter",
    values =  my_eqtl$variant_id,
    mart = ensembl,
    uniqueRows = TRUE
  )
  
  message("Filtering SNPs ...")
  snp_lookup_hg38_filt <- snp_lookup_hg38 %>%
    filter(chr_name %in% c(1:22, "X", "Y")) %>%
    mutate(chr_name = paste0("chr", chr_name)) %>%
    distinct()
  
  message(nrow(snp_lookup_hg38) - nrow(snp_lookup_hg38_filt), 
          ' eQTL removed due to non-standard chr annotation.')
  
  message("Number of SNPs retained: ", nrow(snp_lookup_hg38_filt))
  message("Writing SNP list ...")
  write_rds(snp_lookup_hg38_filt, out_file)
  message('Done.')

}

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------