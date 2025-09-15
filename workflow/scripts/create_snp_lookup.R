#--------------------------------------------------------------------------------------
#
#    Enrichment of eQTL in ATAC-seq peaks
#
#--------------------------------------------------------------------------------------

# Need to generate snp lookup tables for Ziffra enrichments in a seperate rule
# Need to run on head node net access is required

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
qtl_perm <- snakemake@input[["qtl_perm"]]
out_file <- snakemake@output[["out_file"]]

message("\n--- Generating SNP lookup file ---")
message("\nProcessing: ", qtl_perm)

# Load significant eQTL SNPs
my_eqtl <- read_tsv(qtl_perm) %>%
  filter(qval < 0.05) %>%
  dplyr::select(variant_id)

if (nrow(my_eqtl) == 0) {
  message("No significant eQTLs found for ", cell_type)
  quit(save = "no", status = 0)
}

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

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------