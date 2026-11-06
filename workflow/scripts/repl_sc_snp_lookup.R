#--------------------------------------------------------------------------------------
#
#    Create hg38 SNP lookup table for sig. eQTL
#
#--------------------------------------------------------------------------------------
#
# Pipeline:  07REPLICATION-SC | Rule: snp_lookup (one job per cell type; localrule)
#            Upstream:   05TENSORQTL tensorqtl_perm
#            Downstream: atac_enrich
#
# Purpose:   - Load the cell type's sig. eQTL (qval < 0.05)
#            - Query Ensembl (biomaRt) for the hg38 coordinates of each lead SNP
#            - Keep SNPs on standard chromosomes (1-22, X, Y)
#
# Inputs:    qtl_perm   tensorQTL permutation output for the cell type
#
# Outputs:   out_file   rds: refsnp_id, chr_name, chrom_start, chrom_end
#
# Notes:     Runs on the head node (localrule) because biomaRt needs network access;
#            a temporary biomaRt cache is used for the same reason.
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

message("\n\nObtain SNP lookup file for cell-specific sig. eQTL ...")

##  Load packages, functions and variables  -------------------------------------------
library(tidyverse)
library(biomaRt)
library(GenomicRanges)

# Input and output paths
qtl_perm  <- snakemake@input[["qtl_perm"]]
out_file  <- as.character(snakemake@output)
cell_type <- snakemake@wildcards[["cell_type"]]

# Print what each variable is set to
message("\nVariables")
message("============================")
tibble(
  variable = c("qtl_perm", "out_file", "cell_type"),
  value    = c(qtl_perm, out_file, cell_type)) |>
  knitr::kable(format = "simple") |>
  cat(sep = "\n")
message("============================\n")

message("\n--- Generating SNP lookup file ---")
message("\nProcessing: ", qtl_perm)

##  Load sig. eQTL  ---------------------------------------------------------------------
my_eqtl <- read_tsv(qtl_perm) %>%
  filter(qval < 0.05) %>%
  dplyr::select(variant_id)

if (nrow(my_eqtl) == 0) {
  message("No significant eQTLs found for ", cell_type)
  quit(save = "no", status = 0)
}

##  Look up hg38 coordinates  -----------------------------------------------------------
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

##  Write output  ----------------------------------------------------------------------
message("Number of SNPs retained: ", nrow(snp_lookup_hg38_filt))
message("Writing SNP list ...")
write_rds(snp_lookup_hg38_filt, out_file)
message('Done.')

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
