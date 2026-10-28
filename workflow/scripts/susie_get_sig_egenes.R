#--------------------------------------------------------------------------------------
#
#    Extract significant eGenes for SuSiE fine-mapping
#
#--------------------------------------------------------------------------------------
#
# Pipeline:  08SUSIE | Rule: get_sig_eGenes (one job per cell type)
#            Upstream:   05TENSORQTL (tensorQTL permutation output)
#            Downstream: none at present (run_susie applies the same filter itself)
#
# Purpose:   - Read the tensorQTL permutation output for one cell type
#            - Keep genes with qval < 0.05
#            - Write one row per gene
#
# Inputs:    eqtl_in    tensorQTL permutation output (*_perm.cis_qtl.txt.gz)
#
# Outputs:   genes_out  TSV with a single phenotype_id column
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

message("\n\nExtracting table of sig. eGenes for SuSiE ...")

##  Load packages, functions and variables  -------------------------------------------
library(tidyverse)

# Input and output paths
eqtl_in   <- snakemake@input[[1]]
genes_out <- snakemake@output[[1]]
cell_type <- snakemake@wildcards[["cell_type"]]

# Print what each variable is set to
message("\nVariables")
message("============================")
tibble(
  variable = c("eqtl_in", "genes_out", "cell_type"),
  value    = c(eqtl_in, genes_out, cell_type)) |>
  knitr::kable(format = "simple") |>
  cat(sep = "\n")
message("============================\n")

##  Extract eGenes  ---------------------------------------------------------------------
message("Loading eQTL file ...")
eqtl_tbl <- read_tsv(eqtl_in, show_col_types = FALSE)

message("Extracting genes at FDR < 0.05 ...")
sig_genes_tbl <- eqtl_tbl |>
  dplyr::filter(qval < 0.05) |>
  distinct(phenotype_id)

if (nrow(sig_genes_tbl) == 0) stop("No genes at FDR < 0.05 for ", cell_type)

write_tsv(sig_genes_tbl, genes_out)
message(nrow(sig_genes_tbl), " genes written to file.")

message("Done.")

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
