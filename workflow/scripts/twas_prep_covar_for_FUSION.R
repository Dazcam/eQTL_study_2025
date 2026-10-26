#--------------------------------------------------------------------------------------
#
#    Prepare covariate file for FUSION
#
#--------------------------------------------------------------------------------------
#
# Pipeline:  11TWAS | Rule: prepare_covar (one job per cell type)
#            Upstream:   05TENSORQTL split_covariates (final genotype and expression PC
#                        settings for each cell type)
#            Downstream: compute_weights
#
# Purpose:   - Transpose the tensorQTL covariate file (one row per covariate) to
#              PLINK covariate format (FID, IID, one column per covariate) for
#              FUSION.compute_weights.R --covar
#
# Inputs:    covar_in   tensorQTL covariate file: first column is the covariate
#                       name (unnamed in the header), then one column per sample
#
# Outputs:   covar_out  PLINK-format covariate file: FID, IID, one column per covariate
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

message("\n\nPreparing covariate file for FUSION weight calculation ...")

##  Load packages, functions and variables  -------------------------------------------
library(tidyverse)

# Load functions
source(file.path(snakemake@scriptdir, "twas_functions.R"))

# Input and output paths
covar_in <- snakemake@input[[1]]
covar_out <- snakemake@output[[1]]

# Print what each variable is set to
message("\nVariables")
message("============================")
tibble(
  variable = c("covar_in", "covar_out"),
  value    = c(covar_in, covar_out)) |>
  knitr::kable(format = "simple") |>
  cat(sep = "\n")
message("============================\n")

##  Read covariates  ------------------------------------------------------------------
# Name the first (unnamed) column 'covariate'
message("Reading covariate file ...")
covar <- read_tsv(covar_in,
                  col_names = c("covariate", readLines(covar_in, n = 1) %>% str_split("\t") %>% .[[1]] %>% .[-1]),
                  skip = 1)

##  Transpose and save  ---------------------------------------------------------------
message("Transposing covariate file ...")
covar_t <- to_plink_pheno(covar, id_col = "covariate")

message("Writing transposed covariate file ...")
write_tsv(covar_t, covar_out)

message("Done.")

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
