#--------------------------------------------------------------------------------------
#
#    Divide covriates for tensorQTL tensorQTL
#
#--------------------------------------------------------------------------------------

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

message("\n\nPulling out genotype and expression PC covarite combinations ...")

##  Load Packages, functions and variables  -------------------------------------------
# Install and load required libraries
library(tidyverse)

cov_input <- snakemake@input
cov_output <- snakemake@output
cell_type <- snakemake@wildcards[["cell_type"]]
geno_pc <- snakemake@wildcards[["geno_pc"]]
exp_pc <- snakemake@wildcards[["exp_pc"]]
    
  
# Pull out genotype and expression PC covarite combinations  
cov_tbl <- read_tsv(cov_input)
message("Input covariate matrix for ", cell_type, ":\n")
cov_tbl

cov_final_tbl <- cov_tbl |>
  dplyr::select(
    id, PCW, Sex = sex_code,
    all_of(paste0("genPC", 1:geno_pc)),
    all_of(paste0("expPC", 1:exp_pc))
  ) 

message("\nOutput covariate matrix:\n")
cov_final_tbl

message("\nWriting covarite matrix ...\n")
write.table(cov_final_tbl, 
            file = paste0(out_dir, cell_type, "_genPC_", geno_pc, "_expPC_", exp_pc, "_covariates.txt"), 
            sep = "\t", 
            quote = FALSE, 
            col.names = NA,    # Write sample IDs as column names
            row.names = TRUE,    # Write covariate names as row names
            append = FALSE)

message("Done.")

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------