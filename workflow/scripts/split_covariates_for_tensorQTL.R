#--------------------------------------------------------------------------------------
#
#    Split covariates for tensorQTL tensorQTL
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

message("\n\nPulling out genotype and expression PC covariate combinations ...")

##  Load Packages, functions and variables  -------------------------------------------
# Install and load required libraries
library(tidyverse)

cov_input <- as.character(snakemake@input)
cov_output <- as.character(snakemake@output)
cell_type <- snakemake@wildcards[["cell_type"]]
geno_pc <- snakemake@wildcards[["geno_pc"]]
exp_pc <- snakemake@wildcards[["exp_pc"]]

# Make a tibble showing what each variable is set to
message("============================")
tibble(
  variable = c("cov_input", "cov_output", "cell_type", "geno_pc", "exp_pc"),
  value    = c(cov_input, cov_output, cell_type, geno_pc, exp_pc)) |> 
  knitr::kable(format = "simple", align = "l") |>
  print()
message("============================\n")
    
  
# Pull out genotype and expression PC covarite combinations  
cov_tbl <- read_tsv(cov_input, show_col_types = FALSE) |>
  rename(covariate = `...1`)
message("Input covariate matrix for ", cell_type, ":\n")
print(cov_tbl)

cov_final_tbl <- cov_tbl |>
  filter(covariate %in% c("PCW", "sex_code", paste0("genPC", 1:geno_pc), paste0("expPC", 1:exp_pc))) |>
  mutate(covariate = if_else(covariate == "sex_code", "Sex", covariate))

message("\nOutput covariate matrix:\n")
print(cov_final_tbl)

# Create matrix with covariates as rows and samples as columns
cov_matrix <- as.matrix(cov_final_tbl[, -1])
rownames(cov_matrix) <- cov_final_tbl$covariate
colnames(cov_matrix) <- colnames(cov_final_tbl)[-1]
message("\nCovariate matrix (first 5 rows and columns):\n")
print(cov_matrix[1:5, 1:5])

message("\nWriting covariate matrix ...\n")
write.table(cov_matrix, 
            file = cov_output, 
            sep = "\t", 
            quote = FALSE, 
            col.names = NA,    # Write sample IDs as column names
            row.names = TRUE,    # Write covariate names as row names
            append = FALSE)

message("Done.")

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------