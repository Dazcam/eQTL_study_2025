#--------------------------------------------------------------------------------------
#
#    Prep covariate file for FUSION 
#
#--------------------------------------------------------------------------------------

## Set up logging for smk  ------------------------------------------------------------
if (exists("snakemake")) {
  log_smk <- function() {
    if (exists("snakemake") & length(snakemake@log) != 0) {
      log <- file(snakemake@log[1][[1]], open = "wt")
      sink(log, append = TRUE)
      sink(log, append = TRUE, type = "message")
    }
  }
}
log_smk()
message('\n\nMunging covariate file format for FUSION weight calculation ...')

# Load packages
library(tidyverse)

# Input and output paths
covar_in <- snakemake@input[[1]]
covar_out <- snakemake@output[[1]]

# Read the covariate file
# Specify the first column header as 'covariate' and treat it as tab-separated
message('Reading covariate file ...')
covar <- read_tsv(covar_in,
                  col_names = c("covariate", readLines(covar_in, n = 1) %>% str_split("\t") %>% .[[1]] %>% .[-1]),
                  skip = 1)

# Transpose the data
# Pivot the data so sample IDs become rows and covariates become columns
message('Transposing covariate file ...')
covar_t <- covar %>%
  pivot_longer(cols = -covariate, names_to = "IID", values_to = "value") %>%
  pivot_wider(names_from = covariate, values_from = value) %>%
  mutate(FID = 0) %>%
  select(FID, IID, everything())

# Write the corrected file
message('Writing transposed covariate file ...')
write_tsv(covar_transposed, covar_out)

message('Done.')
#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------

