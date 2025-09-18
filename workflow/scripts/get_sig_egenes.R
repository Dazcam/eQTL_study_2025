#--------------------------------------------------------------------------------------
#
#    Extract sig. eGenes - Pass as a list to run_susie.R for each cell type
#
#--------------------------------------------------------------------------------------

## Load libraries and variables -------------------------------------------------------
library(tidyverse)

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
message('\n\nExtracting table of sig. eGenes to use as input for run_susie.R ...')
eqtl_file <- snakemake@input[[1]]
out_file <- snakemake@output[[1]]
cell_type <- snakemake@wildcards[["cell_type"]]

# Check variable assignment
message("\nVariables")
cat("============================")
tibble(
  variable = c("eqtl_file", "out_file", "cell_type"),
  value    = c(eqtl_file, out_file, cell_type)
) |> 
  knitr::kable(format = "simple", align = "l") |>
  print()
message("\n============================\n")

# Load eQTL results from tensorQTL output
message("Loading eQTL file ...")
eqtl_data <- read_tsv(eqtl_file)

# Extract unique genes with q-value < 0.05
message("Extract genes @ FDR < 0.05 ...")
sig_genes_tbl <- eqtl_data %>%
  filter(qval < 0.05) %>%
  distinct(phenotype_id) %>%
  select(phenotype_id) %>%
  write_tsv(out_file)
message(nrow(sig_genes_tbl), " genes written to file ...")
message("All done.")

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
