#--------------------------------------------------------------------------------------
#
#    Create gene lists for SMR plotting
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
message('\n\nCreating the gene lists file for SMR plotting  ...')

## Load libraries and variables -------------------------------------------------------
library(tidyverse)
library(biomaRt)

gene_list_file <- snakemake@input[[1]]
ensembl_out_file <- snakemake@output[['ensembl']]
symbol_out_file <- snakemake@output[['symbol']]
cell_type <- snakemake@wildcards[["cell_type"]]

# Check variable assignment
message("\nVariables")
cat("============================")
tibble(
  variable = c("gene_list_file", "out_file", "cell_type"),
  value    = c(gene_list_file, out_file, cell_type)
) |> 
  knitr::kable(format = "simple", align = "l") |>
  print()
message("\n============================\n")

# Load gex results from cell-specific pseudobulk file
message("Loading gene expression file ...\n")
gene_list <- read_tsv(gene_list_file) |>
  pull(Gene)

# Connect to Ensembl (hg38) using biomaRt
message("Use biomaRt to get start / stop coords and strand and for ", length(gene_list), " genes ...")
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Fetch start / stop coords, strand and symbols for gene list
biomart_tbl <- getBM(
  attributes = c("chromosome_name", "start_position", "end_position", "ensembl_gene_id", "external_gene_name", "strand"),
  filters = "ensembl_gene_id",
  values = gene_list,
  mart = mart) |>
  as_tibble() |>
  mutate(strand = ifelse(strand == 1, "+", "-"))
  distinct()

message(nrow(biomart_tbl), " distinct genes remain after running BiomaRt ...")

# Write tsvs
message("Writing ensembl id tsv file ...")
ensembl_tbl <- biomart_tbl |>
  select(chromosome_name, start_position, end_position, ensembl_gene_id, strand) |>
  write_tsv(ensembl_out_file, col_names = FALSE)

message("Writing symbol id tsv file ...")
symbol_tbl <- biomart_tbl |>
  select(chromosome_name, start_position, end_position, external_gene_name, strand) |>
  write_tsv(symbol_out_file, col_names = FALSE)


message("All done.")
#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
