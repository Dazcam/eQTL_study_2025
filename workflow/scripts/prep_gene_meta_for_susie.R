#--------------------------------------------------------------------------------------
#
#    Create gene metadata file for run_susie.R
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
message('\n\nCreating the gene metadata file for all genes in GeX bed file  ...')

## Load libraries and variables -------------------------------------------------------
library(tidyverse)
library(biomaRt)

gex_file <- snakemake@input[[1]]
out_file <- snakemake@output[[1]]
cell_type <- snakemake@wildcards[["cell_type"]]

# Check variable assignment
message("\nVariables")
cat("============================")
tibble(
  variable = c("gex_file", "out_file", "cell_type"),
  value    = c(gex_file, out_file, cell_type)
) |> 
  knitr::kable(format = "simple", align = "l") |>
  print()
message("\n============================\n")

# Load gex results from cell-specific pseudobulk file
message("Loading gene expression file ...\n")
gex_data <- read_tsv(gex_file)

message("\nPulling out Ensembl IDs ...")
sig_genes <- gex_data %>%
  pull(TargetID)

# Connect to Ensembl (hg38) using biomaRt
message("Use biomaRt to get TSS and for ", length(sig_genes), " genes ...")
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Fetch TSS information for significant genes
tss_tbl <- getBM(
  attributes = c("ensembl_gene_id", "chromosome_name", "start_position", "end_position", "strand"),
  filters = "ensembl_gene_id",
  values = sig_genes,
  mart = mart
) %>%
  as_tibble() %>%
  dplyr::mutate(
    phenotype_id = ensembl_gene_id,
    group_id = ensembl_gene_id,
    gene_id = ensembl_gene_id,
    chromosome = chromosome_name,
    phenotype_pos = if_else(strand == 1, start_position, end_position),
    strand = strand
  ) %>%
  distinct()
message(nrow(tss_tbl), " distinct genes remain after running BiomaRt ...")

# Save the annotated cis-windows
message("Writing tsv file ...")
tss_tbl %>%
  dplyr::select(phenotype_id, group_id, gene_id, chromosome, phenotype_pos, strand) %>%
  write_tsv(out_file)

message("All done.")
#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
