#--------------------------------------------------------------------------------------
#
#    Extract sig. eGenes with cis window
#
#--------------------------------------------------------------------------------------

## Load libraries and variables -------------------------------------------------------
library(tidyverse)
library(biomaRt)

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
message('\n\nExtracting table of sig. eGenes with cis-window ...')
eqtl_file <- snakemake@input[[1]]
out_file <- snakemake@output[[1]]

# Load eQTL results from tensorQTL output
message("Loading eQTL file ...")
eqtl_data <- read_tsv(eqtl_file)

# Extract unique genes with q-value < 0.05
message("Extract genes @ FDR < 0.05 ...")
sig_genes <- eqtl_data %>%
  filter(qval < 0.05) %>%
  distinct(phenotype_id) %>%
  pull(phenotype_id)

# Connect to Ensembl (hg38) using biomaRt
message("Use biomaRt to get gene window cordinates for ", length(sig_genes), " sig. eGenes on hg38 ...")
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Fetch TSS information for significant genes
tss_info <- getBM(
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
  )

# Save the annotated cis-windows
message("Writing tsv file ...")
tss_info %>%
  dplyr::select(phenotype_id, group_id, gene_id, chromosome, phenotype_pos, strand) %>%
  write_tsv(out_file)

message("All done.")
#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
