#--------------------------------------------------------------------------------------
#
#    Create gene metadata file for SuSiE fine-mapping
#
#--------------------------------------------------------------------------------------
#
# Pipeline:  08SUSIE | Rule: prep_susie_gene_meta (one job per cell type; localrule,
#            needs internet access for biomaRt)
#            Upstream:   05TENSORQTL (prep_input expression bed)
#            Downstream: run_susie, susie_report
#
# Purpose:   - Read all genes in the cell type's expression bed
#            - Query Ensembl (biomaRt) for chromosome, gene start/end and strand
#            - Set phenotype_pos to the TSS (start on + strand, end on - strand)
#
# Inputs:    gex_in    tensorQTL expression bed (TargetID column holds Ensembl IDs)
#
# Outputs:   meta_out  TSV: phenotype_id, group_id, gene_id, chromosome,
#                      phenotype_pos, strand
#
# Notes:     useMart("ensembl") queries the current Ensembl release, so a rerun can
#            give different positions. Genes Ensembl does not return are dropped.
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

message("\n\nCreating the gene metadata file for all genes in the expression bed ...")

##  Load packages, functions and variables  -------------------------------------------
library(tidyverse)
library(biomaRt)

# Input and output paths
gex_in    <- snakemake@input[[1]]
meta_out  <- snakemake@output[[1]]
cell_type <- snakemake@wildcards[["cell_type"]]

# Print what each variable is set to
message("\nVariables")
message("============================")
tibble(
  variable = c("gex_in", "meta_out", "cell_type"),
  value    = c(gex_in, meta_out, cell_type)) |>
  knitr::kable(format = "simple") |>
  cat(sep = "\n")
message("============================\n")

##  Get gene IDs  -----------------------------------------------------------------------
message("Loading gene expression file ...")
gex_tbl <- read_tsv(gex_in, show_col_types = FALSE)

gene_ids <- gex_tbl |> pull(TargetID)

##  Query Ensembl for TSS  --------------------------------------------------------------
message("Using biomaRt to get TSS for ", length(gene_ids), " genes ...")
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

tss_tbl <- getBM(
  attributes = c("ensembl_gene_id", "chromosome_name", "start_position", "end_position", "strand"),
  filters = "ensembl_gene_id",
  values = gene_ids,
  mart = mart) |>
  as_tibble() |>
  dplyr::mutate(
    phenotype_id  = ensembl_gene_id,
    group_id      = ensembl_gene_id,
    gene_id       = ensembl_gene_id,
    chromosome    = chromosome_name,
    phenotype_pos = if_else(strand == 1, start_position, end_position)) |>
  distinct()
message(nrow(tss_tbl), " distinct genes returned by biomaRt.")

##  Write output  -----------------------------------------------------------------------
message("Writing gene metadata file ...")
tss_tbl |>
  dplyr::select(phenotype_id, group_id, gene_id, chromosome, phenotype_pos, strand) |>
  write_tsv(meta_out)

message("Done.")

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
