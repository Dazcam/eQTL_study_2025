#--------------------------------------------------------------------------------------
#
#    Prepare expression and gene coordinate files for FUSION
#
#--------------------------------------------------------------------------------------
#
# Pipeline:  11TWAS | Rule: prep_exp_data (one job per cell type)
#            Upstream:   05TENSORQTL prep_input (quantile-normalised expression bed)
#            Downstream: compute_weights, make_pos_file
#
# Purpose:   - Transpose the tensorQTL expression bed to PLINK phenotype format
#              (FID, IID, one column per gene) for FUSION.compute_weights.R
#            - Get gene coordinates from Ensembl BioMart for every gene in the
#              expression file
#
# Inputs:    expr_in    Quantile-normalised expression bed: #Chr, start, end, TargetID,
#                       then one column per sample. Same input as tensorQTL. Not
#                       covariate-corrected; covariates go to FUSION via --covar
#
# Outputs:   expr_out   PLINK-format expression file: FID, IID, ~18K Ensembl IDs
#            coord_out  Gene coordinate file: chr, start, end, gene_id
#
# Notes:     run_test = TRUE keeps a random subset of 50 genes for testing.
#            Coordinates come from the current Ensembl release (www.ensembl.org).
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

message("\n\nCreating PLINK-ready expression data and gene coordinate file for FUSION ...")

##  Load packages, functions and variables  -------------------------------------------
library(tidyverse)
library(biomaRt)

# Load functions
source(file.path(snakemake@scriptdir, "twas_functions.R"))

# Input and output paths
expr_in <- snakemake@input[[1]]
expr_out <- snakemake@output[["exp"]]
coord_out <- snakemake@output[["coord"]]
run_test <- isTRUE(snakemake@params[["run_test"]]) # Option to run on 50 genes for testing

# Print what each variable is set to
message("\nVariables")
message("============================")
tibble(
  variable = c("expr_in", "expr_out", "coord_out", "run_test"),
  value    = c(expr_in, expr_out, coord_out, run_test)) |>
  knitr::kable(format = "simple") |>
  cat(sep = "\n")
message("============================\n")

##  Read expression data  -------------------------------------------------------------
message("Reading expression file ...")
expr_tbl <- read_tsv(expr_in, show_col_types = FALSE)

if (run_test) {
  message("Running with only a subset of genes for testing ...")
  message("Genes before: ", nrow(expr_tbl))
  expr_tbl <- dplyr::slice_sample(expr_tbl, n = 50)
  message("Genes after: ", nrow(expr_tbl))
}

##  Convert to PLINK format  ----------------------------------------------------------
coord_cols <- c("#Chr", "start", "end", "TargetID")
sample_cols <- setdiff(colnames(expr_tbl), coord_cols)

message("Converting expression data to PLINK format ...")
expr_plink_tbl <- expr_tbl |>
  dplyr::select(TargetID, all_of(sample_cols)) |>
  to_plink_pheno(id_col = "TargetID")

message("Saving PLINK phenotype file ...")
write_tsv(expr_plink_tbl, expr_out)

##  Get gene coordinates from BioMart  ------------------------------------------------
message("Setting BioMart mart ...")
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "www.ensembl.org")

message("Getting coordinates for genes ...")
gene_ids <- expr_tbl$TargetID
gene_coord <- getBM(
  attributes = c("chromosome_name", "start_position", "end_position", "ensembl_gene_id"),
  filters = "ensembl_gene_id",
  values = gene_ids,
  mart = mart
) |>
  as_tibble() |>
  dplyr::rename(
    chr = chromosome_name,
    start = start_position,
    end = end_position,
    gene_id = ensembl_gene_id
  ) |>
  mutate(chr = str_remove(chr, "^chr")) |>
  dplyr::select(chr, start, end, gene_id)

# Check for unmatched genes
missing_genes <- setdiff(gene_ids, gene_coord$gene_id)
if (length(missing_genes) > 0) {
  message(str_glue("{length(missing_genes)} genes not found in BioMart: {str_c(head(missing_genes, 5), collapse = ', ')}..."))
}

##  Save gene coordinates  ------------------------------------------------------------
message("Saving gene coordinate file ...")
write_tsv(gene_coord, coord_out)

message(str_glue("Generated {expr_out} with {nrow(expr_plink_tbl)} samples and {length(gene_ids)} genes\n"))
message(str_glue("Generated {coord_out} with {nrow(gene_coord)} genes\n"))
message("Done.")

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
