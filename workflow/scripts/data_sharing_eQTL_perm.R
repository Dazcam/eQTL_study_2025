#--------------------------------------------------------------------------------------
#
#    eQTL permutation files for data sharing
#
#--------------------------------------------------------------------------------------
#
# Pipeline:  18DATA-SHARING | Rule: data_perm
#            Upstream:   05TENSORQTL tensorqtl_perm (final norm method, genPC and
#                        exp_pc_map expPC per cell type); convert_genotypes (pvar)
#            Downstream: mk_eqtl_tar
#
# Purpose:   For each cell type, take the tensorQTL permutation results (top SNP
#            per gene), then:
#            - Drop SNP IDs of '.' and SNPs absent from the pvar
#            - Add CHROM, POS, REF, ALT from the pvar and gene symbols
#            - Count (but keep) rows with no gene symbol and duplicate gene-SNP pairs
#            - Name each output with the manuscript cell type label
#
# Inputs:    eqtl         Permutation results, one per cell type (cell_types order)
#            pvar         Genotype pvar
#            gene_lookup  BioMart gene annotation (hg38)
#
# Outputs:   tsv          {cell_type}_cis_eQTL_perm.tsv.gz, one per cell type
#                         (manuscript labels, cell_types order)
#            summary      Rows in, removed and counted per cell type, plus TOTAL
#            dup_snps     SNP IDs listed more than once in the pvar
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

message("\n\neQTL permutation files for data sharing ...")

##  Load packages, functions and variables  -------------------------------------------
library(tidyverse)

# Load functions
source(file.path(snakemake@scriptdir, "functions.R"))
source(file.path(snakemake@scriptdir, "data_sharing_functions.R"))

# Input and output paths
eqtl_files       <- snakemake@input[["eqtl"]]
pvar_file        <- snakemake@input[["pvar"]]
gene_lookup_file <- snakemake@input[["gene_lookup"]]
out_files        <- snakemake@output[["tsv"]]
summary_file     <- snakemake@output[["summary"]]
dup_snps_file    <- snakemake@output[["dup_snps"]]

# Settings
cell_types <- snakemake@params[["cell_types"]]
labels     <- snakemake@params[["labels"]]

# Print what each variable is set to
message("\nVariables")
message("============================")
tibble(
  variable = c("eqtl_files", "pvar_file", "gene_lookup_file", "out_files",
               "summary_file", "dup_snps_file", "cell_types", "labels"),
  value    = c(paste(eqtl_files, collapse = ", "), pvar_file, gene_lookup_file,
               paste(out_files, collapse = ", "), summary_file, dup_snps_file,
               paste(cell_types, collapse = ", "),
               paste(names(labels), unlist(labels), sep = " -> ", collapse = ", "))) |>
  knitr::kable(format = "simple") |>
  cat(sep = "\n")
message("============================\n")

## Check inputs and outputs line up  ---------------------------------------------------
cell_type_labels <- relabel_cell_types(cell_types, labels)

if (anyDuplicated(cell_type_labels)) {
  stop("Relabelling produced duplicate cell type labels: ",
       paste(cell_type_labels[duplicated(cell_type_labels)], collapse = ", "))
}
if (length(eqtl_files) != length(cell_types) || length(out_files) != length(cell_types)) {
  stop("Expected one input and one output per cell type (", length(cell_types), "); got ",
       length(eqtl_files), " inputs and ", length(out_files), " outputs")
}
if (!all(startsWith(basename(out_files), paste0(cell_type_labels, "_")))) {
  stop("Output files are not in cell_types order")
}

## Load pvar and gene lookup  ----------------------------------------------------------
message("Loading pvar file ...")
pvar_list <- read_pvar(pvar_file)
write_tsv(pvar_list$dup_snps, dup_snps_file)

gene_lookup <- read_gene_lookup(gene_lookup_file)

## Annotate and write one file per cell type  ------------------------------------------
message("eQTL permutation export: ", length(cell_types), " cell types")
drop_summary <- tibble()

for (i in seq_along(cell_types)) {

  message("\n--- ", cell_types[i], " -> ", cell_type_labels[i], " ---")
  message("Loading eQTL perm file: ", eqtl_files[i])
  eqtl_raw <- read_tsv(eqtl_files[i], show_col_types = FALSE) |>
    dplyr::select(ensembl_id = phenotype_id, SNP = variant_id, af,
                  pval_nominal, pval_beta, qval, slope, slope_se)

  annotated <- annotate_eqtl(eqtl_raw, pvar_list$pvar, gene_lookup,
                             extra_cols = c("pval_beta", "qval"))

  message("Writing: ", out_files[i])
  write_tsv(annotated$tbl, out_files[i])

  drop_summary <- bind_rows(drop_summary,
                            mutate(annotated$counts, cell_type = cell_type_labels[i], .before = 1))
}

## Summary  ----------------------------------------------------------------------------
write_tsv(summarise_row_drops(drop_summary), summary_file)
message("Wrote row-drop summary: ", summary_file)

message("Done.")

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
