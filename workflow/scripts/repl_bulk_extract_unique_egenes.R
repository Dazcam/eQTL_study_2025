#--------------------------------------------------------------------------------------
#
#    Extract FDR-significant eGenes across all cell types
#
#--------------------------------------------------------------------------------------
#
# Pipeline:  21REPLICATION-BULK | Rule: extract_unique_egenes (single job)
#            Upstream:   05TENSORQTL tensorqtl_perm (final expPC per cell type)
#            Downstream: overlap_obrien; isoform_eqtl_comparison; repl_bulk_report;
#                        14 bulk_egene_overlaps_plt (Supp Fig 9)
#
# Purpose:   - Read the tensorQTL permutation output for every cell type in
#              config (L1, L2 and pseudotime bins)
#            - Keep eGenes at qval < eqtl_fdr, one lead eSNP per eGene
#            - Tag each cell type as L1, L2 or pseudotime
#            - Report the per-cell-type and total unique eGene counts
#
# Inputs:    perm_files   tensorQTL perm results, one per cell type (same order
#                         as cell_types)
#            cell_types   config['cell_types']
#            l1_types     config['fetal_l1_types'] (L1 level tag)
#            eqtl_fdr     qval threshold for an eGene
#
# Outputs:   output       TSV, one row per (cell_type, eGene): phenotype_id,
#                         variant_id, slope, slope_se, pval_nominal, qval,
#                         cell_type, level
#
# Notes:     Rows are kept per cell type (not deduplicated across cell types);
#            downstream steps take the unique set from this table.
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

message("\n\nExtracting FDR-significant eGenes across all cell types ...")

##  Load packages, functions and variables  -------------------------------------------
library(tidyverse)

# Input and output paths
perm_files    <- snakemake@input[["perm_files"]]
output_file   <- snakemake@output[[1]]
cell_types    <- snakemake@params[["cell_types"]]
geno_pc       <- snakemake@params[["geno_pc"]]
norm_method   <- snakemake@params[["norm_method"]]
cell_types_L1 <- unlist(snakemake@params[["l1_types"]])
eqtl_fdr      <- as.numeric(snakemake@params[["eqtl_fdr"]])

names(perm_files) <- cell_types

# Print what each variable is set to
message("\nVariables")
message("============================")
tibble(
  variable = c("n_cell_types", "geno_pc", "norm_method", "cell_types_L1", "eqtl_fdr",
               "output_file"),
  value    = c(length(cell_types), geno_pc, norm_method,
               paste(cell_types_L1, collapse = ", "), eqtl_fdr, output_file)) |>
  knitr::kable(format = "simple") |>
  cat(sep = "\n")
message("============================\n")

##  Functions  -------------------------------------------------------------------------
is_pseudotime <- function(ct) str_detect(ct, "-to-")

get_level <- function(ct, cell_types_L1) {
  if (ct %in% cell_types_L1) return("L1")
  if (is_pseudotime(ct)) return("pseudotime")
  "L2"
}

# FDR-significant eGenes for one cell type (tensorQTL perm output already has
# one lead eSNP per eGene)
read_sig_egenes <- function(cell_type, file_path) {
  message("Reading sig. eGenes for: ", cell_type, " (", file_path, ")")

  if (!file.exists(file_path)) {
    message("  File not found, skipping: ", file_path)
    return(tibble())
  }

  dat <- read_tsv(file_path, show_col_types = FALSE) |>
    filter(!is.na(qval) & qval < eqtl_fdr) |>
    dplyr::select(phenotype_id, variant_id, slope, slope_se, pval_nominal, qval) |>
    mutate(cell_type = cell_type)

  message("  Sig. eGenes (qval < ", eqtl_fdr, "): ", nrow(dat))
  dat
}

##  Read significant eGenes  -----------------------------------------------------------
per_celltype_egenes <- map2_dfr(names(perm_files), perm_files, read_sig_egenes)

n_before <- nrow(per_celltype_egenes)
per_celltype_egenes <- per_celltype_egenes |> distinct(cell_type, phenotype_id, .keep_all = TRUE)
n_after <- nrow(per_celltype_egenes)

if (n_before != n_after) {
  message("Removed ", n_before - n_after, " duplicate (cell_type, phenotype_id) rows.")
} else {
  message("No duplicate (cell_type, phenotype_id) rows found.")
}

per_celltype_egenes <- per_celltype_egenes |>
  mutate(level = map_chr(cell_type, get_level, cell_types_L1 = cell_types_L1))

##  Summarise  -------------------------------------------------------------------------
message("\neGene counts per cell type:\n")
per_celltype_egenes |>
  count(cell_type, level, name = "n_egenes") |>
  arrange(match(cell_type, cell_types)) |>
  knitr::kable(format = "simple", align = "l") |>
  print()

message("\nCell types per level:")
per_celltype_egenes |>
  distinct(cell_type, level) |>
  count(level, name = "n_cell_types") |>
  print()

n_unique_egenes <- n_distinct(per_celltype_egenes$phenotype_id)
message("\nTotal unique eGenes across ", length(cell_types), " cell types: ", n_unique_egenes)

##  Write output  ----------------------------------------------------------------------
write_tsv(per_celltype_egenes, output_file)
message("\nPer-cell-type significant eGene rows written to: ", output_file)
message("Done.")

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
