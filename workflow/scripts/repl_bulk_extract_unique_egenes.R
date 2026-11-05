#--------------------------------------------------------------------------------------
#
#    Replication (bulk): unique eGene extraction across L1 + L2 + pseudotime
#
#--------------------------------------------------------------------------------------

#  Reads FDR<0.05 (qval<0.05) eGenes from every cell type's TensorQTL permutation
#  output -- L1, L2, AND pseudotime bins (27 cell types total)
#
#  Retains every per-cell-type significant row (not deduplicated) alongside the
#  final unique eGene set
#
#  Step 1 of the replication_bulk pipeline.
#
#--------------------------------------------------------------------------------------

# Set up logging for Snakemake
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

message("\n\nExtracting unique FDR-sig eGenes across all L1+L2+pseudotime cell types ...")

suppressPackageStartupMessages({
  library(tidyverse)
})

# Input and output paths ---------------------------------------------------------------
perm_files  <- snakemake@input[["perm_files"]]
output_file <- snakemake@output[[1]]
cell_types  <- snakemake@params[["cell_types"]]
geno_pc     <- snakemake@params[["geno_pc"]]
norm_method <- snakemake@params[["norm_method"]]

names(perm_files) <- cell_types

# L1 cell types 
cell_types_L1 <- c("Glu-UL", "Glu-DL", "GABA", "NPC", "MG", "OPC", "Endo-Peri")

# Pseudotime bins 
is_pseudotime <- function(ct) str_detect(ct, "-to-")

get_level <- function(ct, cell_types_L1) {
  if (ct %in% cell_types_L1) return("L1")
  if (is_pseudotime(ct)) return("pseudotime")
  "L2"
}

# Check variable assignment
message("\nVariables")
cat("============================")
tibble(
  variable = c("n_cell_types", "geno_pc", "norm_method", "output_file"),
  value    = c(length(cell_types), geno_pc, norm_method, output_file)
) |>
  knitr::kable(format = "simple", align = "l") |>
  print()
message("\n============================\n")

# Functions ----------------------------------------------------------------------------

# Read FDR<0.05 eGenes for a single cell type (one row per gene -- TensorQTL
# perm output already reports one lead eSNP per eGene)
read_sig_egenes <- function(cell_type, file_path) {
  message("Reading sig. eGenes for: ", cell_type, " (", file_path, ")")

  if (!file.exists(file_path)) {
    message("  File not found, skipping: ", file_path)
    return(tibble())
  }

  dat <- read_tsv(file_path, show_col_types = FALSE) %>%
    filter(!is.na(qval) & qval < 0.05) %>%
    select(phenotype_id, variant_id, slope, slope_se, pval_nominal, qval) %>%
    mutate(cell_type = cell_type)

  message("  Sig. eGenes (qval<0.05): ", nrow(dat))
  dat
}

# Main -----------------------------------------------------------------------------------

per_celltype_egenes <- map2_dfr(names(perm_files), perm_files, read_sig_egenes)

# Dedup 
n_before <- nrow(per_celltype_egenes)
per_celltype_egenes <- per_celltype_egenes %>% distinct(cell_type, phenotype_id, .keep_all = TRUE)
n_after <- nrow(per_celltype_egenes)

if (n_before != n_after) {
  message("WARNING: removed ", n_before - n_after,
          " duplicate (cell_type, phenotype_id) rows.")
} else {
  message("No duplicate (cell_type, phenotype_id) rows found.")
}

# Tag level (L1 / L2 / pseudotime)
per_celltype_egenes <- per_celltype_egenes %>%
  mutate(level = map_chr(cell_type, get_level, cell_types_L1 = cell_types_L1))

# --- Per-cell-type eGene counts (as requested: log each cell type's total) ------------
message("\n--- eGene counts per cell type ---\n")

per_celltype_summary <- per_celltype_egenes %>%
  count(cell_type, level, name = "n_egenes") %>%
  arrange(match(cell_type, cell_types))

per_celltype_summary %>%
  knitr::kable(format = "simple", align = "l") %>%
  print()

message("\nTotals by level:")
per_celltype_egenes %>%
  distinct(cell_type, level) %>%
  count(level, name = "n_cell_types") %>%
  print()

# --- Unique eGene set across ALL cell types (L1 + L2 + pseudotime pooled) -------------
unique_egenes <- per_celltype_egenes %>%
  distinct(phenotype_id) %>%
  arrange(phenotype_id)

n_unique_egenes <- nrow(unique_egenes)

message("\n========================================")
message("TOTAL unique eGenes across all ", length(cell_types),
        " cell types (L1+L2+pseudotime): ", n_unique_egenes)
message("========================================\n")

# --- Write TSV --------------
write_tsv(per_celltype_egenes, output_file)
message("\nPer-cell-type significant eGene rows (n_unique_egenes=", n_unique_egenes,
        ") written to: ", output_file)
