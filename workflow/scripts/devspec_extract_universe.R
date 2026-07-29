#--------------------------------------------------------------------------------------
#
#    Developmental specificity: fetal eGene/eSNP universe extraction
#
#--------------------------------------------------------------------------------------

#  Builds the master long-format table of FDR<0.05 eGene-eSNP pairs across all
#  L1 + L2 fetal cell types, tagging each row with:
#    - level (L1 / L2)
#    - parent_cell_type (L2 -> parent L1 label; L1 -> itself)
#    - fugita_label_match (a priori adult Fujita cell type, via parent for L2;
#      NA where no defensible adult equivalent exists, e.g. NPC)
#
#  Step 1 of the developmental-specificity pipeline (independent for now;
#  designed to fold into replication_pi1_enrichment.R later if retained for the paper).

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

message("\n\nExtracting fetal eGene/eSNP universe (L1 + L2) ...")

suppressPackageStartupMessages({
  library(tidyverse)
})

# Input and output paths ---------------------------------------------------------------
qtl_perm_files  <- snakemake@input[["qtl_perm_files"]]
output_universe <- snakemake@output[[1]]
cell_types      <- snakemake@params[["cell_types"]]
geno_pc         <- snakemake@params[["geno_pc"]]
norm_method     <- snakemake@params[["norm_method"]]

names(qtl_perm_files) <- cell_types  # preserve cell_type <-> file correspondence

# L1 cell types -- fixed by cell-type annotation, not a tunable pipeline parameter.
cell_types_L1 <- c("Glu-UL", "Glu-DL", "GABA", "NPC", "MG", "OPC", "Endo-Peri")

# A priori biological correspondence between fetal L1 cell types and Fujita adult
# cell types. Fixed mapping -- not a tunable pipeline parameter, hence hardcoded
# here rather than in config. NPC has no defensible adult equivalent (no progenitor
# population in the adult cortex), so it is excluded from label/best-match
# comparisons and only included in the "any" matching tier downstream.
fugita_label_match <- c(
  "Glu-UL"    = "Exc",
  "Glu-DL"    = "Exc",
  "GABA"      = "Inh",
  "NPC"       = NA_character_,
  "MG"        = "Mic",
  "OPC"       = "OPC",
  "Endo-Peri" = "End"
)

# Check variable assignment
message("\nVariables")
cat("============================")
tibble(
  variable = c("n_cell_types", "geno_pc", "norm_method", "output_universe"),
  value    = c(length(cell_types), geno_pc, norm_method, output_universe)
) |>
  knitr::kable(format = "simple", align = "l") |>
  print()
message("\n============================\n")

# Functions ----------------------------------------------------------------------------

# Strip L2 sub-cluster suffix to get parent L1 cell type
# (same convention used in tensorqtl_report.Rmd's egenes_L2 %>% mutate(parent = ...))
get_parent_cell_type <- function(cell_type, cell_types_L1) {
  if (cell_type %in% cell_types_L1) return(cell_type)
  sub("-[^-]+$", "", cell_type)
}

# Read FDR<0.05 eGene-eSNP pairs for a single cell type
# (tensorQTL perm output already reports one lead eSNP per eGene)
read_sig_eqtl <- function(cell_type, file_path) {
  message("Reading sig. eQTL for: ", cell_type, " (", file_path, ")")

  if (!file.exists(file_path)) {
    message("  File not found, skipping: ", file_path)
    return(tibble())
  }

  read_tsv(file_path, show_col_types = FALSE) %>%
    filter(!is.na(qval) & qval < 0.05) %>%
    select(variant_id, phenotype_id, slope, qval) %>%
    mutate(cell_type = cell_type)
}

# Main -----------------------------------------------------------------------------------

eqtl_universe <- map2_dfr(names(qtl_perm_files), qtl_perm_files, read_sig_eqtl)

# Dedup safeguard: one row per (cell_type, phenotype_id)
n_before <- nrow(eqtl_universe)
eqtl_universe <- eqtl_universe %>% distinct(cell_type, phenotype_id, .keep_all = TRUE)
n_after <- nrow(eqtl_universe)

if (n_before != n_after) {
  message("WARNING: removed ", n_before - n_after,
          " duplicate (cell_type, phenotype_id) rows.")
} else {
  message("No duplicate (cell_type, phenotype_id) rows found.")
}

# Tag level, parent cell type, and a priori Fujita label-match
eqtl_universe <- eqtl_universe %>%
  mutate(
    level = ifelse(cell_type %in% cell_types_L1, "L1", "L2"),
    parent_cell_type = map_chr(cell_type, get_parent_cell_type, cell_types_L1 = cell_types_L1),
    fugita_label_match = fugita_label_match[parent_cell_type]
  )

# Sanity check against manuscript's reported ~3,121 unique eGene figure
n_unique_egenes <- n_distinct(eqtl_universe$phenotype_id)
message("\nTotal unique eGenes across all ", length(cell_types),
        " cell types (L1+L2): ", n_unique_egenes,
        " (manuscript reports ~3,121 -- confirm these match)")

message("\nRows per cell type:")
print(eqtl_universe %>% count(cell_type, level))

# Write single tracked output ----------------------------------------------------------
write_rds(eqtl_universe, output_universe)
message("\nEgene/eSNP universe written to: ", output_universe)
