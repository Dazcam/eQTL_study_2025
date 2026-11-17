#--------------------------------------------------------------------------------------
#
#    Subcluster-specific eGenes (Supp Table 3)
#
#--------------------------------------------------------------------------------------
#
# Pipeline:  16MANUSCRIPT-TABLES | Rule: subcluster_specific_tbl
#            Upstream:   15 eqtl_effects (eqtl_effect_sizes.rds)
#            Downstream: Supp Table 3
#
# Purpose:   - Read the subcluster-specific eGenes computed by 15 (FDR < 0.05 in an
#              L2 subcluster but not in its L1 parent), with L1 and L2 slope, SE
#              and nominal P side by side, as in eqtl_sharing_report.Rmd
#            - Round numeric columns to 3 significant figures
#            - Remove exact repeats of (ENSG, eSNP, cell_type_L2)
#            - Log how many genes are specific to one subcluster vs more
#            - Manuscript cell type labels
#
# Inputs:    eqtl_effects      15 eqtl_effects RDS ($subcluster_specific)
#
# Outputs:   out_file          xlsx, one sheet
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

message("\n\nBuilding subcluster-specific eGene table ...")

##  Load packages, functions and variables  -------------------------------------------
suppressMessages({
  library(tidyverse)
  library(openxlsx)
})

# Load functions (relabel_cell_types)
source(file.path(snakemake@scriptdir, "functions.R"))

# Inputs, outputs and settings
eqtl_effects_file <- snakemake@params[["eqtl_effects"]]
manuscript_labels <- snakemake@params[["labels"]]
out_file          <- snakemake@output[[1]]

# Print what each variable is set to
message("\nVariables")
message("============================")
tibble(
  variable = c("eqtl_effects_file", "labels", "out_file"),
  value    = c(eqtl_effects_file,
               paste(names(manuscript_labels), unlist(manuscript_labels), sep = " -> ", collapse = ", "),
               out_file)) |>
  knitr::kable(format = "simple") |>
  cat(sep = "\n")
message("============================\n")

##  Subcluster-specific eGenes  -------------------------------------------------------
message("Reading: ", eqtl_effects_file)
eqtl_effects   <- read_rds(eqtl_effects_file)
subcluster_tbl <- eqtl_effects$subcluster_specific

n_before <- nrow(subcluster_tbl)
message("Rows in subcluster_specific (raw): ", n_before)

final_tbl <- subcluster_tbl |>
  dplyr::select(Gene = gene_symbol, ENSG = phenotype_id, eSNP = variant_id,
         cell_type_L1 = parent_L1, Slope_L1 = slope_parent, SE_L1 = se_parent, Pval_L1 = pval_nom_parent,
         cell_type_L2 = cell_type_L2, Slope_L2 = slope_L2, SE_L2 = se_L2, Pval_L2 = pval_nom_L2) |>
  mutate(across(where(is.numeric), ~ signif(.x, 3))) |>
  mutate(cell_type_L1 = relabel_cell_types(cell_type_L1, manuscript_labels),
         cell_type_L2 = relabel_cell_types(cell_type_L2, manuscript_labels))

# Deduplicate on (ENSG, eSNP, cell_type_L2): a gene can have more than one
# significant eSNP, or be specific in more than one subcluster
n_dupe_keys <- final_tbl |> count(ENSG, eSNP, cell_type_L2) |> filter(n > 1) |> nrow()
final_tbl <- final_tbl |> distinct(ENSG, eSNP, cell_type_L2, .keep_all = TRUE)
n_after <- nrow(final_tbl)

if (n_before == n_after) {
  message("No duplicate (ENSG, eSNP, cell_type_L2) rows found (", n_before, " rows unchanged).")
} else {
  message("Removed ", n_before - n_after, " duplicate row(s) (", n_dupe_keys,
          " distinct duplicated (ENSG, eSNP, cell_type_L2) keys): ",
          n_before, " -> ", n_after, " rows.")
}

##  Subpopulation specificity summary  ------------------------------------------------
gene_subpop_counts <- final_tbl |>
  distinct(ENSG, cell_type_L2) |>
  count(ENSG, name = "n_subpops")

n_genes_total     <- nrow(gene_subpop_counts)
n_genes_single    <- sum(gene_subpop_counts$n_subpops == 1)
n_genes_multi     <- sum(gene_subpop_counts$n_subpops > 1)

message("")
message("=== Subpopulation specificity summary ===")
message("Distinct subcluster-specific genes: ", n_genes_total)
message("  Detected in exactly ONE subpopulation: ", n_genes_single,
        " (", round(100 * n_genes_single / n_genes_total, 1), "%)")
message("  Detected in MORE THAN ONE subpopulation: ", n_genes_multi,
        " (", round(100 * n_genes_multi / n_genes_total, 1), "%)")
if (n_genes_multi > 0) {
  breakdown <- gene_subpop_counts |> filter(n_subpops > 1) |> count(n_subpops)
  for (i in seq_len(nrow(breakdown))) {
    message("    - in exactly ", breakdown$n_subpops[i], " subpopulations: ", breakdown$n[i], " gene(s)")
  }
}

##  Write xlsx  -----------------------------------------------------------------------
message("")
message("Writing: ", out_file)
write.xlsx(final_tbl, file = out_file, overwrite = TRUE,
           headerStyle = createStyle(textDecoration = "bold"))
message("Done.")

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
