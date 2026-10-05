#--------------------------------------------------------------------------------------
#
#    Generate Supplementary Table 3: subcluster-specific eGenes for manuscript
#
#--------------------------------------------------------------------------------------

# Ports the "Subcluster-specific gene table" from devspec_report.Rmd's
# "Subcluster-specific eQTL" tab: FDR-sig eGenes in an L2 subcluster but not
# their L1 parent, with L1-vs-L2 slope/SE/p-value shown side by side.
#
# Adds, per manuscript request:
#   - Defensive deduplication on the natural key (ENSG, eSNP, cell_type_L2)
#     -- a gene can legitimately have more than one significant eSNP (not a
#     duplicate), so dedup is keyed on the full (gene, SNP, subcluster)
#     triple, not just the gene.
#   - Logs how many distinct genes are subcluster-specific in exactly one
#     subpopulation (cell_type_L2) vs more than one.
#   - Cell types relabelled Glu-UL/Glu-DL -> Glu-A/Glu-B, matching every
#     other manuscript table/plot in this pipeline (the Rmd itself doesn't
#     do this -- flagging in case it should too).
#
# Independent of devspec_report.Rmd (same reasoning as smr_sig_tbl vs
# smr_report.Rmd / ctwas_sig_tbl vs ctwas_report.Rmd): reads directly from
# the eqtl_effects RDS rather than the report's output, so it doesn't matter
# which pipeline runs first.

## Info  ------------------------------------------------------------------------------

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

message("\n\nGenerating Supplementary Table 3 (subcluster-specific eGenes) for the manuscript ...")

# -------------------------------------------------------------------------------------
suppressMessages({
  library(tidyverse)
  library(openxlsx)
})

# --- Set variables
eqtl_effects_file <- snakemake@params[["eqtl_effects"]]
out_file           <- snakemake@output[[1]]

# eqtl_effects_file <- "../results/19DEV-SPECIFICITY/eqtl_effects/eqtl_effect_sizes.rds"
# out_file           <- "../results/13MANUSCRIPT_TABLES/subcluster_specific_eqtl_tbl.xlsx"

# --- Relabel cell type names for display: Glu-UL -> Glu-A, Glu-DL -> Glu-B
relabel_cell_type <- function(x) {
  x <- str_replace(x, "Glu-UL", "Glu-A")
  x <- str_replace(x, "Glu-DL", "Glu-B")
  x
}

# --- Read and select/rename, exactly matching the Rmd's chunk
message("Reading: ", eqtl_effects_file)
eqtl_effects   <- read_rds(eqtl_effects_file)
subcluster_tbl <- eqtl_effects$subcluster_specific

n_before <- nrow(subcluster_tbl)
message("Rows in subcluster_specific (raw): ", n_before)

final_tbl <- subcluster_tbl |>
  select(Gene = gene_symbol, ENSG = phenotype_id, eSNP = variant_id,
         cell_type_L1 = parent_L1, Slope_L1 = slope_parent, SE_L1 = se_parent, Pval_L1 = pval_nom_parent,
         cell_type_L2 = cell_type_L2, Slope_L2 = slope_L2, SE_L2 = se_L2, Pval_L2 = pval_nom_L2) |>
  mutate(across(where(is.numeric), ~ signif(.x, 3))) |>
  mutate(cell_type_L1 = relabel_cell_type(cell_type_L1),
         cell_type_L2 = relabel_cell_type(cell_type_L2))

# --- Deduplicate on the natural key: (gene, eSNP, subcluster). A gene can
# legitimately have more than one significant eSNP, or be specific in more
# than one subcluster - neither of those is a duplicate. Only an exact
# repeat of the same (ENSG, eSNP, cell_type_L2) triple is.
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

# --- Log: how many distinct genes are subcluster-specific in exactly one
# subpopulation (cell_type_L2) vs more than one
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

# --- Export
message("")
message("Writing Supplementary Table 3: ", out_file)
write.xlsx(final_tbl, file = out_file, overwrite = TRUE,
           headerStyle = createStyle(textDecoration = "bold"))
message("Export complete.")

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
