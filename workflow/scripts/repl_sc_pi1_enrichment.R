#--------------------------------------------------------------------------------------
#
#    pi1 replication of sig. eQTL between datasets
#
#--------------------------------------------------------------------------------------
#
# Pipeline:  07REPLICATION-SC | Rules: pi1_enrich_internal (cell type x cell type),
#                                       pi1_enrich_jang (cell type x Jang L1 type)
#            Upstream:   05TENSORQTL tensorqtl_perm; cat_nom_qtl;
#                        Jang et al. 2026 SingleBrain L1 full + top association files
#            Downstream: 07 report (internal and Jang pi1 tabs; Fig 3D)
#
# Purpose:   For one query cell type and one reference, run compute_pi1() both ways:
#            - Forward: sig. eQTL of our cell type (qval < 0.05) looked up in the
#              reference's nominal associations
#            - Reverse: sig. eQTL of the reference looked up in our cell type's
#              nominal associations
#            Modes, chosen from the output path (ref_name):
#            - Internal: reference is another of our cell types; self-comparisons
#              write NA dummy outputs so Snakemake can track them
#            - Jang: reference is a Jang et al. 2026 adult L1 cell type
#
# Inputs:    public_all   reference nominal associations
#            public_top   reference top associations (sig. eQTL)
#            qtl_all      our nominal associations for the cell type (cat_nom_qtl)
#            qtl_top      our tensorQTL permutation output for the cell type
#
# Outputs:   enrich       one-row tibble: pi1, prop_replicating, prop_same_direction
#                         and overlap_count, forward and reverse (_cell_ref)
#            pi1          rds list(forward, reverse) of full compute_pi1() results
#
# Notes:     Jang reference:
#            - Forward p-value is Random_P (qvalue() needs raw p-values to estimate
#              pi0; random-effects side used for significance throughout)
#            - Beta is relative to `Allele`, not always alt; sign-normalised at read
#              time (align_jang_beta()), which prop_same_direction depends on
#            - Placeholder-cohort rows (beta = 0, sd = 0 in a cohort, giving
#              Fixed_P == Random_P == 1 exactly) are removed before qvalue(), as the
#              spike at p = 1 would bias pi0 towards 1 (filter_jang_artifacts())
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

message("\n\nCalculating pi1 enrichments: sig. sn-eQTL & public reference eQTL ...")

##  Load packages, functions and variables  -------------------------------------------
suppressPackageStartupMessages({
  library(tidyverse)
  library(qvalue)
})

source(file.path(snakemake@scriptdir, "repl_sc_functions.R"))

# Input and output paths
public_all_qtl <- snakemake@input[['public_all']]
public_top_qtl <- snakemake@input[['public_top']]
qtl_all <- snakemake@input[['qtl_all']]
qtl_top <- snakemake@input[['qtl_top']]
output_enrich <- snakemake@output[["enrich"]]
output_pi1 <- snakemake@output[["pi1"]]
cell_type <- snakemake@wildcards[["cell_type"]]

ref_cell_type <- if ("ref_cell_type" %in% names(snakemake@wildcards)) {
  snakemake@wildcards[["ref_cell_type"]]
} else if ("jang_cell_type" %in% names(snakemake@wildcards)) {
  snakemake@wildcards[["jang_cell_type"]]
} else {
  NA_character_
}

exp_pc <- if ("exp_pc" %in% names(snakemake@wildcards)) snakemake@wildcards[["exp_pc"]] else NA_character_
geno_pc <- if ("geno_pc" %in% names(snakemake@wildcards)) snakemake@wildcards[["geno_pc"]] else NA_character_
norm_method <- if ("norm_method" %in% names(snakemake@wildcards)) snakemake@wildcards[["norm_method"]] else NA_character_
ref_name <- str_split_i(output_pi1, "/", 4)

# Print what each variable is set to
message("\nVariables")
message("============================")
tibble(
  variable = c("public_all_qtl", "public_top_qtl", "qtl_all", "qtl_top", "output_enrich",
               "output_pi1", "cell_type", "ref_cell_type", "exp_pc", "geno_pc", "norm_method",
               "ref_name"),
  value    = c(public_all_qtl, public_top_qtl, qtl_all, qtl_top, output_enrich,
               output_pi1, cell_type, ref_cell_type, exp_pc, geno_pc, norm_method,
               ref_name)) |>
  knitr::kable(format = "simple") |>
  cat(sep = "\n")
message("============================\n")

##  Main function  ---------------------------------------------------------------------
run_pi1_enrichment <- function(cell_type, public_all_qtl, public_top_qtl,
                               qtl_all, qtl_top, output_enrich, output_pi1,
                               ref_name, ref_cell_type = NA_character_) {

  if (!is.na(ref_cell_type) && cell_type == ref_cell_type &&
      str_detect(ref_name, "internal")) {
    message("Skipping identical cell_type/ref_cell_type comparison: ", cell_type)

    enrichment_results <- tibble(
      cell_type = cell_type,
      ref_cell_type = ref_cell_type,
      pi1 = NA_real_,
      prop_replicating = NA_real_,
      prop_same_direction = NA_real_,
      overlap_count = NA_integer_,
      pi1_cell_ref = NA_real_,
      prop_replicating_cell_ref = NA_real_,
      prop_same_direction_cell_ref = NA_real_,
      overlap_count_cell_ref = NA_integer_
    )

    pi1_results <- list(forward = NULL, reverse = NULL)

    write_rds(enrichment_results, output_enrich)
    write_rds(pi1_results, output_pi1)
    message("Wrote dummy outputs for skipped comparison.")
    return(invisible(NULL))
  }

  # Load reference
  if (str_detect(ref_name, "jang")) {
    message("Loading Jang et al. 2026 (SingleBrain) full-association eQTL...")

    public_full_raw <- read_tsv(
      public_all_qtl,
      col_select = c(feature, variant_id, ref, alt, Allele, fixed_beta, Random_P, Fixed_P),
      show_col_types = FALSE
    )

    # Drop placeholder-cohort rows (Fixed_P == Random_P == 1); see Notes
    public_full_raw <- filter_jang_artifacts(public_full_raw, "Jang full-association")

    # Jang beta is relative to Allele; align to alt (needed for prop_same_direction)
    public_full <- public_full_raw %>%
      mutate(
        phenotype_id = str_extract(feature, "ENSG[0-9]+"),
        slope_ref = align_jang_beta(fixed_beta, Allele, alt),
        pval = suppressWarnings(as.numeric(Random_P))
      ) %>%
      select(variant_id, phenotype_id, pval, slope_ref)

    message("Loading Jang top-association (lead SNP) eQTL...")
    public_top_raw <- read_tsv(
      public_top_qtl,
      col_select = c(feature, variant_id, ref, alt, Allele, fixed_beta, Random_P, Fixed_P, qval),
      show_col_types = FALSE
    )

    # Same filter on the top file (expected to remove ~0 rows; logged as a check)
    public_top_raw <- filter_jang_artifacts(public_top_raw, "Jang top-association")

    # top table = gene's lead SNP with qval<0.05 (FDR across variants within
    # gene, then Storey's q across genes -- Jang's own two-stage correction)
    public_top <- public_top_raw %>%
      filter(!is.na(qval) & qval < 0.05) %>%
      mutate(
        phenotype_id = str_extract(feature, "ENSG[0-9]+"),
        slope_my = align_jang_beta(fixed_beta, Allele, alt)
      ) %>%
      select(variant_id, phenotype_id, slope_my)

    message("Jang: number of top (qval<0.05) eQTL retained: ", nrow(public_top))
    message("Jang: number of nominal (full) rows retained: ", nrow(public_full))
  } else {
    message("Internal mode: comparing ", cell_type, " ↔ ", ref_cell_type)

    public_full <- read_tsv(public_all_qtl, show_col_types = FALSE) %>%
      select(variant_id, phenotype_id, pval = pval_nominal, slope_ref = slope)

    public_top <- read_tsv(public_top_qtl, show_col_types = FALSE) %>%
      filter(qval < 0.05) %>%
      select(variant_id, phenotype_id, slope_my = slope)
  }

  # Load our cell type (query)
  message("Loading sig. eQTL (query) for cell type: ", cell_type)
  query_eqtl <- read_tsv(qtl_top, show_col_types = FALSE) %>%
    filter(qval < 0.05) %>%
    select(variant_id, phenotype_id, slope_my = slope)

  message("Loading all eQTL (query full) for cell type: ", cell_type)
  full_cell <- read_tsv(qtl_all, show_col_types = FALSE) %>%
    select(variant_id, phenotype_id,
           pval = pval_nominal, slope_ref = slope)

  enrichment_results <- tibble(
    cell_type = cell_type,
    ref_cell_type = ifelse(is.na(ref_cell_type), NA_character_, ref_cell_type),
    pi1 = NA_real_,
    prop_replicating = NA_real_,
    prop_same_direction = NA_real_,
    overlap_count = NA_integer_,
    pi1_cell_ref = NA_real_,
    prop_replicating_cell_ref = NA_real_,
    prop_same_direction_cell_ref = NA_real_,
    overlap_count_cell_ref = NA_integer_
  )

  pi1_results <- list(forward = NULL, reverse = NULL)

  # Forward: query sig -> reference full
  if (nrow(query_eqtl) > 0) {
    message("\n--- Forward enrichment: ", cell_type, " sig → ", ifelse(is.na(ref_cell_type), ref_name, paste0("(", ref_name, ":", ref_cell_type, ")")), " full ---")
    pi1_result <- compute_pi1(query_eqtl, public_full, ref_name)
    pi1_results$forward <-  pi1_result
    enrichment_results <- enrichment_results %>%
      mutate(
        pi1 = pi1_result$pi1,
        prop_replicating = pi1_result$prop_replicating,
        prop_same_direction = pi1_result$prop_same_direction,
        overlap_count = pi1_result$overlap_count
      )
  } else {
    message("No significant eQTLs found for ", cell_type, " (forward).")
  }

  # Reverse: reference sig -> query full
  if (nrow(full_cell) > 0) {
    message("\n--- Reverse enrichment: reference sig → ", cell_type, " full ---")
    pi1_result <- compute_pi1(public_top, full_cell, ifelse(is.na(ref_cell_type), ref_name, paste0(ref_name, ":", ref_cell_type)))
    pi1_results$reverse <-  pi1_result
    enrichment_results <- enrichment_results %>%
      mutate(
        pi1_cell_ref = pi1_result$pi1,
        prop_replicating_cell_ref = pi1_result$prop_replicating,
        prop_same_direction_cell_ref = pi1_result$prop_same_direction,
        overlap_count_cell_ref = pi1_result$overlap_count
      )
  } else {
    message("No full eQTL data found for ", cell_type, " (reverse).")
  }

  write_rds(enrichment_results, output_enrich)
  write_rds(pi1_results, output_pi1)
  message("\nEnrichment tbl written to: ", output_enrich)
  message("\nPi1 results written to: ", output_pi1)
}

##  Run with Snakemake params  ---------------------------------------------------------
if (exists("snakemake")) {
  run_pi1_enrichment(
    cell_type       = cell_type,
    public_all_qtl  = public_all_qtl,
    public_top_qtl  = public_top_qtl,
    qtl_all         = qtl_all,
    qtl_top         = qtl_top,
    output_enrich   = output_enrich,
    output_pi1      = output_pi1,
    ref_name        = ref_name,
    ref_cell_type   = ref_cell_type
  )
}
#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
