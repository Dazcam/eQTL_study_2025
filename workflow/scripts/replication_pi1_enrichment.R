#--------------------------------------------------------------------------------------
#
#    pi1 enrichments script
#
#--------------------------------------------------------------------------------------

#  Built to handle 6 options
#   - Our eQTL vs. O'Brien bulk eQTL 
#   - Our eQTL vs. Wen bulk eQTL
#   - Our eQTL vs. Bryois
#   - Our eQTL vs. Fugita
#   - Our eQTL vs. Jang et al. 2026 (SingleBrain)
#   - Our eQTL vs. Our eQTL (internal comparison)
#     - For latter R gens dummy file if cell_type == ref_cell_type for smk to track

#   - Tests: 
#     - Forward: pi1 enrichment of our sig. sn-eQTL in public nominal eQTL data
#     - Reverse: pi1 of public sig. sn-eQTL in our nominal sn-eQTL data

#--------------------------------------------------------------------------------------
#
#  ADDITION -- Jang et al. 2026 (SingleBrain) branch:
#
#  Jang publishes one file per cell type already (unlike Fugita's single
#  combined file with a `celltype` column to filter on) -- the right cell type
#  is selected purely via the jang_cell_type wildcard's file path, no
#  in-script filtering needed.
#
#  Forward test p-value: Random_P (not Fixed_P) -- consistent with using the
#  random-effects side for significance/evidence throughout the Fujita->Jang
#  migration elsewhere in this pipeline (classify_specificity.R,
#  heterogeneity_test.R, beta_correlation.R all use Random_FDR/qval).
#  qvalue() needs a raw p-value to estimate pi0 itself, so this must be
#  Random_P, not the already-FDR-adjusted Random_FDR/qval.
#
#  Sign normalization: Jang's beta is relative to the `Allele` column, not
#  always `alt` (unlike Fugita, where beta was implicitly always
#  ALT-relative). This directly affects prop_same_direction (computed by
#  compute_pi1() below via sign(slope_my) == sign(slope_ref)) -- skipping
#  this normalization would silently corrupt that output for any row where
#  Allele != alt. Applied at read time for both public_full and public_top.
#
#  Placeholder-cohort artifact filter: some cohorts encode "not tested for
#  this pair" as beta=0/sd=0 rather than NA, which gives that cohort
#  infinite weight (1/SE^2) in the meta-analysis and forces the pooled
#  result to Random_P == Fixed_P == 1 exactly -- a data artifact, not a
#  genuine null, confirmed by manual inspection (see git history / analysis
#  log for the diagnostic session). Rows matching this exact signature are
#  removed from public_full before it feeds compute_pi1(), with the removed
#  count logged -- otherwise this spike at p=1 biases qvalue()'s pi0
#  estimate toward 1 and masks real signal elsewhere in the distribution.
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

message("\n\nCalculating pi1 enrichments: sig. sn-eQTL & public reference eQTL ...")

suppressPackageStartupMessages({
  library(tidyverse)
  library(qvalue)
  library(readxl)
  library(cowplot)
})

# Input and output paths)
public_all_qtl <- snakemake@input[['public_all']]
public_top_qtl <- snakemake@input[['public_top']]
qtl_all <- snakemake@input[['qtl_all']]
qtl_top <- snakemake@input[['qtl_top']]
output_enrich <- snakemake@output[["enrich"]]
output_pi1 <- snakemake@output[["pi1"]]
cell_type <- snakemake@wildcards[["cell_type"]]

ref_cell_type <- if ("ref_cell_type" %in% names(snakemake@wildcards)) {
  snakemake@wildcards[["ref_cell_type"]]
} else if ("fugita_cell_type" %in% names(snakemake@wildcards)) {
  snakemake@wildcards[["fugita_cell_type"]]
} else if ("jang_cell_type" %in% names(snakemake@wildcards)) {
  snakemake@wildcards[["jang_cell_type"]]
} else {
  NA_character_
}

exp_pc <- if ("exp_pc" %in% names(snakemake@wildcards)) snakemake@wildcards[["exp_pc"]] else NA_character_
geno_pc <- if ("geno_pc" %in% names(snakemake@wildcards)) snakemake@wildcards[["geno_pc"]] else NA_character_
norm_method <- if ("norm_method" %in% names(snakemake@wildcards)) snakemake@wildcards[["norm_method"]] else NA_character_
ref_name <- str_split_i(output_pi1, "/", 4)

# Map ref_cell_type to Excel cell_type column
bryois_cell_map <- c(
  "Astrocytes"          = "Astrocytes",
  "Endothelial_cells"   = "Endothelial cells",
  "Excitatory_neurons"  = "Excitatory neurons",
  "Inhibitory_neurons"  = "Inhibitory neurons",
  "Microglia"           = "Microglia",
  "Oligodendrocytes"    = "Oligodendrocytes",
  "OPCs"                = "OPCs / COPs",
  "Pericytes"           = "Pericytes"
)




# Check variable assignment
message("\nVariables")
cat("============================")
tibble(
  variable = c("public_all_qtl", "public_top_qtl", "qtl_all", "qtl_top", "output_enrich",
               "output_pi1", "cell_type", "ref_cell_type", "exp_pc", "geno_pc", "norm_method", "ref_name"),
  value    = c(public_all_qtl, public_top_qtl, qtl_all, qtl_top, output_enrich,
               output_pi1, cell_type, ref_cell_type, exp_pc, geno_pc, norm_method, ref_name)
) |>
  knitr::kable(format = "simple", align = "l") |>
  print()
message("\n============================\n")


# Functions --------------------------------------------------------------------------
# Helper
compute_pi1 <- function(query_eqtl, full_eqtl, ref_name, min_overlap = 50) {
  message("Computing pi1 statistic for ", ref_name)
  
  merged <- left_join(query_eqtl, full_eqtl, by = c("variant_id", "phenotype_id"))
  valid_merged <- merged %>%
    filter(!is.na(pval), is.finite(pval), pval >= 0, pval <= 1)
  valid_pvals <- valid_merged$pval
  overlap_count <- length(valid_pvals)
  
  message("Number of overlapping eQTLs: ", overlap_count)
  
  if (overlap_count < min_overlap) {
    message("Too few overlapping eQTLs (<", min_overlap, ") for pi1 calculation.")
    return(list(pi1 = NA, prop_replicating = NA,
                prop_same_direction = NA,
                overlap_count = overlap_count,
                pvals = valid_pvals))
  }
  
  qobj <- tryCatch(
    qvalue(p = valid_pvals, pi0.method = "smoother"),
    error = function(e) {
      message("Smoother failed, trying bootstrap ...")
      tryCatch(
        qvalue(p = valid_pvals, pi0.method = "bootstrap"),
        error = function(e2) {
          message("Bootstrap also failed: ", e2$message)
          NULL
        }
      )
    }
  )
  
  if (!is.null(qobj)) {
    pi1 <- 1 - qobj$pi0
    prop_replicating <- mean(qobj$qvalues < 0.05)
    rep_idx <- which(qobj$qvalues < 0.05)
    replicating <- valid_merged[rep_idx, ]
    prop_same_direction <- if (nrow(replicating) > 0) {
      mean(sign(replicating$slope_my) == sign(replicating$slope_ref), na.rm = TRUE)
    } else {
      NA
    }
    message("pi1 calculated: ", round(pi1, 3))
  } else {
    pi1 <- NA
    prop_replicating <- NA
    prop_same_direction <- NA
  }
  
  return(list(pi1 = pi1,
              prop_replicating = prop_replicating,
              prop_same_direction = prop_same_direction,
              overlap_count = overlap_count,
              pvals = valid_pvals))
}

#------------------------------------------
# Main function
run_pi1_enrichment <- function(cell_type, public_all_qtl, public_top_qtl,
                               qtl_all, qtl_top, output_enrich, output_pi1,
                               ref_name, ref_cell_type = NA_character_) {
  
  # Skip internal self-comparisons but still write dummy outputs 
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
  
  # Load datasets ---------------------------------------------------------------
  if (str_detect(ref_name, 'obrien')) {
    message("Loading O'Brien bulk all eQTL ...")
    public_full <- read_tsv(public_all_qtl, show_col_types = FALSE) %>%
      select(variant_id, phenotype_id = gene_id,
             pval = pval_nominal, slope_ref = slope)
    
    message("Loading O'Brien bulk top eQTL (FDR < 0.05) ...")
    public_top <- read_tsv(public_top_qtl, show_col_types = FALSE) %>%
      filter(qval < 0.05) %>%
      select(variant_id, phenotype_id = gene_id, slope_my = slope)
  }
  
  if (str_detect(ref_name, 'wen')) {
    message("Loading Wen bulk all eQTL ...")
    public_full <- read_tsv(public_all_qtl, show_col_types = FALSE) %>%
      select(variant_id = sid, phenotype_id = pid,
             pval = npval, slope_ref = slope)
    
    message("Loading Wen bulk top eQTL (FDR < 0.05) ...")
    public_top <- read_tsv(public_top_qtl, show_col_types = FALSE) %>%
      filter(qval < 0.05) %>%
      select(variant_id = sid, phenotype_id = pid, slope_my = slope)
  }
  
  if (str_detect(ref_name, 'bryois')) {
    message("Loading Bryois nominal all-eQTL (no header) ...")
    public_nom <- read.table(public_all_qtl, header = T, sep = "", stringsAsFactors = FALSE) %>%
      as_tibble()
    public_full <- public_nom %>%
      select(variant_id, phenotype_id, pval = pval_nominal, slope_ref = slope) %>%
      mutate(phenotype_id = str_extract(phenotype_id, "ENSG[0-9]+"))
    
    message("Bryois nominal all-eQTL dimensions:", dim(public_full))
    message("Check that phenotype_ids are ENSG:\n")
    print(public_full[1:5, 1:4])
    
    # Map ref_cell_type (from Snakemake wildcard) to Excel cell_type string
    if (!(ref_cell_type %in% names(bryois_cell_map))) {
      stop("ref_cell_type '", ref_cell_type, "' not recognized. Must be one of: ",
           paste(names(bryois_cell_map), collapse = ", "))
    }
    excel_cell_type <- bryois_cell_map[[ref_cell_type]]
    
    message("\nLoading Bryois permutation / Table S2 (Excel) and extracting significant eQTLs for ref cell type: ", ref_cell_type)
    public_table <- read_excel(public_top_qtl, sheet = "Table S2", skip = 3)  # header start on row 4
    
    public_top <- public_table %>%
      filter(cell_type == excel_cell_type) %>%
      mutate(adj_p = suppressWarnings(as.numeric(adj_p))) %>%
      filter(!is.na(adj_p) & adj_p < 0.05) %>%
      transmute(variant_id = SNP,
                phenotype_id = ensembl,
                slope_my = as.numeric(beta))
    
    message("Bryois: number of top (FDR<0.05) eQTL retained for ", excel_cell_type, ": ", nrow(public_top))
  }
  
  if (str_detect(ref_name, "fugita")) {
    message("Loading Fugita snRNA-seq reference eQTL...")
    
    # read entire table once
    public_nom <- read_tsv(public_all_qtl, show_col_types = FALSE) %>% as_tibble()
    
    # Keep only rows for the requested fugita cell type, if provided
    if (!is.na(ref_cell_type)) {
      public_nom <- public_nom %>% filter(celltype == ref_cell_type)
      message("Filtering Fugita table to celltype: ", ref_cell_type, "  (rows retained: ", nrow(public_nom), ")")
    } else {
      message("No ref_cell_type wildcard provided; using entire Fugita table (rows: ", nrow(public_nom), ")")
    }
  
    public_full <- public_nom %>%
      transmute(
        variant_id = snps,
        phenotype_id = gene_id,
        pval = suppressWarnings(as.numeric(pvalue)),
        slope_ref = suppressWarnings(as.numeric(beta))
      ) %>%
      mutate(phenotype_id = str_extract(phenotype_id, "ENSG[0-9]+"))
    
    # top table = those marked as significant in two-step permutation (significant_by_2step_FDR == "Yes")
    public_top <- public_nom %>%
      filter(!is.na(significant_by_2step_FDR) & significant_by_2step_FDR == "Yes") %>%
      transmute(
        variant_id = snps,
        phenotype_id = gene_id,
        slope_my = suppressWarnings(as.numeric(beta))
      ) %>%
      mutate(phenotype_id = str_extract(phenotype_id, "ENSG[0-9]+"))
    
    message("Fugita: number of top (two-step Yes) eQTL retained: ", nrow(public_top))
    message("Fugita: number of nominal (full) rows retained: ", nrow(public_full))
  }

  if (str_detect(ref_name, "jang")) {
    message("Loading Jang et al. 2026 (SingleBrain) full-association eQTL...")

    # One file per cell type already (unlike Fugita's single combined file
    # with a `celltype` column) -- the right cell type is selected via the
    # jang_cell_type wildcard's file path, no in-script filtering needed.
    public_full_raw <- read_tsv(
      public_all_qtl,
      col_select = c(feature, variant_id, ref, alt, Allele, fixed_beta, Random_P, Fixed_P),
      show_col_types = FALSE
    )

    # ADDITION -- placeholder-cohort artifact filter (confirmed via manual
    # diagnosis): some cohorts encode "not tested for this pair" as
    # beta=0/sd=0 rather than NA. An SE of exactly 0 gives that cohort
    # infinite weight (1/SE^2) in the inverse-variance-weighted meta-
    # analysis, forcing the pooled beta to exactly 0 regardless of what
    # other cohorts show, and producing Random_P == Fixed_P == 1 exactly
    # (not just a large p-value -- an exact tie at the boundary). This is a
    # data-encoding artifact, not a genuine null result, and left in place
    # it creates a spike at p=1 that biases qvalue()'s pi0 (smoother)
    # estimate toward 1 regardless of real signal elsewhere in the
    # distribution. Excluded here on the signature Random_P==1 & Fixed_P==1
    # (both checked -- an artifact row shows this identical exact tie;
    # genuine strong nulls essentially never produce two different
    # meta-analysis models landing on the identical floating-point value).
    n_before_artifact_filter <- nrow(public_full_raw)
    public_full_raw <- public_full_raw %>%
      filter(!(Random_P == 1 & Fixed_P == 1))
    n_removed_artifact <- n_before_artifact_filter - nrow(public_full_raw)
    message("Jang full-association: removed ", n_removed_artifact,
            " placeholder-cohort artifact rows (Random_P==1 & Fixed_P==1 exactly) out of ",
            n_before_artifact_filter, " (",
            round(100 * n_removed_artifact / n_before_artifact_filter, 3), "%)")

    # Sign-normalize against Allele: Jang's beta is relative to Allele, not
    # always alt (see ADDITION note at top of file). This directly affects
    # prop_same_direction below, so skipping it would silently corrupt that
    # output for any row where Allele != alt.
    public_full <- public_full_raw %>%
      mutate(
        phenotype_id = str_extract(feature, "ENSG[0-9]+"),
        slope_ref = if_else(Allele == alt, fixed_beta, -fixed_beta),
        pval = suppressWarnings(as.numeric(Random_P))
      ) %>%
      select(variant_id, phenotype_id, pval, slope_ref)

    message("Loading Jang top-association (lead SNP) eQTL...")
    public_top_raw <- read_tsv(
      public_top_qtl,
      col_select = c(feature, variant_id, ref, alt, Allele, fixed_beta, Random_P, Fixed_P, qval),
      show_col_types = FALSE
    )

    # Same artifact filter applied defensively to the top file too. In
    # practice this should remove ~0 rows here, since qval is itself
    # derived from Random_P and an artifact row's Random_P==1 should already
    # fail the qval<0.05 filter below -- logged anyway so that assumption is
    # visible rather than silently relied upon.
    n_before_artifact_filter_top <- nrow(public_top_raw)
    public_top_raw <- public_top_raw %>%
      filter(!(Random_P == 1 & Fixed_P == 1))
    n_removed_artifact_top <- n_before_artifact_filter_top - nrow(public_top_raw)
    message("Jang top-association: removed ", n_removed_artifact_top,
            " placeholder-cohort artifact rows (Random_P==1 & Fixed_P==1 exactly) out of ",
            n_before_artifact_filter_top, " (expected ~0, since these should already fail qval<0.05)")

    # top table = gene's lead SNP with qval<0.05 (FDR across variants within
    # gene, then Storey's q across genes -- Jang's own two-stage correction)
    public_top <- public_top_raw %>%
      filter(!is.na(qval) & qval < 0.05) %>%
      mutate(
        phenotype_id = str_extract(feature, "ENSG[0-9]+"),
        slope_my = if_else(Allele == alt, fixed_beta, -fixed_beta)
      ) %>%
      select(variant_id, phenotype_id, slope_my)

    message("Jang: number of top (qval<0.05) eQTL retained: ", nrow(public_top))
    message("Jang: number of nominal (full) rows retained: ", nrow(public_full))
  }
  
  
  if (!(str_detect(ref_name, "wen") | str_detect(ref_name, "obrien") | str_detect(ref_name, "bryois") | str_detect(ref_name, "fugita") | str_detect(ref_name, "jang"))) {
    message("Internal mode: comparing ", cell_type, " ↔ ", ref_cell_type)
    
    # Load ref cell type all-eQTL
    public_full <- read_tsv(public_all_qtl, show_col_types = FALSE) %>%
      select(variant_id, phenotype_id, pval = pval_nominal, slope_ref = slope)
    
    # Load ref cell type top-eQTL
    public_top <- read_tsv(public_top_qtl, show_col_types = FALSE) %>%
      filter(qval < 0.05) %>%
      select(variant_id, phenotype_id, slope_my = slope)
  }
  
  # Load our dataset (query) -----------------------------------------------------
  message("Loading sig. eQTL (query) for cell type: ", cell_type)
  query_eqtl <- read_tsv(qtl_top, show_col_types = FALSE) %>%
    filter(qval < 0.05) %>%
    select(variant_id, phenotype_id, slope_my = slope)
  
  message("Loading all eQTL (query full) for cell type: ", cell_type)
  full_cell <- read_tsv(qtl_all, show_col_types = FALSE) %>%
    select(variant_id, phenotype_id,
           pval = pval_nominal, slope_ref = slope)
  
  # Initialize results tibble
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
  
  # Initialize list to store both pi1 results
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
  
  # Always write result (even if NA)
  write_rds(enrichment_results, output_enrich)
  write_rds(pi1_results, output_pi1)
  message("\nEnrichment tbl written to: ", output_enrich)
  message("\nPi1 results written to: ", output_pi1)
}
#------------------------------------------

# Run with Snakemake params
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
