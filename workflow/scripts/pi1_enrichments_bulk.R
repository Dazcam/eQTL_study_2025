#--------------------------------------------------------------------------------------
#
#    pi1 enrichments between sig. sn-eQTL & O'brien 2018 bulk brain eQTL
#   
#--------------------------------------------------------------------------------------

# Questions: What is the lower threshold for overlaps to calc pi1?
# MG and OPC currently failing as too few SNPs to judge pi1

# pi1 process
# For each cell type:
#   1. Extract sig. eQTL @ FDR < 0.05
#   2. Extract pvals in ref eQTL data for SNPs ID'ed in step 1
#   3. Compute pi1 stat using the pvals from step 2 (i.e. ref data)

## Info  ------------------------------------------------------------------------------

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

message("\n\nCalculating pi1 enrichments: sig. sn-eQTL & public bulk brain eQTL ...")

# -------------------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(tidyverse)
  library(qvalue)
  library(cowplot)
})

# Input and output paths
public_all_qtl <- snakemake@input[['public_all']]
public_top_qtl <- snakemake@input[['public_top']]
qtl_all <- snakemake@input[['qtl_all']]
qtl_top <- snakemake@input[['qtl_top']]
output_enrich <- snakemake@output[["enrich"]]
output_pi1 <- snakemake@output[["pi1"]]
cell_type <- snakemake@wildcards[["cell_type"]]
exp_pc <- snakemake@wildcards[["exp_pc"]]
geno_pc <- snakemake@wildcards[["geno_pc"]]
norm_method <- snakemake@wildcards[["norm_method"]]
ref_name <- str_split_i(output_pi1, "/", 4)

# Check variable assignment
message("\nVariables")
cat("============================")
tibble(
  variable = c("public_all_qtl", "public_top_qtl", "qtl_all", "qtl_top", "output_enrich", 
               "output_pi1", "cell_type", "exp_pc", "geno_pc", "norm_method", "ref_name"),
  value    = c(public_all_qtl, public_top_qtl, qtl_all, qtl_top, output_enrich, 
               output_pi1, cell_type, exp_pc, geno_pc, norm_method, ref_name)
) |> 
  knitr::kable(format = "simple", align = "l") |>
  print()
message("\n============================\n")


# Functions  --------------------------------------------------------------------------
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
                               ref_name) {
  
  # Load datasets
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
  
  message("Loading sig. eQTL: ", cell_type)
  query_eqtl <- read_tsv(qtl_top, show_col_types = FALSE) %>%
    filter(qval < 0.05) %>%
    select(variant_id, phenotype_id, slope_my = slope)
  
  message("Loading all eQTL: ", cell_type)
  full_cell <- read_tsv(qtl_all, show_col_types = FALSE) %>%
    select(variant_id, phenotype_id,
           pval = pval_nominal, slope_ref = slope)
  
  # Initialize results tibble
  enrichment_results <- tibble(
    cell_type = cell_type,
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
  
  # Forward: cell sig > O’Brien full
  if (nrow(query_eqtl) > 0) {
    message("\n--- Forward enrichment: cell sig → ", ref_name, " full ---")
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
  
  # Reverse: O’Brien sig > cell full
  if (nrow(full_cell) > 0) {
    message("\n--- Reverse enrichment: ", ref_name, " sig → cell full ---")
    pi1_result <- compute_pi1(public_top, full_cell, cell_type)
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
    ref_name        = ref_name
  )
}
#------------------------------------------

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------



