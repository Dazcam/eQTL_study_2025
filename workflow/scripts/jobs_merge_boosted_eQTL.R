#--------------------------------------------------------------------------------------
#
#    MERGE JOBS-BOOSTED eQTLs WITH TENSORQTL NOMINAL RESULTS (PER CELL TYPE)
#
#--------------------------------------------------------------------------------------

#--------------------------------------------------------------------------------------
#
# 1. Combine original TensorQTL nominal eQTLs with JOBS-boosted effect estimates
# 2. Retain all TensorQTL associations (left join)
# 3. Replace beta, SE, and p-value with JOBS values where available
# 4. Fixes JOBS SE = 0 artifact via gene-specific mean SE imputation (excluding zeros in calc)
# 5. Recalculates ominal P values from β/SE
# 6. Compute final per-gene FDR across *all* hybrid associations
# 6. Output clean hybrid table (no FDR) + separate FDR-augmented file
#
#
# INPUTS (via Snakemake):
# - tensor_file:  TensorQTL nominal output
# - jobs_file:    JOBS output: jobs_{cell_type}_beta_se_p_fdr.tsv.gz
#
# OUTPUT:
# - Hybrid eQTL table: original + JOBS boosted, with fixed SEs and recalculated p-values
#
#--------------------------------------------------------------------------------------

## Set up logging for smk  ------------------------------------------------------------
if (exists("snakemake")) {
  log_smk <- function() {
    if (exists("snakemake") & length(snakemake@log) != 0) {
      log <- file(snakemake@log[1][[1]], open = "wt")
      sink(log, append = TRUE)
      sink(log, append = TRUE, type = "message")
    }
  }
}
log_smk()
message('\n\nMerging JOBS boosted eQTL with TensorQTL nominal results ...')

## Load libraries --------------------------------------------------------------------
library(data.table)
library(tidyverse)

# Inputs via snakemake
tensor_file <- snakemake@input[['tensor']] 
jobs_dir <- snakemake@params[['jobs_dir']]
out_file <- snakemake@output[[1]]
cell_type <- snakemake@wildcards[['cell_type']] 
jobs_file <- paste0(jobs_dir, "/jobs_", cell_type, "_beta_se_p.tsv.gz")


message(paste("\nProcessing cell:", cell_type))
message(paste("Tensor file:", tensor_file))
message(paste("JOBS file:", jobs_file))
message(paste("Output:", out_file))

# Load TensorQTL nominals
message("\nLoading TensorQTL data for ", cell_type, "...")
tensor <- fread(tensor_file, header = TRUE, 
                colClasses = c(slope = "numeric", slope_se = "numeric", pval_nominal = "numeric", af = "numeric"))
message("Tensor loaded:", nrow(tensor), " rows")
message("Zeros in slope (beta): ", sum(tensor$slope == 0, na.rm = TRUE))
message("Zeros in slope_se (SE): ", sum(tensor$slope_se == 0, na.rm = TRUE))
message("Zeros in pval_nominal: ", sum(tensor$pval_nominal == 0, na.rm = TRUE))

# Create ID for merging
tensor[, ID := paste(phenotype_id, variant_id, sep = "-")]

# Load JOBS cell-specific
message("\nLoading JOBS for ", cell_type, "...")
jobs <- fread(jobs_file, header = TRUE,
              colClasses = c(beta = "numeric", se = "numeric", pval_nominal = "numeric", fdr = "numeric"))
message("JOBS loaded:", nrow(jobs), " rows")
message("Zeros in slope (beta): ", sum(jobs$beta == 0, na.rm = TRUE))
message("Zeros in slope_se (SE): ", sum(jobs$se == 0, na.rm = TRUE))
message("Zeros in pval_nominal: ", sum(jobs$pval_nominal == 0, na.rm = TRUE))

# Ensure ID in jobs (fallback to gene-snp if needed)
if (!"ID" %in% names(jobs)) {
  jobs[, ID := paste(gene, snp, sep = "-")]
}
message("Unique JOBS IDs:", uniqueN(jobs$ID))

## Impute SEs for genes with SE values of zero in JOBS --------------------------------
message("Imputing SE = 0 with gene-specific mean SE and recalculating pval/FDR...")

# Compute mean SE *per gene* (excluding zeros and NAs)
gene_mean_se <- jobs[se > 0 & !is.na(se), .(mean_se = mean(se, na.rm = TRUE)), by = gene]

# Identify rows where se == 0
zero_se_mask <- jobs$se == 0 | is.na(jobs$se)

if (any(zero_se_mask)) {
  n_zero <- sum(zero_se_mask)
  message(paste(" Found", n_zero, "rows with SE == 0 or NA. Imputing with gene mean SE..."))
  
  # Merge mean SE back
  jobs_imputed <- merge(jobs, gene_mean_se, by = "gene", all.x = TRUE)
  
  # Replace se == 0 or NA with mean_se (only where applicable)
  jobs_imputed[zero_se_mask, se := ifelse(is.na(mean_se), se, mean_se)]
  
  # Clean up
  jobs_imputed[, mean_se := NULL]
  
  # Overwrite jobs
  jobs <- jobs_imputed
  rm(jobs_imputed, gene_mean_se)
} else {
  message(" No SE == 0 found. Skipping imputation.")
}

# Recalculate pval_nominal from beta / se
message(" Recalculating pval_nominal from beta / se ...")
jobs[, pval_nominal := 2 * pnorm(-abs(beta / se))]

# Handle any remaining Inf/-Inf or NA p-values (e.g., beta=0, se=0 still)
jobs[is.infinite(pval_nominal) | is.na(pval_nominal), pval_nominal := 1]

message(paste(" After imputation: Zeros in SE =", sum(jobs$se == 0, na.rm = TRUE)))
message(paste(" Min SE =", round(min(jobs$se[jobs$se > 0], na.rm = TRUE), 6)))
message(paste(" Max pval_nominal =", round(max(jobs$pval_nominal, na.rm = TRUE), 6)))

# Left merge: tensor + jobs on ID
message("\nMerging on ID ...")
merged <- merge(tensor, jobs[, .(ID, jobs_beta = beta, jobs_se = se, jobs_pval = pval_nominal)], 
                by = "ID", all.x = TRUE)
message(paste("Merged:", nrow(merged), "rows"))

# Hybrid update: replace with JOBS where non-NA
message("Applying hybrid updates ...")
merged[, slope := ifelse(!is.na(jobs_beta), jobs_beta, slope)]
merged[, slope_se := ifelse(!is.na(jobs_se), jobs_se, slope_se)]
merged[, pval_nominal := ifelse(!is.na(jobs_pval), jobs_pval, pval_nominal)]

# Drop temp columns
merged[, c("ID", "jobs_beta", "jobs_se", "jobs_pval", "fdr") := NULL]

# Reorder columns to match original
setcolorder(merged, c("phenotype_id", "variant_id", "slope", "slope_se", "pval_nominal", "af"))

message(paste("\nMerged has:", nrow(merged), "rows"))
message("Zeros in slope (beta): ", sum(merged$slope == 0, na.rm = TRUE))
message("Zeros in slope_se (SE): ", sum(merged$slope_se == 0, na.rm = TRUE))
message("Zeros in pval_nominal: ", sum(merged$pval_nominal == 0, na.rm = TRUE))

# Write output
message("\nWriting hybrid output ...")
fwrite(merged, out_file, sep = "\t", na = "NA", quote = FALSE, row.names = FALSE)
message(paste("Wrote", nrow(merged), "rows to", out_file))

message("Computing per-gene FDR across all hybrid eQTLs...")
merged[, gene := sub("^(ENSG[0-9]+).*", "\\1", phenotype_id)]
merged[, fdr := p.adjust(pval_nominal, method = "BH"), by = gene]

# Write FDR-augmented version
fdr_file <- sub("\\.tsv$", "_fdr.tsv", out_file)
fwrite(merged, fdr_file, sep = "\t", na = "NA", quote = FALSE, row.names = FALSE)
message(paste("Wrote FDR-augmented file:", fdr_file))
message(paste("  → Significant eQTLs (FDR < 0.05):", sum(merged$fdr < 0.05, na.rm = TRUE)))

# Clean up temporary column
merged[, gene := NULL]

message("All Done!")

message("All Done!")

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------