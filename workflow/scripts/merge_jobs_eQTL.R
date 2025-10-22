#--------------------------------------------------------------------------------------
#
#    Merge JOBS-boosted eQTL with original TensorQTL nom ouput per cell type
#
#--------------------------------------------------------------------------------------

# - 

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
message('\n\nMerging JOBS with TensorQTL nominals ...')

## Load libraries --------------------------------------------------------------------
library(data.table)
library(tidyverse)

# Inputs via snakemake
# Note: jobs_file is the jobs all file for smk rule order, need to change here for
# cell type specific data or you will get dups in ouput
tensor_file <- snakemake@input[['tensor']] 
jobs_file <- snakemake@input[['jobs']]
out_file <- snakemake@output[[1]]
cell_type <- snakemake@wildcards[['cell_type']] 
jobs_dir <- dirname(jobs_file)
jobs_file <- paste0(jobs_dir, "/jobs_", cell_type, "_beta_se_p_fdr.tsv.gz")

message(paste("Processing cell:", cell_type))
message(paste("Tensor file:", tensor_file))
message(paste("JOBS file:", jobs_file))
message(paste("Output:", out_file))

# Load TensorQTL nominals
message("Loading TensorQTL data ...")
tensor <- fread(tensor_file, header = TRUE, 
                colClasses = c(slope = "numeric", slope_se = "numeric", pval_nominal = "numeric", af = "numeric"))
message("Tensor loaded:", nrow(tensor), "rows")
message("Zeros in slope (beta): ", sum(tensor$slope == 0, na.rm = TRUE))
message("Zeros in slope_se (SE): ", sum(tensor$slope_se == 0, na.rm = TRUE))
message("Zeros in pval_nominal: ", sum(tensor$pval_nominal == 0, na.rm = TRUE))

# Create ID for merging
tensor[, ID := paste(phenotype_id, variant_id, sep = "-")]

# Load JOBS cell-specific
message("Loading JOBS data ...")
jobs <- fread(jobs_file, header = TRUE,
              colClasses = c(beta = "numeric", se = "numeric", pval_nominal = "numeric", fdr = "numeric"))
message("JOBS loaded:", nrow(jobs), "rows")
message("Zeros in slope (beta): ", sum(jobs$beta == 0, na.rm = TRUE))
message("Zeros in slope_se (SE): ", sum(jobs$se == 0, na.rm = TRUE))
message("Zeros in pval_nominal: ", sum(jobs$pval_nominal == 0, na.rm = TRUE))

# Ensure ID in jobs (fallback to gene-snp if needed)
if (!"ID" %in% names(jobs)) {
  jobs[, ID := paste(gene, snp, sep = "-")]
}
message("Unique JOBS IDs:", uniqueN(jobs$ID))

# Left merge: tensor + jobs on ID
message("Merging on ID ...")
merged <- merge(tensor, jobs[, .(ID, jobs_beta = beta, jobs_se = se, jobs_pval = pval_nominal, fdr)], 
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

message(paste("Merged has:", nrow(merged), "rows"))
message("Zeros in slope (beta): ", sum(merged$slope == 0, na.rm = TRUE))
message("Zeros in slope_se (SE): ", sum(merged$slope_se == 0, na.rm = TRUE))
message("Zeros in pval_nominal: ", sum(merged$pval_nominal == 0, na.rm = TRUE))

# Write output
message("Writing hybrid output ...")
fwrite(merged, out_file, sep = "\t", na = "NA", quote = FALSE, row.names = FALSE)
message(paste("Wrote", nrow(merged), "rows to", out_file))

# Quick QC
n_boosted <- sum(!is.na(merged$beta))
message(paste("Pairs with JOBS boost:", n_boosted, "(", round(100 * n_boosted / nrow(merged), 1), "%)"))

message("Done!")

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------