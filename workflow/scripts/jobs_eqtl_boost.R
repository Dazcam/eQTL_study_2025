#--------------------------------------------------------------------------------------
#
#    Running JOBS to boost sc-eQTL signal
#
#--------------------------------------------------------------------------------------

# Paper: https://www.cell.com/cell-genomics/fulltext/S2666-979X(25)00076-X
# GitHub: https://github.com/LidaWangPSU/JOBS

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
message('\n\nRunning JOBS to boost sc-eQTL signal  ...')

## Load libraries and variables -------------------------------------------------------
library(JOBS)
library(data.table)
library(tidyverse)

bulk_file <- snakemake@input[[1]]
sc_dir <- snakemake@params[['sc_dir']]
out_dir <- snakemake@params[['out_dir']]
cell_types <- c("ExN-UL", "ExN-DL", "RG", "InN", "Endo-Peri", "OPC", "MG")
dups_dir <- paste0(out_dir, "duplicates/")  
dir.create(dups_dir, recursive = TRUE)

# Load full bulk, add SE and ID
message("Loading bulk data and adding SEs ...")
bulk <- fread(bulk_file, header = TRUE, colClasses = c(npval = "numeric", slope = "numeric"))
bulk[, se := abs(slope) / qnorm(pmax(1e-10, 1 - npval / 2))]  # Avoid p=1 issues
bulk[, ID := paste(pid, sid, sep = "-")]

# Dedup bulk 
n_dups_bulk <- anyDuplicated(bulk, by = "ID")
if (n_dups_bulk > 0) {
  message(paste("Found", n_dups_bulk, "duplicates in bulk, deduping..."))
  bulk <- unique(bulk, by = "ID")
}
message(paste("Bulk loaded:", nrow(bulk), "unique pairs"))

# Prepare wide sc beta/se by stepwise outer merge (memory-efficient)
message("\nMerging sc data wide...")
sc_beta_wide <- data.table(ID = character(0))  # Start with empty
sc_se_wide <- data.table(ID = character(0))

for (cell in cell_types) {
  
  message("\nLoading sc data for ", cell ," ...")
  sc_file <- paste0(sc_dir, cell, "/", cell, "_nom.cis_qtl_pairs.tsv")
  
  if (!file.exists(sc_file)) { 
    warning(paste("Missing file:", sc_file)); next 
  }
  
  message("Creating unique eQTL IDs and pulling out betas and SEs for ", cell ," ...")
  sc <- fread(sc_file, header = TRUE, 
              colClasses = c(slope = "numeric", slope_se = "numeric"))
  sc[, ID := paste(phenotype_id, variant_id, sep = "-")]
  
  # Dedup sc data (to avoid merge errors)
  n_dups_sc <- anyDuplicated(sc, by = "ID")
  if (n_dups_sc > 0) {
    message(paste("Found", n_dups_sc, "duplicates in", cell, "- deduping to", nrow(unique(sc, by = "ID")), "rows"))
    
    # Write dups table to check  
    dups <- sc[duplicated(sc, by = "ID"), ]
    fwrite(dups, paste0(dups_dir, cell, "_dups.tsv.gz"), sep = "\t", na = "NA")
    message(paste("Exported", nrow(dups), "dup rows to", dups_dir, cell, "_dups.tsv.gz"))
    
    # Peek at first 10 dup IDs
    message(paste("Sample dup IDs:", paste(head(dups$ID, 10), collapse = ", ")))

    sc <- unique(sc, by = "ID")  # Keeps first occurrence
  }
  
  sc_beta_cell <- sc[, .(ID, beta = slope)]
  sc_se_cell <- sc[, .(ID, se = slope_se)]
  
  # Outer merge to growing wide
  sc_beta_wide <- merge(sc_beta_wide, sc_beta_cell, by = "ID", all = TRUE)
  sc_se_wide <- merge(sc_se_wide, sc_se_cell, by = "ID", all = TRUE)
  colnames(sc_beta_wide)[ncol(sc_beta_wide)] <- cell
  colnames(sc_se_wide)[ncol(sc_se_wide)] <- cell
  
  # Post-merge dedup (double-check)
  n_dups_wide <- anyDuplicated(sc_beta_wide, by = "ID")
  if (n_dups_wide > 0) {
    message(paste("Post-merge dups in wide for", cell, "- deduping"))
    sc_beta_wide <- unique(sc_beta_wide, by = "ID")
    sc_se_wide <- unique(sc_se_wide, by = "ID")
  }
  
  message(paste("Merged", cell, "- rows:", nrow(sc_beta_wide)))
  rm(sc, sc_beta_cell, sc_se_cell); gc()  # Cleanup per cell
}

# Add bulk: outer merge
message("\nMerging with bulk...")
beta_all <- merge(bulk[, .(ID, bulk = slope)], sc_beta_wide, by = "ID", all = TRUE)
se_all <- merge(bulk[, .(ID, bulk = se)], sc_se_wide, by = "ID", all = TRUE)

# Reorder: ID, bulk, cells (drop cols for missing cells if any)
avail_cells <- setdiff(cell_types, c("bulk"))  # Assume all loaded
beta_all <- beta_all[, c("ID", "bulk", avail_cells), with = FALSE]
se_all <- se_all[, c("ID", "bulk", avail_cells), with = FALSE]

# Drop pairs with all sc NA (JOBs skips them anyway, but reduces size)
sc_mask <- rowSums(is.na(beta_all[, avail_cells, with = FALSE])) < length(avail_cells)
beta_all <- beta_all[sc_mask]
se_all <- se_all[sc_mask]
message(paste("Final pairs after cleaning:", nrow(beta_all)))

if (nrow(beta_all) == 0) stop("No overlapping data!")

# Estimate weights with optional subsampling (recommended for large data)
subsample_n <- 0  # Adjust: 0 for full (risky), 50000-200000 stable
message(paste("\nEstimating weights (subsample:", ifelse(subsample_n > 0, subsample_n, "full"), ")..."))
if (subsample_n > 0 && nrow(beta_all) > subsample_n) {
  sub_idx <- sample(nrow(beta_all), subsample_n)
  sub_beta <- beta_all[sub_idx]
  sub_se <- se_all[sub_idx]
} else {
  sub_beta <- beta_all
  sub_se <- se_all
}
weight <- jobs.nnls.weights(sub_beta, sub_se)
message(paste("Weights:", paste(round(weight, 3), collapse = ", ")))

# Refine eQTLs on full data
message("\nRefining eQTLs...")
jobs_eqtls <- jobs.eqtls(beta_all, se_all, weight, COR = FALSE)
ref_beta <- jobs_eqtls$eqtls_new
ref_se <- jobs_eqtls$eqtls_se_new

head(jobs_eqtls$eqtls_new)
head(jobs_eqtls$eqtls_se_new)

# Write uncorrected output files
message("\nSaving JOBS output uncorrected files ...")
write_rds(jobs_eqtls, paste0(out_dir, "jobs_eqtlsput.rds"))
fwrite(ref_beta, paste0(out_dir, "jobs_ref_beta_genomewide.tsv.gz"), sep = "\t", na = "NA")
fwrite(ref_se, paste0(out_dir, "jobs_ref_se_genomewide.tsv.gz"), sep = "\t", na = "NA")


# Ensure ID is character
ref_beta$ID <- as.character(ref_beta$ID)
ref_se$ID <- as.character(ref_se$ID)

# Coerce all cell columns to numeric
numeric_cols <- c("bulk", avail_cells)
ref_beta[numeric_cols] <- lapply(ref_beta[numeric_cols], as.numeric)
ref_se[numeric_cols] <- lapply(ref_se[numeric_cols], as.numeric)

# Check structure
message("Structure of ref_beta:")
print(str(ref_beta))
if (!"ID" %in% colnames(ref_beta)) stop("ID column missing in ref_beta - check JOBs output")

# Compute nominal p-values and per-gene FDR with dplyr (cleaner)
message("\nComputing p-values and per-gene FDR ...")
# Add gene column once
ref_beta$gene <- sub("-.*", "", ref_beta$ID)

pval_dt <- data.frame()  # Start empty
n_cells <- length(avail_cells)
for (i in seq_along(avail_cells)) {
  cell <- avail_cells[i]
  message(paste("Processing cell", i, "/", n_cells, ":", cell))
  
  if (!cell %in% colnames(ref_beta)) { 
    message(paste("Skipping", cell, "- col missing"))
    next 
  }
  
  # Compute p_nom
  beta_col <- ref_beta[[cell]]
  se_col <- ref_se[[cell]]
  p_nom <- 2 * pnorm(-abs(beta_col / se_col))
  
  # Skip if all NA
  if (all(is.na(p_nom))) {
    message(paste("Skipping FDR for", cell, "- all pvals NA"))
    next
  }
  
  # Add p to ref_beta (mutate style)
  ref_beta[[paste0(cell, "_p")]] <- p_nom
  
  # FDR per gene: group and adjust
  ref_beta_temp <- ref_beta %>%
    group_by(gene) %>%
    mutate(fdr = p.adjust(.data[[paste0(cell, "_p")]], method = "BH")) %>%
    ungroup() %>%
    select(ID, gene, fdr)
  
  # Add cell label and bind
  fdr_dt <- ref_beta_temp
  fdr_dt$cell <- cell
  pval_dt <- rbind(pval_dt, fdr_dt)
  
  # Warn on NAs
  if (any(is.na(fdr_dt$fdr))) {
    message(paste("Warning: NAs in FDR for", cell, "- check small gene groups"))
  }
}

# Convert pval_dt back to data.table for fwrite
setDT(pval_dt)

# Export
message("\nSaving JOBS FDR corrected output files ...")
fwrite(pval_dt, paste0(out_dir, "jobs_pval_fdr_genomewide.tsv.gz"), sep = "\t", na = "NA")

message("Done! Check outputs in", out_dir)
message("For significant hits: e.g., sig_hits <- pval_dt[fdr < 0.05]")

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------