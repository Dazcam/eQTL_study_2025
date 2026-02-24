#--------------------------------------------------------------------------------------
#
#    Generate eQTL supplementary table for manuscript
#
#--------------------------------------------------------------------------------------

# For each cell type extract sig. eQTL
# Add chr, position and allele info from genotype VCF
# Check for OCR overlap of each eQTL in Ziffra union of OCRs file
# Note: Input tbls were generated in smr_report.Rmd

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

message("\n\neQTL supplementary table for manuscript ...")

# -------------------------------------------------------------------------------------

library(tidyverse)
library(openxlsx)
g
smr_dir <- snakemake@params[["smr_dir"]]
eqtl_nom_dir <- snakemake@params[["eqtl_nom_dir"]]
out_file <- snakemake@output[[1]]

cell_types <- c("Glu-UL", "Glu-DL", "NPC", "GABA",
                "Endo-Peri", "OPC", "MG",
                "Glu-UL-0", "Glu-UL-1", "Glu-UL-2",
                "Glu-DL-0", "Glu-DL-1", "Glu-DL-2",
                "GABA-0", "GABA-1", "GABA-2",
                "NPC-0", "NPC-1", "NPC-2")

disorders <- c('scz', 'bpd', 'mdd', 'ocd') # No sig. SMR eQTL for ADHD

# ── Collect all significant SMR associations ─────────────────────────────────────
message("Collecting significant SMR eQTLs across disorders...")

sig_smr_list <- list()
for (gwas in disorders) {
  
  message('Extracting sig. SMR eQTL for: ', gwas)
  smr_tbl <- read_tsv(paste0(smr_dir, gwas, '_smk.tsv'),
                      show_col_types = FALSE) |>
    select(Chr, Symbol, `Ensembl ID`, SNP, A1, A2) |>
    mutate(key = paste(SNP, `Ensembl ID`, sep = '_'),
           GWAS = gwas)
  
  sig_smr_list[[gwas]] <- smr_tbl
  
}

# Bind rows - May want to omit duplicates here
sig_smr_tbl <- bind_rows(sig_smr_list)
message('Total sig. SMR: ', nrow(sig_smr_tbl))

# ── Add nominal eQTL slopes from every cell type ────────────────────────────────
message("Extracting nominal eQTL betas for significant pairs...")
for (cell_type in cell_types) {
  
  eqtl_tbl <- read_tsv(paste0(eqtl_nom_dir, cell_type, '/', 
                              cell_type, '_nom.cis_qtl_pairs.tsv'),
                       show_col_types = FALSE) |>
    mutate(key = paste(variant_id, phenotype_id, sep = '_')) |>
    filter(key %in% sig_smr_tbl$key) |>
    select(key, !!cell_type := slope)
    
  
  message('No. eQTL found for ', cell_type, ': ', nrow(eqtl_tbl))
  
  sig_smr_tbl <- sig_smr_tbl |>
    left_join(eqtl_tbl, by = "key")
  
  rm(eqtl_tbl)
  
  }

# How many cell types have the eQTL for each row
sig_smr_tbl <- sig_smr_tbl |>
  mutate(
    n_cell_types_with_eqtl = rowSums(!is.na(select(., all_of(cell_types))))
  )

# ── Export ───────────────────────────────────────────────────────────────────────
message("Writing to Excel file...")
write.xlsx(sig_smr_tbl, 
           file = out_file, 
           overwrite = TRUE,
           headerStyle = createStyle(textDecoration = "bold"))

message("Export Complete.")

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------