#--------------------------------------------------------------------------------------
#
#    Generate eQTL supplementary table for manuscript
#
#--------------------------------------------------------------------------------------

# For each cell type extract sig. eQTL
# Add chr, position and allele info from genotype VCF
# Check for OCR overlap of each eQTL in Ziffra union of OCRs file
# Note: Input tbls were generated in smr_report.Rmd
#
# NOTE: the two lines above (VCF chr/position/allele annotation, OCR overlap
# check) describe steps that are NOT implemented below -- flagging in case
# this is a stale docstring rather than the intended current behaviour.

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

smr_dir <- snakemake@params[["smr_dir"]]
eqtl_nom_dir <- snakemake@params[["eqtl_nom_dir"]]
cell_types <- snakemake@params[["cell_types"]]
disorders <- snakemake@params[["disorders"]]
out_file <- snakemake@output[[1]]

# smr_dir <- "../results/13MANUSCRIPT_TABLES/"
# eqtl_nom_dir <- "../results/10SMR/smr_input/"
# cell_types <- c("Glu-UL", "Glu-DL", "NPC", "GABA", "Endo-Peri", "OPC", "MG", ...)
# disorders <- c("scz", "bpd", "mdd", "adhd", "ocd")

message("Cell types (", length(cell_types), "): ", paste(cell_types, collapse = ", "))

# Relabel for display only: Glu-UL -> Glu-A, Glu-DL -> Glu-B (covers L2
# subtypes too, e.g. Glu-UL-0 -> Glu-A-0). File paths below still use the
# raw cell_type value, since that's what's actually on disk.
relabel_cell_type <- function(x) {
  x <- str_replace(x, "Glu-UL", "Glu-A")
  x <- str_replace(x, "Glu-DL", "Glu-B")
  x
}
message("Disorders (", length(disorders), "): ", paste(disorders, collapse = ", "))

# ── Collect all significant SMR associations ─────────────────────────────────────
message("Collecting significant SMR eQTLs across disorders...")

# Forces these 6 columns to a consistent type on every read, regardless of
# how many data rows a given {gwas}_smr.tsv has. Without this, read_tsv()
# infers types per file independently: a populated file with a purely
# numeric Chr column infers double, while a file with 0 data rows (nothing
# passed the significance filter for that GWAS, e.g. ADHD) has nothing to
# guess a type from and can default to character instead - bind_rows()
# below then fails combining a double Chr against a character Chr.
smr_col_types <- cols(
  Chr = col_character(), Symbol = col_character(), `Ensembl ID` = col_character(),
  SNP = col_character(), A1 = col_character(), A2 = col_character(),
  .default = col_guess()
)

sig_smr_list <- list()
for (gwas in disorders) {
  
  message('Extracting sig. SMR eQTL for: ', gwas)
  smr_tbl <- read_tsv(paste0(smr_dir, gwas, '_smr.tsv'),
                      show_col_types = FALSE, col_types = smr_col_types) |>
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
  
  display_name <- relabel_cell_type(cell_type)
  
  eqtl_tbl <- read_tsv(paste0(eqtl_nom_dir, cell_type, '/',
                              cell_type, '_nom.cis_qtl_pairs.tsv'),
                       show_col_types = FALSE) |>
    mutate(key = paste(variant_id, phenotype_id, sep = '_')) |>
    filter(key %in% sig_smr_tbl$key) |>
    select(key, !!display_name := slope)
    
  
  message('No. eQTL found for ', cell_type, ': ', nrow(eqtl_tbl))
  
  sig_smr_tbl <- sig_smr_tbl |>
    left_join(eqtl_tbl, by = "key")
  
  rm(eqtl_tbl)
  
  }

# How many cell types have the eQTL for each row
# sig_smr_tbl <- sig_smr_tbl |>
#   mutate(
#     n_cell_types_with_eqtl = rowSums(!is.na(select(., all_of(cell_types))))
#   )

# ── Export ───────────────────────────────────────────────────────────────────────
message("Writing to Excel file...")
write.xlsx(sig_smr_tbl,
           file = out_file,
           overwrite = TRUE,
           headerStyle = createStyle(textDecoration = "bold"))

message("Export Complete.")

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
