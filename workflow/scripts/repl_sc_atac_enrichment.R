#--------------------------------------------------------------------------------------
#
#    Enrichment of sig. eQTL in Ziffra 2021 snATAC-seq peaks
#
#--------------------------------------------------------------------------------------
#
# Pipeline:  06QTL-REPLICATION | Rule: atac_enrich (one job per cell type)
#            Upstream:   05TENSORQTL tensorqtl_perm; snp_lookup; dwnld_ziffra
#            Downstream: 06 report (Ziffra tab); 13 eqtl_tbl / 14 eqtl_qc_plt (Fig 2E)
#
# Purpose:   - Take the cell type's sig. eQTL (qval < 0.05) and their hg38
#              coordinates from the SNP lookup
#            - Map the cell type to its prenatal snATAC-seq peak set(s) (Ziffra
#              2021, ST3 MACS peaks by cell type; GABA types tested against both
#              CGE and MGE interneuron peaks)
#            - Test overlap enrichment against a block-bootstrap null
#              (perform_permutation_test(): 1,000 iterations, 500 kb blocks)
#
# Inputs:    qtl_perm   tensorQTL permutation output for the cell type
#            snp_file   SNP lookup (hg38 coordinates of the sig. eQTL)
#            peaks      Ziffra 2021 supplementary tables (xlsx)
#            peak_dir   directory holding the Ziffra xlsx
#
# Outputs:   out_file   enrichment table: one row per peak set tested
#
# Notes:     bootRanges: https://academic.oup.com/bioinformatics/article/39/5/btad190/7115835
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

message("\n\nRunning enrichment analysis for cell-specific sig. eQTL in Ziffra snATAC-seq peaks ...")

##  Load packages, functions and variables  -------------------------------------------
library(tidyverse)
library(GenomicRanges)
library(readxl)

source(file.path(snakemake@scriptdir, "repl_sc_functions.R"))

# Input and output paths
qtl_perm <- snakemake@input[["qtl_perm"]]
snp_file <- snakemake@input[["snp_file"]]
peaks <- snakemake@input[["peaks"]]
out_file <- as.character(snakemake@output)
peak_dir <- snakemake@params[["peak_dir"]]
cell_type <- snakemake@wildcards[["cell_type"]]
gen_pc <- snakemake@wildcards[["geno_pc"]]
exp_pc <- snakemake@wildcards[["exp_pc"]]
norm_method <- snakemake@wildcards[["norm_method"]]

# Our cell type -> Ziffra 2021 ST3 peak set(s)
ziffra_mapping_all <- list(
  "NPC"       = "RG_MACSpeaks",
  "NPC-0"     = "RG_MACSpeaks",
  "NPC-1"     = "RG_MACSpeaks",
  "NPC-2"     = "RG_MACSpeaks",
  "Glu-UL"    = "ulEN_MACSpeaks",
  "Glu-UL-0"  = "ulEN_MACSpeaks",
  "Glu-UL-1"  = "ulEN_MACSpeaks",
  "Glu-UL-2"  = "ulEN_MACSpeaks",
  "Glu-DL"    = "dlEN_MACSpeaks",
  "Glu-DL-0"  = "dlEN_MACSpeaks",
  "Glu-DL-1"  = "dlEN_MACSpeaks",
  "Glu-DL-2"  = "dlEN_MACSpeaks",
  "GABA"      = c("IN_CGE_MACSpeaks", "IN_MGE_MACSpeaks"),
  "GABA-0"    = c("IN_CGE_MACSpeaks", "IN_MGE_MACSpeaks"),
  "GABA-1"    = c("IN_CGE_MACSpeaks", "IN_MGE_MACSpeaks"),
  "GABA-2"    = c("IN_CGE_MACSpeaks", "IN_MGE_MACSpeaks"),
  "Endo-Peri" = "EndoMural_MACSpeaks",
  "MG"        = "Microglia_MACSpeaks",
  "OPC"       = "AstroOligo_MACSpeaks"
)
ziffra_sets <- ziffra_mapping_all[[cell_type]]

# Print what each variable is set to
message("\nVariables")
message("============================")
tibble(
  variable = c("qtl_perm", "snp_file", "out_file", "peaks", "peak_dir", "cell_type",
               "gen_pc", "exp_pc", "norm_method", "ziffra_sets"),
  value    = c(qtl_perm, snp_file, out_file, peaks, peak_dir, cell_type,
               gen_pc, exp_pc, norm_method, paste(ziffra_sets, collapse = ", "))) |>
  knitr::kable(format = "simple") |>
  cat(sep = "\n")
message("============================\n")

##  Load sig. eQTL  ---------------------------------------------------------------------
my_eqtl <- read_tsv(qtl_perm) %>%
  filter(qval < 0.05) %>%
  dplyr::select(variant_id)

if (nrow(my_eqtl) == 0) {
  message("No significant eQTLs found for ", cell_type, "; writing empty table.")
  write_rds(tibble(), out_file)
  quit(save = "no", status = 0)
}

sig_snps_gr <- make_snp_gr(my_eqtl$variant_id, read_rds(snp_file))

##  Load Ziffra peaks  -----------------------------------------------------------------
# ST2: coordinates of all peaks; ST3: 0/1 membership of each peak per cell type
message("\nLoading Ziffra 2021 snATAC-seq fetal brain peaks ...")
options(scipen = 999) # stop peak coordinates being written in scientific notation
peak_list <- read_excel(file.path(peak_dir, "Ziffra_2021_supp_tables_2_13.xlsx"), sheet = "ST2 AllPrimaryPeaks") %>%
  dplyr::select(seqnames, start, end, peak_name)
lookup_peaks <- read_excel(file.path(peak_dir, "Ziffra_2021_supp_tables_2_13.xlsx"), sheet = "ST3 MACSpeaks_byCelltype")

ziffra_peaks <- map(ziffra_sets, function(ziffra_ct) {
  peaks_in_set <- lookup_peaks %>%
    filter(!!sym(ziffra_ct) == 1) %>%
    pull(peak_name)

  peak_list %>%
    filter(peak_name %in% peaks_in_set) %>%
    mutate(Cell_Type = ziffra_ct)
}) |>
  bind_rows()

##  Enrichment tests  ------------------------------------------------------------------
message("\n--- Starting Ziffra ATAC-seq Peak Enrichment Analysis ---")
ziffra_enrichment_results <- tibble()

for (ziffra_cell in ziffra_sets) {
  message("\nPerforming permutation test for ", cell_type, " vs ", ziffra_cell)

  peaks_gr <- ziffra_peaks |>
    filter(Cell_Type == ziffra_cell) |>
    makeGRangesFromDataFrame(keep.extra.columns = TRUE)

  perm_result <- perform_permutation_test(sig_snps_gr, peaks_gr, n_perm = 1000)

  # perm_overlaps (vector of null overlaps) is not kept in the output table
  ziffra_enrichment_results <- bind_rows(
    ziffra_enrichment_results,
    tibble(
      eqtl_cell_type = cell_type,
      peak_cell_type = ziffra_cell,
      p_value = perm_result$p_value,
      observed_overlap = perm_result$observed_overlap,
      mean_permuted_overlap = perm_result$mean_permuted_overlap,
      fold_enrichment = perm_result$fold_enrichment
    )
  )
}
message("\n---------------------------------------------------------")

##  Write output  ----------------------------------------------------------------------
message("\nSaving enrichment results to: ", out_file)
write_rds(ziffra_enrichment_results, out_file)

message("Done.")

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
