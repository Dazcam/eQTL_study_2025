#--------------------------------------------------------------------------------------
#
#    Enrichment of eQTL in ATAC-seq peaks
#
#--------------------------------------------------------------------------------------

# See bootRanges for SNP-Peak overlap: https://academic.oup.com/bioinformatics/article/39/5/btad190/7115835

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

message("\n\nRunning enrichment analysis for cell-specific sig. eQTL in Ziffra snATAC-seq peaks ...")

# -------------------------------------------------------------------------------------

library(tidyverse)
library(qvalue)
library(GenomicRanges)
library(nullranges)
library(GenomeInfoDb)
library(readxl)

# Load functions
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

# Define directories and cell types (from your script)
#cell_types <- c("ExN-UL", "ExN-DL", "InN", "RG", "MG", "OPC", "Endo-Peri")
# exp_pc <- 40
# gen_pc <- 4
# in_dir <- "~/Desktop/eQTL_study_2025/workflow/reports/05TENSORQTL/logs/"
# public_dir <- "~/Desktop/eQTL_study_2025/resources/public_datasets/"
# tensorqtl_dir <- paste0(in_dir, 'results/05TENSORQTL/')
# in_dir <- paste0(public_dir, 'ziffra_2021/')

# Make a tibble showing what each variable is set to
message("\nVariables")
cat("============================")
tibble(
  variable = c("qtl_perm", "snp_file", "out_file", "peaks", "peak_dir", "cell_type", 
               "gen_pc", "exp_pc", "norm_method"),
  value    = c(qtl_perm, snp_file, out_file, peaks, peak_dir, cell_type, 
               gen_pc, exp_pc, norm_method)) |> 
  knitr::kable(format = "simple", align = "l") |>
  print()
message("\n============================\n")

##  Load Ziffra peaks  -----------------------------------------------------------------
# peak list = all peaks key; mac2 all peaks per cell type; specific; specific peaks per cell type
cat('\nLoading Ziffra 2021 snATAC-seq fetal brain data ... \n')
options(scipen = 999) # required to prevent peak coords. being abbr. in sci' notation
peak_list <- read_excel(file.path(peak_dir, "Ziffra_2021_supp_tables_2_13.xlsx"), sheet = "ST2 AllPrimaryPeaks") %>%
  dplyr::select(seqnames, start, end, peak_name)
#enhancer_peaks <- read_excel(file.path(peak_dir, "Ziffra_2021_supp_tables_2_13.xlsx"), sheet = "ST5 EnhancerPeaks") %>%
#  dplyr::select(peak_name) # Might use these later
lookup_peaks <- read_excel(file.path(peak_dir, "Ziffra_2021_supp_tables_2_13.xlsx"), sheet = "ST3 MACSpeaks_byCelltype")

cell_types <- c("Glu-UL", "Glu-DL", "GABA", "RG", "MG", "OPC", "Endo-Peri",
                "Glu-UL-0", "Glu-UL-1", "Glu-UL-2", 
                "Glu-DL-0", "Glu-DL-1", "Glu-DL-2",
                "GABA-0", "GABA-1", "GABA-2", 
                "NPC-0", "NPC-1", "NPC-2")
ziffra_mapping_all <- list(
  "NPC"        = "RG_MACSpeaks",
  "NPC-0"      = "RG_MACSpeaks", 
  "NPC-1"      = "RG_MACSpeaks", 
  "NPC-2"      = "RG_MACSpeaks", 
  "Glu-UL"    = "ulEN_MACSpeaks",
  "Glu-UL-0"  = "ulEN_MACSpeaks", 
  "Glu-UL-1"  = "ulEN_MACSpeaks", 
  "Glu-UL-2"  = "ulEN_MACSpeaks",
  "Glu-DL"    = "dlEN_MACSpeaks",
  "Glu-DL-0"  = "dlEN_MACSpeaks",
  "Glu-DL-1"  = "dlEN_MACSpeaks",
  "Glu-DL-2"  = "dlEN_MACSpeaks",
  "GABA"       = c("IN_CGE_MACSpeaks", "IN_MGE_MACSpeaks"),
  "GABA-0"     = c("IN_CGE_MACSpeaks", "IN_MGE_MACSpeaks"),
  "GABA-1"     = c("IN_CGE_MACSpeaks", "IN_MGE_MACSpeaks"),
  "GABA-2"     = c("IN_CGE_MACSpeaks", "IN_MGE_MACSpeaks"),
  "Endo-Peri" = "EndoMural_MACSpeaks",
  "MG"        = "Microglia_MACSpeaks",
  "OPC"       = "AstroOligo_MACSpeaks"
)

# Only run relevant cell type mappings
ziffra_mapping <- ziffra_mapping_all[cell_type]

### Ziffra peak overlaps
ziffra_cell_types <- names(lookup_peaks)[-1]  # Remove 'peak_name'

# Create a list where each element is a filtered lookup table
results_list <- lapply(ziffra_cell_types, function(cell_type) {
  peaks_with_1 <- lookup_peaks %>%
    filter(!!sym(cell_type) == 1) %>%  # select peaks where 1
    pull(peak_name)
  
  peak_list  %>%
    filter(peak_name %in% peaks_with_1) %>%
    mutate(Cell_Type = cell_type)  # Optional: keep track of which cell type
})

# Bind into one table if you want:
ziffra_peaks <- bind_rows(results_list)

# Map eQTL cell types to Ziffra cell types
cell_type_mapping <- ziffra_mapping

# Results
ziffra_enrichment_results <- tibble()

message("\n--- Starting Ziffra ATAC-seq Peak Enrichment Analysis ---")
for (cell_type in names(cell_type_mapping)) {
  message("\nProcessing cell type: ", cell_type)
  
  # Load significant eQTL SNPs
  my_eqtl <- read_tsv(qtl_perm) %>%
    filter(qval < 0.05) %>%
    dplyr::select(variant_id)
  
  if (nrow(my_eqtl) == 0) {
    message("No significant eQTLs found for ", cell_type)
    next
  }
  
  # Convert to GRanges
  sig_snps_gr <- make_snp_gr(my_eqtl$variant_id, read_rds(snp_file))
  
  # Match peaks
  matching_ziffra_cells <- cell_type_mapping[[cell_type]]
  
  for (ziffra_cell in matching_ziffra_cells) {
    message("\nPerforming permutation test for ", cell_type, " vs ", ziffra_cell)
    
    peaks_gr <- ziffra_peaks |> 
      filter(Cell_Type == ziffra_cell) |>
      makeGRangesFromDataFrame(
        keep.extra.columns = TRUE
      )
    
    message('Check unique cell types in peaks: ', unique(peaks_gr$Cell_Type))
    
    # Permutation test
    perm_pvalue_ziffra <- perform_permutation_test(sig_snps_gr, peaks_gr, n_perm = 1000)
    
    # Record result
    ziffra_enrichment_results <- bind_rows(
      ziffra_enrichment_results,
      tibble(
        eqtl_cell_type = cell_type,
        peak_cell_type = ziffra_cell,
        p_value = perm_pvalue_ziffra$p_value,
        observed_overlap = perm_pvalue_ziffra$observed_overlap,
        mean_permuted_overlap = perm_pvalue_ziffra$mean_permuted_overlap,
        fold_enrichment = perm_pvalue_ziffra$fold_enrichment
        # (perm_overlaps is NOT included here because it's a vector, not a scalar)
      )
    )
  }
}
message("\n---------------------------------------------------------")

message("\nSaving enrichment results to: ", out_file)
write_rds(ziffra_enrichment_results, out_file)

message("Done.")

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
