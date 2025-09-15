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

# Input and output paths
qtl_perm <- snakemake@input[["qtl_perm"]]
snp_file <- snakemake@input[["snp_file"]]
out_file <- snakemake@output[["out_file"]]
peak_dir <- snakemake@params[["peak_dir"]]
cell_type <- snakemake@wildcards[["cell_type"]]
gen_pc <- snakemake@wildcards[["gen_pc"]]
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
  variable = c("qtl_perm", "snp_file", "out_file", "peak_dir", "cell_type", 
               "gen_pc", "exp_pc", "norm_method"),
  value    = c(qtl_perm, snp_file, out_file, peak_dir, cell_type, 
               gen_pc, exp_pc, norm_method)) |> 
  knitr::kable(format = "simple", align = "l") |>
  print()
message("\n============================\n")

##  Load Ziffra peaks  -----------------------------------------------------------------
# peak list = all peaks key; mac2 all peaks per cell type; specific; specific peaks per cell type
cat('\nLoading Ziffra 2021 snATAC-seq fetal brain data ... \n')
options(scipen = 999) # required to prevent peak coords. being abbr. in sci' notation
peak_list <- read_excel(paste0(peak_dir, 'Ziffra_2021_supp_tables_2_13.xlsx'), sheet = 'ST2 AllPrimaryPeaks') %>%
  dplyr::select(seqnames, start, end, peak_name) 
enhancer_peaks <- read_excel(paste0(peak_dir, 'Ziffra_2021_supp_tables_2_13.xlsx'), sheet = 'ST5 EnhancerPeaks') %>%
  dplyr::select(peak_name) 
lookup_peaks <- read_excel(paste0(peak_dir, 'Ziffra_2021_supp_tables_2_13.xlsx'), sheet = 'ST3 MACSpeaks_byCelltype')
cell_types <- c("ExN-UL", "ExN-DL", "InN", "RG", "MG", "OPC", "Endo-Peri")
ziffra_mapping <- list(
  "RG" = "RG_MACSpeaks",
  "ExN-UL" = "dlEN_MACSpeaks",
  "ExN-DL" = "dlEN_MACSpeaks",
  "InN" = c("IN_CGE_MACSpeaks", "IN_MGE_MACSpeaks"),
  "Endo-Peri" = "EndoMural_MACSpeaks",
  "MG" = "Microglia_MACSpeaks",
  "OPC" = "AstroOligo_MACSpeaks"
)

# Load your SNP coordinates (assuming snp_lookup_hg38_filt is available from your script)
# For simplicity, assume it's pre-loaded or regenerated as needed

# Function for permutation test with ATAC-seq peaks
perform_permutation_test <- function(sig_snps_gr, peaks_gr, n_perm = 1000) {
  message("Starting permutation test for ATAC-seq peak enrichment...")
  
  # Verify bootRanges availability
  if (!exists("bootRanges", envir = asNamespace("nullranges"))) {
    stop("bootRanges function not found. Ensure nullranges is installed and loaded.")
  }
  
  # Ensure peaks_gr has proper seqinfo for hg38
  hg38_info <- getChromInfoFromUCSC("hg38")
  seqinfo(peaks_gr) <- Seqinfo(seqnames = hg38_info$chrom,
                               seqlengths = hg38_info$size,
                               genome = "hg38")[seqlevels(peaks_gr)]
  
  observed_overlap <- length(findOverlaps(sig_snps_gr, peaks_gr))
  message("Observed overlaps: ", observed_overlap)
  
  # Set block length to preserve local structure (e.g., 500kb)
  block_length <- 500000
  
  perm_overlaps <- numeric(n_perm)
  for (i in seq_len(n_perm)) {
    if (i %% 100 == 0) message("Permutation ", i, " of ", n_perm)
    # Generate bootstrapped peaks
    boot_peaks <- nullranges::bootRanges(
      y = peaks_gr,
      blockLength = block_length,
      R = 1,
      withinChrom = TRUE
    )
    perm_overlaps[i] <- length(findOverlaps(sig_snps_gr, boot_peaks))
  }
  
  # Calculate corrected p-value
  p_value <- (sum(perm_overlaps >= observed_overlap) + 1) / (n_perm + 1)
  
  # Calculate mean of permuted overlaps
  mean_permuted_overlap <- mean(perm_overlaps)
  
  # Calculate fold enrichment
  fold_enrichment <- observed_overlap / mean_permuted_overlap
  
  message("Permutation p-value: ", signif(p_value, 3))
  message("Fold enrichment: ", signif(fold_enrichment, 3))
  
  return(list(
    p_value = p_value,
    observed_overlap = observed_overlap,
    mean_permuted_overlap = mean_permuted_overlap,
    fold_enrichment = fold_enrichment,
    perm_overlaps = perm_overlaps  # include raw permuted values for plotting if needed
  ))
}


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
  sig_snps_gr <- read_rds(snp_file) %>%
    filter(refsnp_id %in% my_eqtl$variant_id) %>%
    distinct(refsnp_id, .keep_all = TRUE) %>%
    makeGRangesFromDataFrame(
      seqnames.field = "chr_name",
      start.field   = "chrom_start",
      end.field     = "chrom_end",
      keep.extra.columns = TRUE
    )
  
  # Match peaks
  matching_ziffra_cells <- cell_type_mapping[[cell_type]]
  
  for (ziffra_cell in matching_ziffra_cells) {
    message("Performing permutation test for ", cell_type, " vs ", ziffra_cell)
    
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

write_rds(ziffra_enrichment_results, out_file)

# Move plots to Rmd file
# ziffra_enrichment_results <- ziffra_enrichment_results %>%
#   mutate(
#     logp = -log10(p_value),
#     sig = ifelse(p_value < 0.05, "*", "")
#   ) |>
#   complete(eqtl_cell_type, peak_cell_type) |>
#   ggplot(aes(x = peak_cell_type, y = eqtl_cell_type, fill = fold_enrichment)) +
#   geom_tile(color = "black", linewidth = 0.3, na.rm = FALSE) +
#   geom_text(aes(label = ifelse(!is.na(p_value) & p_value < 0.05, "*", "")),
#             vjust = 0.5, size = 4, color = "black") +
#   scale_fill_viridis_c(
#     option = "C",
#     direction = 1,
#     name = "Enrichment",
#     na.value = "white"   # blanks appear white
#   )+
#   labs(
#     x = "ATAC-seq peak cell type",
#     y = "eQTL cell type",
#     title = "eQTL enrichment in ATAC-seq peaks"
#   ) +
#   theme_minimal(base_size = 14) +
#   theme(
#     axis.text.x = element_text(angle = 45, hjust = 1),
#     panel.grid = element_blank()
#   )
# 
# 
# ziffra_enrichment_results <- ziffra_enrichment_results %>%
#   mutate(
#     logp = -log10(p_value),
#     sig = ifelse(p_value < 0.05, "*", "")
#   ) |> ggplot(aes(x = eqtl_cell_type, y = fold_enrichment, fill = peak_cell_type)) +
#   geom_col(position = position_dodge(width = 0.8), width = 0.7, color = "black") +
#   geom_text(aes(label = sig),
#             position = position_dodge(width = 0.8),
#             vjust = -0.5, size = 5, color = "black") +
#   labs(
#     title = "eQTL enrichment in ATAC-seq peaks",
#     x = "eQTL cell type",
#     y = "Fold enrichment",
#     fill = "ATAC-seq peak cell type"
#   ) +
#   scale_fill_brewer(palette = "Set2") +
#   theme_classic(base_size = 14) +
#   theme(
#     axis.text.x = element_text(angle = 45, hjust = 1),
#     legend.position = "right"
#   )

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------