#--------------------------------------------------------------------------------------
#
#    Enrichment of fetal-specific and shared eQTL in adult cortical OCRs
#
#--------------------------------------------------------------------------------------
#
# Pipeline:  EQTL-SHARING | Rule: ocr_enrich (single job)
#            Upstream:   classify_eqtl_specificity; heterogeneity_test; parse_ld_clump;
#                        07REPLICATION-SC snp_lookup; Li et al. 2023 cCREs (S6, S7)
#            Downstream: report (OCR enrichment tab)
#
# Purpose:   Test four LD-independent fetal SNP sets for overlap with the union of
#            Li et al. 2023 adult cortical cCREs (block bootstrap,
#            perform_permutation_test(): 1,000 iterations, 500 kb blocks):
#            - fetal_specific          testable, not significant in any adult cell type
#            - shared                  testable, significant in an adult cell type
#            - specific_heterogeneous  heterogeneity_test heterogeneous_snps
#            - not_expressed           gene not expressed in any Jang cell type
#
# Inputs:    classification_file    classify_eqtl_specificity output
#            heterogeneity_file     heterogeneity_test output
#            independent_snps_file  LD-independent pairs (parse_ld_clump)
#            snp_lookup_files       hg38 SNP coordinates, one file per cell type
#            li_s6_file             Li et al. 2023 cCRE coordinates
#            li_s7_file             Li et al. 2023 cCRE to cell type map
#            seed                   random seed for the bootstrap (params)
#
# Outputs:   output_file  rds: one row per SNP set (n_snps, n_peaks, p_value,
#                         observed_overlap, mean_permuted_overlap, fold_enrichment)
#                         (plus an untracked .tsv copy)
#
# Notes:     specific_heterogeneous is a subset of fetal_specific, so the two are
#            not independent tests.
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

message("\n\nTesting fetal-specific/shared/heterogeneous/not-expressed eQTL enrichment in adult cortical union OCRs ...")

##  Load packages, functions and variables  -------------------------------------------
suppressPackageStartupMessages({
  library(tidyverse)
  library(GenomicRanges)
  library(nullranges)
  library(GenomeInfoDb)
})

# Load functions
source(file.path(snakemake@scriptdir, "repl_sc_functions.R"))

# Input and output paths
classification_file    <- snakemake@input[["classification"]]
heterogeneity_file     <- snakemake@input[["heterogeneity_test"]]
independent_snps_file  <- snakemake@input[["independent_snps"]]
snp_lookup_files        <- snakemake@input[["snp_lookups"]]
li_s6_file             <- snakemake@input[["li_s6"]]
li_s7_file             <- snakemake@input[["li_s7"]]
output_file            <- snakemake@output[[1]]
seed                   <- snakemake@params[["seed"]]

# Print what each variable is set to
message("\nVariables")
message("============================")
tibble(
  variable = c("classification_file", "heterogeneity_file", "independent_snps_file",
               "n_snp_lookups", "li_s6_file", "li_s7_file", "output_file", "seed"),
  value    = c(classification_file, heterogeneity_file, independent_snps_file,
               length(snp_lookup_files), li_s6_file, li_s7_file, output_file, seed)) |>
  knitr::kable(format = "simple") |>
  cat(sep = "\n")
message("============================\n")

# Li et al. 2023 cortical cell types; cCREs in any of these form the union
all_cortical_types <- c(
  "ITL23_1", "ITL23_2", "ITL23_3", "ITL23_4", "ITL23_5", "ITL23_6",
  "ITL34", "ITL4_1", "ITL4_2", "ITL45_1", "ITL45_2",
  "ITV1C_1", "ITV1C_2", "ITV1C_3",
  "ITL5_1", "ITL5_2", "ITL5_3", "ITL5_4",
  "ITL6_1_1", "ITL6_1_2", "ITL6_2_1", "ITL6_2_2",
  "CT_1", "CT_2", "L6B_1", "L6B_2", "NP_1", "NP_2", "NP_3", "ET",
  "SST_1", "SST_2", "SST_3", "SST_4", "SST_5", "SST_CHODL",
  "PVALB_1", "PVALB_2", "PVALB_3", "PVALB_4", "PV_ChCs",
  "VIP_1", "VIP_2", "VIP_3", "VIP_4", "VIP_5", "VIP_6", "VIP_7",
  "SNCG_1", "SNCG_2", "SNCG_3", "SNCG_4", "SNCG_5",
  "LAMP5", "LAMP5_LHX6",
  "MGC_1", "MGC_2",
  "OPC", "COP",
  "EC", "SMC",
  "ASCT_1", "ASCT_2", "ASCT_3",
  "OGC_1", "OGC_2", "OGC_3"
)

##  Load SNP sets (LD-independent pairs only)  ----------------------------------------
message("\nLoading classification output ...")
classification <- read_rds(classification_file)
any_table <- classification$any_table

message("Loading heterogeneity test output ...")
het_test <- read_rds(heterogeneity_file)
heterogeneous_pairs <- het_test$heterogeneous_snps  # already heterogeneous_all AND magnitude_consistent

message("\nRestricting to LD-independent eGene-eSNP pairs ...")
independent_universe <- read_rds(independent_snps_file)

n_before_prune <- nrow(any_table)
any_table <- any_table %>%
  semi_join(independent_universe, by = c("phenotype_id", "variant_id"))
message("  any_table rows before LD pruning: ", n_before_prune)
message("  any_table rows after LD pruning: ", nrow(any_table))

n_het_before_prune <- nrow(heterogeneous_pairs)
heterogeneous_pairs <- heterogeneous_pairs %>%
  semi_join(independent_universe, by = c("phenotype_id", "variant_id"))
message("  heterogeneous_pairs rows before LD pruning: ", n_het_before_prune)
message("  heterogeneous_pairs rows after LD pruning: ", nrow(heterogeneous_pairs))

fetal_specific_snps <- any_table %>%
  filter(final_category == "testable_specific") %>%
  pull(variant_id) %>%
  unique()
message("  Fetal-specific unique SNPs: ", length(fetal_specific_snps))

shared_snps <- any_table %>%
  filter(final_category == "testable_shared") %>%
  pull(variant_id) %>%
  unique()
message("  Shared unique SNPs: ", length(shared_snps))

specific_heterogeneous_snps <- heterogeneous_pairs %>%
  pull(variant_id) %>%
  unique()
message("  Specific + heterogeneous unique SNPs: ", length(specific_heterogeneous_snps))

# not_expressed: gene not expressed in any Jang cell type, whether or not the
# SNP was genotyped
not_expressed_snps <- any_table %>%
  filter(final_category %in% c("untestable_genotyped_not_expressed",
                                "untestable_not_genotyped_not_expressed")) %>%
  pull(variant_id) %>%
  unique()
message("  Not-expressed unique SNPs: ", length(not_expressed_snps))

##  SNP coordinates  ------------------------------------------------------------------
message("\nBuilding SNP coordinate lookup from ", length(snp_lookup_files), " existing files ...")

snp_coords <- map_dfr(snp_lookup_files, read_rds) %>%
  distinct(refsnp_id, .keep_all = TRUE)
message("  Total unique SNPs with coordinates: ", nrow(snp_coords))

specific_gr               <- make_snp_gr(fetal_specific_snps, snp_coords)
shared_gr                 <- make_snp_gr(shared_snps, snp_coords)
specific_heterogeneous_gr <- make_snp_gr(specific_heterogeneous_snps, snp_coords)
not_expressed_gr          <- make_snp_gr(not_expressed_snps, snp_coords)

message("  Fetal-specific SNPs with coordinates: ", length(specific_gr))
message("  Shared SNPs with coordinates: ", length(shared_gr))
message("  Specific+heterogeneous SNPs with coordinates: ", length(specific_heterogeneous_gr))
message("  Not-expressed SNPs with coordinates: ", length(not_expressed_gr))

##  Li et al. 2023 adult cortical cCREs (union)  --------------------------------------
message("\nLoading Li et al. 2023 adult brain cCREs (cortical union only) ...")

options(scipen = 999)

li_s6 <- read_tsv(li_s6_file, col_names = c("chr", "start", "end", "ccre_id"),
                   col_types = "ciic", show_col_types = FALSE)
li_s7 <- read_tsv(li_s7_file, col_names = c("ccre_id", "cell_types"),
                   col_types = "cc", show_col_types = FALSE)

ccre_ids_union <- li_s7 %>%
  filter(map_lgl(cell_types, ~ any(str_split_1(.x, ",") %in% all_cortical_types))) %>%
  pull(ccre_id)

adult_cortical_union_gr <- li_s6 %>%
  filter(ccre_id %in% ccre_ids_union) %>%
  makeGRangesFromDataFrame(
    seqnames.field = "chr",
    start.field    = "start",
    end.field      = "end",
    keep.extra.columns = TRUE
  )

message("  adult_cortical_union: ", length(adult_cortical_union_gr), " peaks")

##  Enrichment tests  -----------------------------------------------------------------
message("\n--- Running enrichment tests (4 sets x adult_cortical_union) ---\n")

snp_sets <- list(
  fetal_specific         = specific_gr,
  shared                 = shared_gr,
  specific_heterogeneous = specific_heterogeneous_gr,
  not_expressed          = not_expressed_gr
)

results <- tibble()
set.seed(seed)

for (snp_set_name in names(snp_sets)) {
  snp_gr <- snp_sets[[snp_set_name]]

  message("\n=== adult_cortical_union vs ", snp_set_name, " SNPs (",
          length(snp_gr), " SNPs, ", length(adult_cortical_union_gr), " peaks) ===")

  perm_result <- perform_permutation_test(snp_gr, adult_cortical_union_gr, n_perm = 1000)

  results <- bind_rows(
    results,
    tibble(
      peak_set                = "adult_cortical_union",
      snp_set                  = snp_set_name,
      n_snps                    = length(snp_gr),
      n_peaks                    = length(adult_cortical_union_gr),
      p_value                      = perm_result$p_value,
      observed_overlap              = perm_result$observed_overlap,
      mean_permuted_overlap          = perm_result$mean_permuted_overlap,
      fold_enrichment                  = perm_result$fold_enrichment
    )
  )
}

##  Summary  --------------------------------------------------------------------------
message("\n--- Results summary ---\n")
print(results, n = Inf)

##  Write output  ---------------------------------------------------------------------
write_rds(results, output_file)
message("\nOCR enrichment results written to: ", output_file)

# Untracked TSV for readability outside R
output_tsv <- sub("\\.rds$", ".tsv", output_file)
write_tsv(results, output_tsv)
message("OCR enrichment results (TSV) written to: ", output_tsv)

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
