#--------------------------------------------------------------------------------------
#
#    Developmental specificity: OCR enrichment of fetal-specific eQTL
#
#--------------------------------------------------------------------------------------

#  Tests whether fetal-specific eQTL SNPs (from the "any" tier of the
#  classification step) are enriched in:
#    (a) fetal open chromatin regions (Ziffra et al. 2021 snATAC-seq)
#    (b) adult open chromatin regions (Li et al. 2023 snATAC-seq)
#
#  Uses the same bootRanges block-bootstrap null-range approach as the
#  existing replication_atac_enrichments.R, with identical parameters
#  (500kb block length, 1000 permutations, within-chromosome).
#
#  The prediction: fetal-specific SNPs are enriched in fetal OCRs but
#  NOT enriched (or significantly less enriched) in adult OCRs.
#
#  Adult peaks from Li et al. 2023 are grouped into L1-equivalent sets
#  matching our fetal cell populations, plus a cortical union set.
#
#  Step 3 of the developmental-specificity pipeline.
#
#--------------------------------------------------------------------------------------
#
#  ADDITION -- LD-independence filtering (reviewer request): before pulling
#  fetal-specific/shared SNP lists, any_table is restricted to the LD-independent
#  (eGene, eSNP) universe produced by the 3-step LD-clumping sub-pipeline
#  (devspec_prep_ld_clump_input.R -> run_ld_clump -> devspec_parse_ld_clump.R). This
#  ensures genes with multiple mapped eSNPs across cell types contribute only their
#  LD-independent representative SNP(s) to the OCR enrichment SNP sets, rather than
#  redundant SNPs tagging the same underlying signal. See parse_ld_clump's log for
#  the final independent-pair counts (3,327 pairs / 3,121 genes at time of writing).
#
#--------------------------------------------------------------------------------------

# Set up logging for Snakemake
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

message("\n\nTesting fetal-specific eQTL enrichment in fetal and adult OCRs ...")

suppressPackageStartupMessages({
  library(tidyverse)
  library(GenomicRanges)
  library(nullranges)
  library(GenomeInfoDb)
  library(readxl)
})

# Input and output paths ---------------------------------------------------------------
classification_file  <- snakemake@input[["classification"]]
independent_snps_file <- snakemake@input[["independent_snps"]]
snp_lookup_files     <- snakemake@input[["snp_lookups"]]
li_s6_file           <- snakemake@input[["li_s6"]]
li_s7_file           <- snakemake@input[["li_s7"]]
output_file          <- snakemake@output[[1]]
ziffra_dir           <- snakemake@params[["ziffra_dir"]]

# Check variable assignment
message("\nVariables")
cat("============================")
tibble(
  variable = c("classification_file", "independent_snps_file", "n_snp_lookups",
               "li_s6_file", "li_s7_file", "ziffra_dir", "output_file"),
  value    = c(classification_file, independent_snps_file, length(snp_lookup_files),
               li_s6_file, li_s7_file, ziffra_dir, output_file)
) |>
  knitr::kable(format = "simple", align = "l") |>
  print()
message("\n============================\n")

# Li et al. 2023 cell type labels → fetal L1 population mapping.
# Cortical neuronal + glial types only. Excludes cerebellar, thalamic,
# midbrain, basal ganglia populations. Hardcoded — not a tunable parameter.
li_to_L1 <- list(
  "Glu-UL" = c("ITL23_1", "ITL23_2", "ITL23_3", "ITL23_4", "ITL23_5", "ITL23_6",
               "ITL34", "ITL4_1", "ITL4_2", "ITL45_1", "ITL45_2",
               "ITV1C_1", "ITV1C_2", "ITV1C_3"),
  "Glu-DL" = c("ITL5_1", "ITL5_2", "ITL5_3", "ITL5_4",
               "ITL6_1_1", "ITL6_1_2", "ITL6_2_1", "ITL6_2_2",
               "CT_1", "CT_2", "L6B_1", "L6B_2", "NP_1", "NP_2", "NP_3", "ET"),
  "GABA"   = c("SST_1", "SST_2", "SST_3", "SST_4", "SST_5", "SST_CHODL",
               "PVALB_1", "PVALB_2", "PVALB_3", "PVALB_4", "PV_ChCs",
               "VIP_1", "VIP_2", "VIP_3", "VIP_4", "VIP_5", "VIP_6", "VIP_7",
               "SNCG_1", "SNCG_2", "SNCG_3", "SNCG_4", "SNCG_5",
               "LAMP5", "LAMP5_LHX6"),
  "MG"     = c("MGC_1", "MGC_2"),
  "OPC"    = c("OPC", "COP"),
  "Endo-Peri" = c("EC", "SMC")
)

# Additional cortical types for the union set (not matched to any fetal L1)
cortical_glial_extra <- c("ASCT_1", "ASCT_2", "ASCT_3",
                          "OGC_1", "OGC_2", "OGC_3")

all_cortical_types <- c(unlist(li_to_L1, use.names = FALSE), cortical_glial_extra)

# Ziffra fetal peak mapping — test each fetal cell type's peaks individually,
# plus a union of all fetal peaks. Uses all Ziffra peak cell types.
ziffra_peak_types <- c("RG_MACSpeaks", "ulEN_MACSpeaks", "dlEN_MACSpeaks",
                       "IN_CGE_MACSpeaks", "IN_MGE_MACSpeaks",
                       "Microglia_MACSpeaks", "AstroOligo_MACSpeaks",
                       "EndoMural_MACSpeaks")

# Functions ----------------------------------------------------------------------------

# Identical to replication_atac_enrichments.R's perform_permutation_test
perform_permutation_test <- function(sig_snps_gr, peaks_gr, n_perm = 1000) {
  message("Starting permutation test for ATAC-seq peak enrichment...")

  if (!exists("bootRanges", envir = asNamespace("nullranges"))) {
    stop("bootRanges function not found. Ensure nullranges is installed and loaded.")
  }

  # Get hg38 chromosome sizes — try UCSC first, fall back to hardcoded
  hg38_seqinfo <- tryCatch({
    message("Fetching hg38 seqinfo from UCSC ...")
    hg38_info <- getChromInfoFromUCSC("hg38")
    Seqinfo(seqnames = hg38_info$chrom,
            seqlengths = hg38_info$size,
            genome = "hg38")
  }, error = function(e) {
    message("UCSC fetch failed (", e$message, "), using hardcoded hg38 chromosome sizes.")
    hg38_seqlengths <- c(
      chr1 = 248956422, chr2 = 242193529, chr3 = 198295559, chr4 = 190214555,
      chr5 = 181538259, chr6 = 170805979, chr7 = 159345973, chr8 = 145138636,
      chr9 = 138394717, chr10 = 133797422, chr11 = 135086622, chr12 = 133275309,
      chr13 = 114364328, chr14 = 107043718, chr15 = 101991189, chr16 = 90338345,
      chr17 = 83257441, chr18 = 80373285, chr19 = 58617616, chr20 = 64444167,
      chr21 = 46709983, chr22 = 50818468, chrX = 156040895, chrY = 57227415
    )
    Seqinfo(seqnames = names(hg38_seqlengths),
            seqlengths = hg38_seqlengths,
            genome = "hg38")
  })

  seqinfo(peaks_gr) <- hg38_seqinfo[seqlevels(peaks_gr)]

  observed_overlap <- length(findOverlaps(sig_snps_gr, peaks_gr))
  message("Observed overlaps: ", observed_overlap)

  block_length <- 500000

  perm_overlaps <- numeric(n_perm)
  for (i in seq_len(n_perm)) {
    if (i %% 100 == 0) message("Permutation ", i, " of ", n_perm)
    boot_peaks <- nullranges::bootRanges(
      y = peaks_gr,
      blockLength = block_length,
      R = 1,
      withinChrom = TRUE
    )
    perm_overlaps[i] <- length(findOverlaps(sig_snps_gr, boot_peaks))
  }

  p_value <- (sum(perm_overlaps >= observed_overlap) + 1) / (n_perm + 1)
  mean_permuted_overlap <- mean(perm_overlaps)
  fold_enrichment <- observed_overlap / mean_permuted_overlap

  message("Permutation p-value: ", signif(p_value, 3))
  message("Fold enrichment: ", signif(fold_enrichment, 3))

  return(list(
    p_value = p_value,
    observed_overlap = observed_overlap,
    mean_permuted_overlap = mean_permuted_overlap,
    fold_enrichment = fold_enrichment,
    perm_overlaps = perm_overlaps
  ))
}

# Main -----------------------------------------------------------------------------------

# 1. Load fetal-specific SNPs from classification output, restricted to the
#    LD-independent universe -----------------------------------------------------------
message("\nLoading fetal-specific SNP list from classification output ...")
classification <- read_rds(classification_file)
any_table <- classification$any_table

message("\nRestricting to LD-independent eGene-eSNP pairs ...")
independent_universe <- read_rds(independent_snps_file)

n_before_prune <- nrow(any_table)
any_table <- any_table %>%
  semi_join(independent_universe, by = c("phenotype_id", "variant_id"))
message("  any_table rows before LD pruning: ", n_before_prune)
message("  any_table rows after LD pruning: ", nrow(any_table))

fetal_specific_snps <- any_table %>%
  filter(final_category == "testable_specific") %>%
  pull(variant_id) %>%
  unique()

message("  Fetal-specific unique SNPs: ", length(fetal_specific_snps))

# Also get shared SNPs for a contrast test
shared_snps <- any_table %>%
  filter(final_category == "testable_shared") %>%
  pull(variant_id) %>%
  unique()

message("  Shared unique SNPs (for contrast): ", length(shared_snps))

# 2. Build SNP coordinate lookup from existing per-cell-type files ----------------------
message("\nBuilding SNP coordinate lookup from ", length(snp_lookup_files), " existing files ...")

snp_coords <- map_dfr(snp_lookup_files, read_rds) %>%
  distinct(refsnp_id, .keep_all = TRUE)

message("  Total unique SNPs with coordinates: ", nrow(snp_coords))

# Create GRanges for fetal-specific and shared SNPs
make_snp_gr <- function(snp_ids, coords) {
  coords %>%
    filter(refsnp_id %in% snp_ids) %>%
    distinct(refsnp_id, .keep_all = TRUE) %>%
    makeGRangesFromDataFrame(
      seqnames.field = "chr_name",
      start.field    = "chrom_start",
      end.field      = "chrom_end",
      keep.extra.columns = TRUE
    )
}

specific_gr <- make_snp_gr(fetal_specific_snps, snp_coords)
shared_gr   <- make_snp_gr(shared_snps, snp_coords)

message("  Fetal-specific SNPs with coordinates: ", length(specific_gr))
message("  Shared SNPs with coordinates: ", length(shared_gr))

n_missing_specific <- length(fetal_specific_snps) - length(specific_gr)
n_missing_shared   <- length(shared_snps) - length(shared_gr)
if (n_missing_specific > 0) message("  WARNING: ", n_missing_specific, " fetal-specific SNPs missing from lookup files")
if (n_missing_shared > 0)   message("  WARNING: ", n_missing_shared, " shared SNPs missing from lookup files")

# 3. Load Li et al. 2023 adult peaks ----------------------------------------------------
message("\nLoading Li et al. 2023 adult brain cCREs ...")

options(scipen = 999)

li_s6 <- read_tsv(li_s6_file, col_names = c("chr", "start", "end", "ccre_id"),
                   col_types = "ciic", show_col_types = FALSE)
message("  Total cCREs in S6: ", nrow(li_s6))

li_s7 <- read_tsv(li_s7_file, col_names = c("ccre_id", "cell_types"),
                   col_types = "cc", show_col_types = FALSE)
message("  Total cCREs in S7: ", nrow(li_s7))

# Build per-L1 and cortical-union peak GRanges
build_li_peaks_gr <- function(li_s6, li_s7, cell_type_labels) {
  ccre_ids <- li_s7 %>%
    filter(map_lgl(cell_types, ~ any(str_split_1(.x, ",") %in% cell_type_labels))) %>%
    pull(ccre_id)

  li_s6 %>%
    filter(ccre_id %in% ccre_ids) %>%
    makeGRangesFromDataFrame(
      seqnames.field = "chr",
      start.field    = "start",
      end.field      = "end",
      keep.extra.columns = TRUE
    )
}

message("\nConstructing adult peak sets ...")

adult_peak_sets <- list()

for (l1_name in names(li_to_L1)) {
  adult_peak_sets[[paste0("adult_", l1_name)]] <- build_li_peaks_gr(li_s6, li_s7, li_to_L1[[l1_name]])
  message("  adult_", l1_name, ": ", length(adult_peak_sets[[paste0("adult_", l1_name)]]), " peaks")
}

adult_peak_sets[["adult_cortical_union"]] <- build_li_peaks_gr(li_s6, li_s7, all_cortical_types)
message("  adult_cortical_union: ", length(adult_peak_sets[["adult_cortical_union"]]), " peaks")

# 4. Load Ziffra 2021 fetal peaks -------------------------------------------------------
message("\nLoading Ziffra 2021 fetal brain snATAC-seq peaks ...")

peak_list <- read_excel(file.path(ziffra_dir, "Ziffra_2021_supp_tables_2_13.xlsx"),
                        sheet = "ST2 AllPrimaryPeaks") %>%
  dplyr::select(seqnames, start, end, peak_name)

lookup_peaks <- read_excel(file.path(ziffra_dir, "Ziffra_2021_supp_tables_2_13.xlsx"),
                           sheet = "ST3 MACSpeaks_byCelltype")

ziffra_cell_types_available <- names(lookup_peaks)[-1]  # drop peak_name column

fetal_peak_sets <- list()

for (zct in ziffra_peak_types) {
  if (!(zct %in% ziffra_cell_types_available)) {
    message("  WARNING: Ziffra cell type '", zct, "' not found in lookup table, skipping.")
    next
  }

  peak_ids <- lookup_peaks %>%
    filter(!!sym(zct) == 1) %>%
    pull(peak_name)

  fetal_peak_sets[[paste0("fetal_", zct)]] <- peak_list %>%
    filter(peak_name %in% peak_ids) %>%
    makeGRangesFromDataFrame(keep.extra.columns = TRUE)

  message("  fetal_", zct, ": ", length(fetal_peak_sets[[paste0("fetal_", zct)]]), " peaks")
}

# Fetal union (all Ziffra peaks)
fetal_peak_sets[["fetal_union"]] <- peak_list %>%
  makeGRangesFromDataFrame(keep.extra.columns = TRUE)
message("  fetal_union: ", length(fetal_peak_sets[["fetal_union"]]), " peaks")

# 5. Run enrichment tests ---------------------------------------------------------------
all_peak_sets <- c(fetal_peak_sets, adult_peak_sets)

message("\n--- Running enrichment tests (", length(all_peak_sets), " peak sets x 2 SNP sets) ---\n")

results <- tibble()

for (peak_name in names(all_peak_sets)) {
  peaks_gr <- all_peak_sets[[peak_name]]
  stage <- ifelse(str_starts(peak_name, "fetal_"), "fetal", "adult")

  for (snp_set_name in c("fetal_specific", "shared")) {
    snp_gr <- if (snp_set_name == "fetal_specific") specific_gr else shared_gr

    message("\n=== ", peak_name, " vs ", snp_set_name, " SNPs (",
            length(snp_gr), " SNPs, ", length(peaks_gr), " peaks) ===")

    perm_result <- perform_permutation_test(snp_gr, peaks_gr, n_perm = 1000)

    results <- bind_rows(
      results,
      tibble(
        peak_set       = peak_name,
        stage          = stage,
        snp_set        = snp_set_name,
        n_snps         = length(snp_gr),
        n_peaks        = length(peaks_gr),
        p_value        = perm_result$p_value,
        observed_overlap     = perm_result$observed_overlap,
        mean_permuted_overlap = perm_result$mean_permuted_overlap,
        fold_enrichment      = perm_result$fold_enrichment
      )
    )
  }
}

# 6. Summary ----------------------------------------------------------------------------
message("\n--- Results summary ---\n")
print(results %>% select(-n_snps, -n_peaks), n = Inf)

message("\n\nFetal-specific SNPs in fetal vs adult peaks (union sets):")
results %>%
  filter(snp_set == "fetal_specific", peak_set %in% c("fetal_union", "adult_cortical_union")) %>%
  select(peak_set, fold_enrichment, p_value) %>%
  print()

# Write single tracked output ----------------------------------------------------------
write_rds(results, output_file)
message("\nOCR enrichment results written to: ", output_file)
