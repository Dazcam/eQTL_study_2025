#--------------------------------------------------------------------------------------
#
#    Developmental specificity: OCR enrichment of fetal-specific/shared/
#    heterogeneous eQTL (adult cortical union only)
#
#--------------------------------------------------------------------------------------

#  Tests whether three fetal SNP sets are enriched in adult open chromatin
#  (Li et al. 2023, cortical union only): testable_specific and testable_shared
#  (from B's "any" tier, unchanged), and the new specific_heterogeneous set (from
#  C -- testable_specific pairs additionally confirmed as formally heterogeneous
#  in >=1 tested adult comparison). All three LD-pruned against the same
#  independent_snps universe (D) before testing.
#
#  NOTE: specific_heterogeneous is a SUBSET of testable_specific, not an
#  independent comparison -- interpret side-by-side results with that nesting
#  in mind, not as three unrelated tests.
#
#  Step 5 of the developmental-specificity pipeline.

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

message("\n\nTesting fetal-specific/shared/heterogeneous eQTL enrichment in adult cortical union OCRs ...")

suppressPackageStartupMessages({
  library(tidyverse)
  library(GenomicRanges)
  library(nullranges)
  library(GenomeInfoDb)
})

# Input and output paths ---------------------------------------------------------------
classification_file    <- snakemake@input[["classification"]]
heterogeneity_file     <- snakemake@input[["heterogeneity_test"]]
independent_snps_file  <- snakemake@input[["independent_snps"]]
snp_lookup_files        <- snakemake@input[["snp_lookups"]]
li_s6_file             <- snakemake@input[["li_s6"]]
li_s7_file             <- snakemake@input[["li_s7"]]
output_file            <- snakemake@output[[1]]

# Check variable assignment
message("\nVariables")
cat("============================")
tibble(
  variable = c("classification_file", "heterogeneity_file", "independent_snps_file",
               "n_snp_lookups", "li_s6_file", "li_s7_file", "output_file"),
  value    = c(classification_file, heterogeneity_file, independent_snps_file,
               length(snp_lookup_files), li_s6_file, li_s7_file, output_file)
) |>
  knitr::kable(format = "simple", align = "l") |>
  print()
message("\n============================\n")

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

# Functions ----------------------------------------------------------------------------

perform_permutation_test <- function(sig_snps_gr, peaks_gr, n_perm = 1000) {
  message("Starting permutation test for ATAC-seq peak enrichment...")

  if (!exists("bootRanges", envir = asNamespace("nullranges"))) {
    stop("bootRanges function not found. Ensure nullranges is installed and loaded.")
  }

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

# 1. Load fetal-specific + shared SNPs from classification, and the new
#    heterogeneous set from C. All restricted to the LD-independent universe. ----------
message("\nLoading classification output (B) ...")
classification <- read_rds(classification_file)
any_table <- classification$any_table

message("Loading heterogeneity test output (C) ...")
het_test <- read_rds(heterogeneity_file)
heterogeneous_pairs <- het_test$heterogeneous_pairs %>% filter(heterogeneous_all)

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

# 2. Build SNP coordinate lookup from existing per-cell-type files ----------------------
message("\nBuilding SNP coordinate lookup from ", length(snp_lookup_files), " existing files ...")

snp_coords <- map_dfr(snp_lookup_files, read_rds) %>%
  distinct(refsnp_id, .keep_all = TRUE)
message("  Total unique SNPs with coordinates: ", nrow(snp_coords))

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

specific_gr             <- make_snp_gr(fetal_specific_snps, snp_coords)
shared_gr               <- make_snp_gr(shared_snps, snp_coords)
specific_heterogeneous_gr <- make_snp_gr(specific_heterogeneous_snps, snp_coords)

message("  Fetal-specific SNPs with coordinates: ", length(specific_gr))
message("  Shared SNPs with coordinates: ", length(shared_gr))
message("  Specific+heterogeneous SNPs with coordinates: ", length(specific_heterogeneous_gr))

# 3. Load Li et al. 2023 adult peaks, cortical union only -------------------------------
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

# 4. Run the three enrichment tests -------------------------------------------------------
message("\n--- Running enrichment tests (3 sets x adult_cortical_union) ---\n")

snp_sets <- list(
  fetal_specific = specific_gr,
  shared         = shared_gr,
  specific_heterogeneous = specific_heterogeneous_gr
)

results <- tibble()

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

# 5. Summary ------------------------------------------------------------------------------
message("\n--- Results summary ---\n")
print(results, n = Inf)

# Write tracked output ----------------------------------------------------------------
write_rds(results, output_file)
message("\nOCR enrichment results written to: ", output_file)

# Untracked TSV for readability outside R -----------------------------------------------
output_tsv <- sub("\\.rds$", ".tsv", output_file)
write_tsv(results, output_tsv)
message("OCR enrichment results (TSV) written to: ", output_tsv)
