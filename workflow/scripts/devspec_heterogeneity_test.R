#--------------------------------------------------------------------------------------
#
#    Developmental specificity: pairwise Cochran's Q heterogeneity test
#
#--------------------------------------------------------------------------------------

#  Tests, for each of the 9 pairwise fetal (L1) - Fujita adult comparisons, whether
#  the fetal and (allele-harmonized) adult beta for every tested-in-both SNP-gene
#  pair differ by more than combined sampling error predicts:
#
#    Q = (beta_my - beta_fugita_harmonized)^2 / (se_my^2 + se_fugita^2) ~ chisq(df=1)
#
#  Pairs with allele_status in {palindromic_ambiguous, mismatched, missing_allele_info}
#  (beta_fugita_harmonized = NA, ~2.55% of the background per
#  devspec_extract_pairwise_universe.R's allele_summary) are dropped here, at the
#  start of this script, not upstream -- so the allele-harmonization QC step and
#  the statistical test remain cleanly separated.
#
#  FDR-corrected within each comparison (BH) -- tested-pair N varies ~5-fold across
#  comparisons, so pooling would let the largest comparisons set the rejection
#  threshold applied to the smallest.
#
#  Classification, from fetal_sig x adult_sig x heterogeneous(FDR<thresh):
#    fetal_sig  adult_sig  heterogeneous  ->  category
#    T          T          F              ->  shared_homogeneous
#    T          T          T              ->  shared_heterogeneous
#    T          F          F              ->  fetal_sig_not_heterogeneous
#    T          F          T              ->  fetal_specific_heterogeneous   (headline)
#    F          T          F              ->  adult_sig_not_heterogeneous
#    F          T          T              ->  adult_specific_heterogeneous
#    F          F          -              ->  background
#
#  Per comparison, tests whether fetal-significant pairs are enriched for
#  heterogeneity relative to the full tested background (Fisher's exact, 2x2).
#
#  Step 2 (revised) of the developmental-specificity pipeline. Standalone -- does
#  not touch classify_specificity.R's existing output.

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

message("\n\nRunning pairwise Cochran's Q heterogeneity test (fetal vs adult) ...")

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
})

# Input and output paths ---------------------------------------------------------------
background_file <- snakemake@input[["background"]]
output_file      <- snakemake@output[[1]]
fdr_thresh        <- snakemake@params[["fdr_thresh"]]

if (is.null(fdr_thresh)) fdr_thresh <- 0.05

# Check variable assignment
message("\nVariables")
cat("============================")
tibble(
  variable = c("background_file", "fdr_thresh", "output_file"),
  value    = c(background_file, fdr_thresh, output_file)
) |>
  knitr::kable(format = "simple", align = "l") |>
  print()
message("\n============================\n")

# ========================================================================================
# 1. Load background
# ========================================================================================

message("\nLoading pairwise tested-pair background ...")
bg <- read_rds(background_file)
background_list <- bg$background
comparison_map   <- bg$comparison_map

message("  Comparisons: ", length(background_list))

# ========================================================================================
# 2. Run Q-test per comparison
# ========================================================================================

message("\n--- Running Cochran's Q per comparison ---\n")

pair_results_list       <- list()
comparison_summary_list <- list()
category_summary_list   <- list()

for (comp in names(background_list)) {

  message("\n=== ", comp, " ===")

  dt <- as.data.table(background_list[[comp]])
  n_tested_total <- nrow(dt)

  # Drop pairs whose adult beta could not be safely harmonized (palindromic
  # ambiguous, mismatched, missing allele info). Logged explicitly here per
  # comparison, mirroring the allele_summary breakdown from step 1.
  n_dropped <- sum(is.na(dt$beta_fugita_harmonized))
  dt <- dt[!is.na(beta_fugita_harmonized)]
  n_usable <- nrow(dt)

  message("  Tested pairs (total): ", n_tested_total)
  message("  Dropped (unresolved allele harmonization): ", n_dropped,
          " (", round(100 * n_dropped / n_tested_total, 2), "%)")
  message("  Usable pairs (entering Q-test): ", n_usable)

  # Cochran's Q (= Z^2), difference, and combined SE -- vectorized
  dt[, delta_beta := beta_my - beta_fugita_harmonized]
  dt[, se_delta   := sqrt(se_my^2 + se_fugita^2)]
  dt[, Q := (delta_beta / se_delta)^2]
  dt[, p := pchisq(Q, df = 1, lower.tail = FALSE)]

  # FDR-correct within this comparison only
  dt[, q := p.adjust(p, method = "BH")]
  dt[, heterogeneous := q < fdr_thresh]

  message("  Heterogeneous pairs (FDR<", fdr_thresh, "): ", sum(dt$heterogeneous),
          " (", round(100 * mean(dt$heterogeneous), 2), "%)")

  # Classification
  dt[, category := fcase(
    fetal_sig  &  adult_sig  & !heterogeneous, "shared_homogeneous",
    fetal_sig  &  adult_sig  &  heterogeneous, "shared_heterogeneous",
    fetal_sig  & !adult_sig  & !heterogeneous, "fetal_sig_not_heterogeneous",
    fetal_sig  & !adult_sig  &  heterogeneous, "fetal_specific_heterogeneous",
    !fetal_sig &  adult_sig  & !heterogeneous, "adult_sig_not_heterogeneous",
    !fetal_sig &  adult_sig  &  heterogeneous, "adult_specific_heterogeneous",
    !fetal_sig & !adult_sig,                   "background"
  )]

  message("  Category counts:")
  cat_tbl <- dt[, .N, by = category][order(-N)]
  print(cat_tbl)

  category_summary_list[[comp]] <- as_tibble(cat_tbl) %>%
    mutate(comparison = comp, prop = round(N / n_usable, 4))

  # Enrichment: are fetal-sig pairs enriched for heterogeneity vs full background?
  fisher_tbl <- table(
    fetal_sig     = dt$fetal_sig,
    heterogeneous = dt$heterogeneous
  )
  message("  2x2 table (fetal_sig x heterogeneous):")
  print(fisher_tbl)

  fisher_res <- fisher.test(fisher_tbl)
  message("  Fisher's exact OR = ", round(fisher_res$estimate, 3),
          ", p = ", signif(fisher_res$p.value, 3))

  comparison_summary_list[[comp]] <- tibble(
    comparison       = comp,
    n_tested_total    = n_tested_total,
    n_dropped          = n_dropped,
    n_usable            = n_usable,
    n_fetal_sig          = sum(dt$fetal_sig),
    n_adult_sig          = sum(dt$adult_sig),
    n_heterogeneous      = sum(dt$heterogeneous),
    n_fetal_specific     = sum(dt$category == "fetal_specific_heterogeneous"),
    n_adult_specific     = sum(dt$category == "adult_specific_heterogeneous"),
    fisher_OR             = unname(fisher_res$estimate),
    fisher_p               = fisher_res$p.value
  )

  pair_results_list[[comp]] <- as_tibble(dt)

  rm(dt); gc()
}

# ========================================================================================
# 3. Combine and summarise
# ========================================================================================

message("\n--- Combined summary across all comparisons ---\n")

pair_results        <- bind_rows(pair_results_list)
comparison_summary  <- bind_rows(comparison_summary_list)
category_summary    <- bind_rows(category_summary_list)

message("\nComparison-level summary:")
comparison_summary %>%
  knitr::kable(format = "simple", digits = 4) %>%
  print()

message("\nCategory counts by comparison:")
category_summary %>%
  select(comparison, category, N, prop) %>%
  knitr::kable(format = "simple") %>%
  print()

message("\nTotal usable rows across all comparisons: ", nrow(pair_results))

# Pooled fetal-specific SNP list -- for ocr_enrich.R to consume in place of the
# old any_table$final_category == "testable_specific" filter
fetal_specific_snps <- pair_results %>%
  filter(category == "fetal_specific_heterogeneous") %>%
  pull(variant_id) %>%
  unique()

message("\nPooled fetal-specific (heterogeneous) unique SNPs: ", length(fetal_specific_snps))

# ========================================================================================
# 4. Write output
# ========================================================================================

write_rds(
  list(
    pair_results         = pair_results,
    comparison_summary    = comparison_summary,
    category_summary      = category_summary,
    fetal_specific_snps    = fetal_specific_snps
  ),
  output_file
)
message("\nHeterogeneity test results written to: ", output_file)
