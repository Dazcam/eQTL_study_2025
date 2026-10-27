#--------------------------------------------------------------------------------------
#
#    test-classify_weights_log.R
#
#--------------------------------------------------------------------------------------
#
# - Used in twas_weights_summary.Rmd
#
#--------------------------------------------------------------------------------------

library(testthat)
library(dplyr)
library(tibble)
library(stringr)

## Toy log  ---------------------------------------------------------------------------
# Same structure as a compute_weights log: rule lines plus FUSION messages, in
# FUSION.compute_weights.R's exact wording. One gene for each outcome.
toy_log <- c(
  "Starting weight computation for Glu-UL",
  # G1: weights computed
  "Processing G1",
  "Loaded 46 covariates",
  "### Estimating heritability",
  "Heritability (se): 0.21 0.08 LRT P-value: 0.001",
  "0 variants removed due to missing genotype data",  # PLINK-style line with spaces
  "Computing full-sample weights",
  "../results/11TWAS/fusion_input/Glu-UL/G1 removed",
  # G2: heritability not significant
  "Processing G2",
  "Loaded 46 covariates",
  "### Estimating heritability",
  "../results/11TWAS/weights/Glu-UL/G2.tmp : heritability 0.02; LRT P-value 0.31 : skipping gene",
  "../results/11TWAS/fusion_input/Glu-UL/G2 removed",
  # G3: GCTA failed
  "Processing G3",
  "Loaded 46 covariates",
  "../results/11TWAS/weights/Glu-UL/G3.tmp does not exist, likely GCTA could not converge, skipping gene",
  "../results/11TWAS/fusion_input/Glu-UL/G3 removed",
  # G4: no cis SNPs, so FUSION never ran
  "Processing G4",
  "../results/11TWAS/fusion_input/Glu-UL/G4 removed",
  # G5: no cis SNPs and the clean-up failed
  "Processing G5",
  "WARNING: Failed to remove ../results/11TWAS/fusion_input/Glu-UL/G5",
  "All genes completed"
)

## Tests  -----------------------------------------------------------------------------
test_that("one row per gene, in log order", {
  out <- classify_weights_log(toy_log)
  expect_equal(out$gene, c("G1", "G2", "G3", "G4", "G5"))
})

test_that("each outcome is classified correctly", {
  out <- classify_weights_log(toy_log)
  expect_equal(out$log_status,
               c("fusion_ran", "hsq_skip", "gcta_failed", "no_fusion_output", "no_fusion_output"))
})

test_that("PLINK lines ending in 'removed' are not mistaken for rule lines", {
  only_plink <- c("Processing G9", "0 variants removed due to missing genotype data")
  expect_equal(classify_weights_log(only_plink)$log_status, "fusion_ran")
})

test_that("a log with no genes returns an empty tibble", {
  out <- classify_weights_log(c("Starting weight computation for Glu-UL"))
  expect_equal(nrow(out), 0)
  expect_equal(colnames(out), c("gene", "log_status"))
})

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
