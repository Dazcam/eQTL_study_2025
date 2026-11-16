#--------------------------------------------------------------------------------------
#
#    test-tensorqtl_parse_perm_log.R
#      - Used in tensorqtl_report.Rmd (function in tensorqtl_functions.R)
#
#--------------------------------------------------------------------------------------

library(testthat)
suppressMessages({
  library(dplyr)
  library(tibble)
})

# Toy tensorQTL cis (permutation) log, in the layout tensorQTL prints
toy_log <- function(path, n_qtl = 1234) {
  writeLines(c(
    "[Nov 01 10:00:00] Running TensorQTL: cis-QTL mapping",
    "cis-QTL mapping: empirical p-values for phenotypes",
    "  * 120 samples",
    "  * 15000 phenotypes",
    "  * 36 covariates",
    "  * 6123456 variants",
    "  * cis-window: ±500,000",
    "  Time elapsed: 42.17 min",
    "Computing q-values",
    "  * Number of phenotypes tested: 14876",
    "  * Correlation between Beta-approximated and empirical p-values: 0.9987",
    "  * Proportion of significant phenotypes (1-pi0): 0.35",
    paste0("  * QTL phenotypes @ FDR 0.05: ", n_qtl),
    "  * min p-value threshold @ FDR 0.05: 0.000123"), path)
}

test_that("extracts every summary number from the log", {
  f <- tempfile(fileext = ".log"); toy_log(f)
  res <- parse_perm_log(f)
  expect_equal(nrow(res), 1)
  expect_equal(res$genes_tested, 14876L)   # 'phenotypes tested', not the 15000 loaded
  expect_equal(res$samples, 120L)
  expect_equal(res$covariates, 36L)
  expect_equal(res$variants, 6123456L)
  expect_equal(res$qtl_fdr_0.05, 1234L)
  expect_equal(res$cis_window_kb, 500L)
  expect_equal(res$beta_cor, 0.9987)
  expect_equal(res$min_p, "1.23e-04")
  expect_equal(res$prop_sig, 35L)
  expect_equal(res$run_time_min, 42L)
})

test_that("takes the last 'covariates' line when there are several", {
  f <- tempfile(fileext = ".log"); toy_log(f)
  writeLines(c("  * 2 covariates (sex, PCW)", readLines(f)), f)
  expect_equal(parse_perm_log(f)$covariates, 36L)
})
