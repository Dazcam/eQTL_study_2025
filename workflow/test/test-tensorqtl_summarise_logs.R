#--------------------------------------------------------------------------------------
#
#    test-tensorqtl_summarise_logs.R
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

make_logs <- function(in_dir, runs) {
  for (i in seq_len(nrow(runs))) {
    d <- tqtl_dir(in_dir, runs$ct[i], runs$nm[i], 4, runs$pc[i])
    dir.create(d, recursive = TRUE, showWarnings = FALSE)
    toy_log(paste0(d, runs$ct[i], "_", runs$nm[i], "_perm.tensorQTL.cis.log"), runs$n[i])
  }
}

test_that("one row per log, cell type fastest, then expPC, then norm method", {
  in_dir <- paste0(tempfile(), "/")
  runs <- tidyr::expand_grid(nm = c("bryois", "quantile"), pc = c(10, 20), ct = c("NPC", "MG")) |>
    mutate(n = row_number() * 10L)
  make_logs(in_dir, runs)
  res <- suppressMessages(summarise_logs(c("NPC", "MG"), in_dir, c(10, 20), 4, c("bryois", "quantile")))
  expect_equal(nrow(res), 8)
  expect_equal(names(res)[1:5], c("cell_type", "genPC", "expPC", "norm_methods", "genes_tested"))
  expect_equal(res$cell_type, rep(c("NPC", "MG"), 4))
  expect_equal(res$expPC, rep(c(10, 10, 20, 20), 2))
  expect_equal(res$norm_methods, rep(c("bryois", "quantile"), each = 4))
  expect_equal(res$qtl_fdr_0.05, runs$n)
})

test_that("missing logs are skipped", {
  in_dir <- paste0(tempfile(), "/")
  make_logs(in_dir, tibble(ct = "NPC", nm = "quantile", pc = 30, n = 99L))
  res <- suppressMessages(summarise_logs(c("NPC", "OPC"), in_dir, 30, 4, "quantile"))
  expect_equal(res$cell_type, "NPC")
  expect_equal(res$qtl_fdr_0.05, 99L)
})
