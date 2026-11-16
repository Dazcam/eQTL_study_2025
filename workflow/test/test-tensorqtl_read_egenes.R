#--------------------------------------------------------------------------------------
#
#    test-tensorqtl_read_egenes.R
#      - Used in tensorqtl_report.Rmd (function in tensorqtl_functions.R)
#
#--------------------------------------------------------------------------------------

library(testthat)
suppressMessages({
  library(dplyr)
  library(tibble)
})

# Toy tensorQTL permutation output (cis_qtl.txt.gz), columns as tensorQTL writes them
toy_perm <- function(path, phenotype_id, start_distance, qval) {
  n <- length(phenotype_id)
  tbl <- tibble(phenotype_id = phenotype_id, num_var = 5000, beta_shape1 = 1, beta_shape2 = 900,
                true_df = 100, pval_true_df = 1e-5, variant_id = paste0("chr1:", seq_len(n), ":A:G"),
                start_distance = start_distance, end_distance = start_distance, ma_samples = 50,
                ma_count = 60, af = 0.3, pval_nominal = 1e-6, slope = 0.3, slope_se = 0.05,
                pval_perm = 1e-4, pval_beta = 1e-4, qval = qval, pval_nominal_threshold = 1e-4)
  con <- gzfile(path, "w")
  write.table(tbl, con, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
  close(con)
}

quiet <- function(expr) {
  out <- NULL
  invisible(capture.output(out <- suppressMessages(expr)))
  out
}

in_dir <- paste0(tempfile(), "/")
d <- tqtl_dir(in_dir, "GABA-1", "quantile", 4, 40); dir.create(d, recursive = TRUE)
toy_perm(paste0(d, "GABA-1_quantile_perm.cis_qtl.txt.gz"),
         c("g1", "g1", "g2", "g3", "g4"), c(10, 20, NA, 5, 7), c(0.01, 0.02, 0.03, NA, 0.9))

test_that("distinct significant eGenes, including those without a TSS distance", {
  res <- read_egenes("GABA-1", in_dir, "quantile", 4, 40, fdr = 0.05)
  expect_equal(res$phenotype_id, c("g1", "g2"))
  expect_equal(unique(res$cell_type), "GABA-1")
})

test_that("a missing file gives an empty result", {
  res <- read_egenes("GABA-2", in_dir, "quantile", 4, 40, fdr = 0.05)
  expect_equal(nrow(res), 0)
})
