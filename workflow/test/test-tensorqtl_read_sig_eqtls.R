#--------------------------------------------------------------------------------------
#
#    test-tensorqtl_read_sig_eqtls.R
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
for (ct in c("NPC", "MG")) {
  pc <- c(NPC = 30, MG = 20)[[ct]]
  d <- tqtl_dir(in_dir, ct, "quantile", 4, pc); dir.create(d, recursive = TRUE)
  if (ct == "NPC") {
    toy_perm(paste0(d, "NPC_quantile_perm.cis_qtl.txt.gz"),
             c("g1", "g2", "g3", "g4"), c(-2500, 1000, NA, 400), c(0.01, 0.2, 0.001, 0.0499))
  } else {
    toy_perm(paste0(d, "MG_quantile_perm.cis_qtl.txt.gz"), "g9", 100, 0.5)  # nothing significant
  }
}

test_that("keeps significant eQTL with a TSS distance, from the mapped expPC run", {
  res <- quiet(read_sig_eqtls(c("NPC", "MG", "OPC"), in_dir, c(NPC = 30, MG = 20, OPC = 30),
                              "quantile", 4, fdr = 0.05))
  expect_equal(res$phenotype_id, c("g1", "g4"))   # g2 not significant, g3 no distance
  expect_equal(res$distance_kb, c(-2.5, 0.4))
  expect_equal(unique(res$cell_type), "NPC")       # MG has none; OPC has no file
})

test_that("the FDR threshold is applied as qval < fdr", {
  res <- quiet(read_sig_eqtls("NPC", in_dir, c(NPC = 30), "quantile", 4, fdr = 0.01))
  expect_equal(nrow(res), 0)
})
