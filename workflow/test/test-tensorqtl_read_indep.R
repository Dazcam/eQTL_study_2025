#--------------------------------------------------------------------------------------
#
#    test-tensorqtl_read_indep.R
#      - Used in tensorqtl_report.Rmd (function in tensorqtl_functions.R)
#
#--------------------------------------------------------------------------------------

library(testthat)
suppressMessages({
  library(dplyr)
  library(tibble)
})

test_that("reads each cell type's independent eQTL from its mapped expPC run", {
  ind_dir <- paste0(tempfile(), "/")
  for (ct in c("NPC", "MG")) {
    d <- tqtl_dir(ind_dir, ct, "quantile", 4, 30); dir.create(d, recursive = TRUE)
    tbl <- tibble(phenotype_id = c("g1", "g1", "g2"), variant_id = c("v1", "v2", "v3"),
                  start_distance = c(100, -200, 300), af = 0.2, pval_perm = 1e-4, rank = c(1, 2, 1))
    con <- gzfile(paste0(d, ct, "_quantile_independent.cis_independent_qtl.txt.gz"), "w")
    write.table(tbl, con, sep = "\t", quote = FALSE, row.names = FALSE); close(con)
  }
  res <- read_indep(c("NPC", "MG"), ind_dir, c(NPC = 30, MG = 30), "quantile", 4)
  expect_equal(nrow(res), 6)
  expect_equal(res$cell_type, rep(c("NPC", "MG"), each = 3))
  expect_equal(res$rank, c(1, 2, 1, 1, 2, 1))
})
