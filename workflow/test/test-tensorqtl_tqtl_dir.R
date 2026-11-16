#--------------------------------------------------------------------------------------
#
#    test-tensorqtl_tqtl_dir.R
#      - Used in tensorqtl_report.Rmd (function in tensorqtl_functions.R)
#
#--------------------------------------------------------------------------------------

library(testthat)
suppressMessages({
  library(dplyr)
  library(tibble)
})

test_that("builds the run directory used by the tensorQTL rules", {
  expect_equal(tqtl_dir("../results/05TENSORQTL/tensorqtl_perm/", "Glu-UL-1", "quantile", 4, 50),
               "../results/05TENSORQTL/tensorqtl_perm/Glu-UL-1_quantile_genPC_4_expPC_50/")
})

test_that("is vectorised over cell types", {
  expect_equal(tqtl_dir("d/", c("NPC", "MG"), "xue", 4, 30),
               c("d/NPC_xue_genPC_4_expPC_30/", "d/MG_xue_genPC_4_expPC_30/"))
})
