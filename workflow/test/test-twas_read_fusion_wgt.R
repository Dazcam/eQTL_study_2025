#--------------------------------------------------------------------------------------
#
#    test-read_fusion_wgt.R
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

## Toy weight files  ------------------------------------------------------------------
# Same objects and layout as FUSION.compute_weights.R saves
toy_dir <- tempfile("wgt_")
dir.create(toy_dir)

make_wgt <- function(gene, wgt.matrix, cv.performance, hsq = c(0.2, 0.05), hsq.pv = 1e-4,
                     N.tot = 134) {
  snps <- data.frame(V1 = 1, V2 = rownames(wgt.matrix), V3 = 0,
                     V4 = seq_len(nrow(wgt.matrix)), V5 = "A", V6 = "G")
  file <- file.path(toy_dir, paste0(gene, ".wgt.RDat"))
  save(wgt.matrix, snps, cv.performance, hsq, hsq.pv, N.tot, file = file)
  file
}

models <- c("top1", "lasso", "enet")

# G1: 4 cis SNPs, 2 non-zero lasso weights, enet has the smallest CV p-value
wgt_g1 <- matrix(c(0, 0, 0, 5,  0.1, 0, 0.2, 0,  0.1, 0.1, 0.1, 0.1), ncol = 3,
                 dimnames = list(paste0("rs", 1:4), models))
cv_g1 <- matrix(c(0.05, 0.02,  0.10, 1e-3,  0.12, 1e-4), nrow = 2,
                dimnames = list(c("rsq", "pval"), models))
f1 <- make_wgt("G1", wgt_g1, cv_g1)

# G2: lasso failed (all NA), only top1 and lasso run
wgt_g2 <- matrix(c(3, 0, NA, NA), ncol = 2, dimnames = list(c("rs1", "rs2"), c("top1", "lasso")))
cv_g2 <- matrix(c(0.03, 0.04, NA, NA), nrow = 2,
                dimnames = list(c("rsq", "pval"), c("top1", "lasso")))
f2 <- make_wgt("G2", wgt_g2, cv_g2, hsq = c(0.1, 0.04), hsq.pv = 0.005, N.tot = 130)

## Tests  -----------------------------------------------------------------------------
test_that("gene, sample size and heritability are read correctly", {
  out <- read_fusion_wgt(f1, models)
  expect_equal(out$gene, "G1")
  expect_equal(out$n_tot, 134)
  expect_equal(out$hsq, 0.2)
  expect_equal(out$hsq_se, 0.05)
  expect_equal(out$hsq_p, 1e-4)
})

test_that("SNP counts and CV performance are read correctly", {
  out <- read_fusion_wgt(f1, models)
  expect_equal(out$n_snps_cis, 4)
  expect_equal(out$n_snps_lasso, 2)
  expect_equal(out$rsq_lasso, 0.10)
  expect_equal(out$pval_lasso, 1e-3)
  expect_equal(out$rsq_enet, 0.12)
})

test_that("best model is the one with the smallest CV p-value", {
  expect_equal(read_fusion_wgt(f1, models)$best_model, "enet")
  expect_equal(read_fusion_wgt(f2, models)$best_model, "top1")
})

test_that("failed lasso gives NA, and a model that was not run gives NA columns", {
  out <- read_fusion_wgt(f2, models)
  expect_true(is.na(out$n_snps_lasso))
  expect_true(is.na(out$rsq_lasso))
  expect_true("rsq_enet" %in% colnames(out))
  expect_true(is.na(out$rsq_enet))
})

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
