#--------------------------------------------------------------------------------------
#
#    test-repl_sc_align_jang_beta.R
#      - Used in repl_sc_pi1_enrichment.R and repl_sc_beta_correlation.R
#        (function in repl_sc_functions.R)
#
#--------------------------------------------------------------------------------------

library(testthat)

beta   <- c(0.5, 0.5, 0.5, 0.5)
allele <- c("T", "C", "G", NA)
ref    <- c("C", "C", "C", "C")
alt    <- c("T", "T", "T", "T")

test_that("without ref: beta kept where Allele is alt, negated otherwise", {
  expect_equal(align_jang_beta(beta, allele, alt), c(0.5, -0.5, -0.5, NA))
})

test_that("with ref: negated only where Allele is ref, NA if neither", {
  expect_equal(align_jang_beta(beta, allele, alt, ref), c(0.5, -0.5, NA, NA))
})
