#--------------------------------------------------------------------------------------
#
#    test-repl_sc_allele_harmonisation.R
#      - Used in repl_sc_beta_correlation.R (get_allele_status() and
#        harmonise_beta() in repl_sc_functions.R)
#
#--------------------------------------------------------------------------------------

library(testthat)

#              aligned  flipped  palindromic  mismatched  missing  pal+flip
ref_my  <- c("A",      "A",     "A",         "A",        NA,      "A")
alt_my  <- c("G",      "G",     "T",         "G",        "G",     "T")
ref_ref <- c("A",      "G",     "A",         "C",        "A",     "T")
alt_ref <- c("G",      "A",     "T",         "T",        "G",     "A")

test_that("each allele combination gets the expected status", {
  expect_equal(get_allele_status(ref_my, alt_my, ref_ref, alt_ref),
               c("aligned", "flipped", "palindromic_unresolved", "mismatched",
                 "missing_allele_info", "palindromic_unresolved"))
})

test_that("a missing reference beta counts as missing only when beta_ref is given", {
  beta <- c(NA, 0.2, 0.2, 0.2, 0.2, 0.2)
  expect_equal(get_allele_status(ref_my, alt_my, ref_ref, alt_ref)[1], "aligned")
  expect_equal(get_allele_status(ref_my, alt_my, ref_ref, alt_ref, beta)[1],
               "missing_allele_info")
})

test_that("harmonise_beta keeps aligned, negates flipped and drops the rest", {
  status <- c("aligned", "flipped", "palindromic_unresolved", "mismatched", "missing_allele_info")
  expect_equal(harmonise_beta(rep(0.3, 5), status), c(0.3, -0.3, NA, NA, NA))
})
