#--------------------------------------------------------------------------------------
#
#    test-repl_sc_compute_pi1.R
#      - Used in repl_sc_pi1_enrichment.R (function in repl_sc_functions.R)
#
#--------------------------------------------------------------------------------------

library(testthat)
suppressMessages({
  library(dplyr)
  library(tibble)
})

pi1_quiet <- function(...) {
  out <- NULL
  invisible(capture.output(out <- suppressMessages(compute_pi1(...))))
  out
}

# Query: 2000 significant pairs (large enough for a stable pi0 estimate)
query <- tibble(variant_id   = paste0("rs", 1:2000),
                phenotype_id = paste0("ENSG", 1:2000),
                slope_my     = rep(c(0.5, -0.5), 1000))

test_that("all-null reference p-values give pi1 close to 0", {
  set.seed(1)
  ref <- tibble(variant_id = query$variant_id, phenotype_id = query$phenotype_id,
                pval = runif(2000), slope_ref = 0.1)
  res <- pi1_quiet(query, ref, "toy")
  expect_equal(res$overlap_count, 2000)
  expect_lt(res$pi1, 0.15)
})

# 1500 strong signals (p < 1e-6) and 500 nulls spread over 0.2-1. pi0 is
# estimated from p > lambda, so these nulls give pi0 ~ 500 / 0.8 / 2000 = 0.31
set.seed(2)
mixed_p <- c(runif(1500, 0, 1e-6), runif(500, 0.2, 1))

test_that("a 75% signal mixture gives pi1 near 0.69 and 75% replicating", {
  ref <- tibble(variant_id = query$variant_id, phenotype_id = query$phenotype_id,
                pval = mixed_p, slope_ref = query$slope_my)
  res <- pi1_quiet(query, ref, "toy")
  expect_gt(res$pi1, 0.55)
  expect_lt(res$pi1, 0.8)
  expect_equal(res$prop_replicating, 0.75)
  expect_equal(res$prop_same_direction, 1)
})

test_that("prop_same_direction counts sign agreement among replicating pairs only", {
  # 300 of the 1500 replicating pairs have the opposite sign; null-pair signs don't count
  slope_ref <- query$slope_my
  slope_ref[1:300]     <- -slope_ref[1:300]
  slope_ref[1501:2000] <- -slope_ref[1501:2000]
  ref <- tibble(variant_id = query$variant_id, phenotype_id = query$phenotype_id,
                pval = mixed_p, slope_ref = slope_ref)
  res <- pi1_quiet(query, ref, "toy")
  expect_equal(res$prop_same_direction, 0.8)
})

test_that("pairs missing from the reference or with invalid p-values are excluded", {
  ref <- tibble(variant_id = query$variant_id[1:120], phenotype_id = query$phenotype_id[1:120],
                pval = c(mixed_p[1:100], rep(NA, 10), rep(1.5, 5), rep(Inf, 5)),
                slope_ref = 0.1)
  res <- pi1_quiet(query, ref, "toy")
  expect_equal(res$overlap_count, 100)
  expect_length(res$pvals, 100)
})

test_that("fewer than min_overlap pairs returns NA statistics", {
  ref <- tibble(variant_id = query$variant_id[1:40], phenotype_id = query$phenotype_id[1:40],
                pval = 1e-10, slope_ref = 0.1)
  res <- pi1_quiet(query, ref, "toy")
  expect_equal(res$overlap_count, 40)
  expect_true(is.na(res$pi1))
  expect_true(is.na(res$prop_replicating))
  expect_true(is.na(res$prop_same_direction))
})
