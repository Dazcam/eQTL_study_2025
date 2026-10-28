# test-susie_extractResults.R
# - Used in susie_run.R

library(testthat)
library(dplyr)
library(tibble)

# Toy SuSiE fit: 4 variants, L = 3 effects. Sets are in SuSiE's order (highest purity
# first): L1 = {v1, v2} high purity; L3 = {v2, v3} overlaps L1; L2 = {v4} low purity
make_fit <- function() {
  ids <- c("chr1_100_A_G", "chr1_200_A_G", "chr1_300_A_G", "chr1_400_A_G")
  fit <- list(
    alpha = rbind(c(0.6, 0.4, 0, 0), c(0, 0, 0, 1), c(0, 0.5, 0.5, 0)),
    mu    = rbind(c(1, 2, 0, 0), c(0, 0, 0, 3), c(0, 1, 1, 0)),
    mu2   = rbind(c(2, 5, 0, 0), c(0, 0, 0, 10), c(0, 2, 2, 0)),
    V = c(1, 1, 1),
    X_column_scale_factors = c(1, 1, 1, 1),
    pip = c(0.6, 0.7, 0.5, 1),
    z = c(4, 3, 2, 5),
    lbf = c(10, 3, 5),
    variant_id = ids,
    sets = list(
      cs = list(L1 = c(1, 2), L3 = c(2, 3), L2 = 4),
      purity = data.frame(min.abs.corr = c(0.9, 0.8, 0.3),
                          mean.abs.corr = c(0.95, 0.85, 0.3),
                          median.abs.corr = c(0.95, 0.85, 0.3),
                          row.names = c("L1", "L3", "L2")),
      cs_index = c(1, 3, 2)))
  # In a real fit the variant IDs are column names (from the genotype matrix)
  colnames(fit$alpha) <- colnames(fit$mu) <- colnames(fit$mu2) <- ids
  fit
}

test_that("overlapping credible sets are dropped and purity summaries are correct", {
  res <- extractResults(make_fit())

  expect_equal(res$cs_df$cs_id, c("L1", "L2"))
  expect_equal(res$cs_df$cs_size, c(2, 1))
  expect_equal(res$cs_df$cs_log10bf, log10(exp(c(10, 3))))
  expect_equal(res$cs_df$cs_min_r2, c(0.81, 0.09))
  expect_equal(res$cs_df$cs_avg_r2, c(0.9025, 0.09))
  expect_equal(res$cs_df$low_purity, c(FALSE, TRUE))
})

test_that("variant table assigns credible sets and keeps all variants", {
  res <- extractResults(make_fit())

  expect_equal(nrow(res$variant_df), 4)
  expect_equal(res$variant_df$cs_id, c("L1", "L1", NA, "L2"))
  expect_equal(res$variant_df$pip, c(0.6, 0.7, 0.5, 1))
  expect_equal(names(res$variant_df)[1:5],
               c("variant_id", "pip", "z", "posterior_mean", "posterior_sd"))
  expect_true(all(c("alpha1", "alpha3", "mean3", "sd3") %in% names(res$variant_df)))
})

test_that("a fit with no credible sets returns NULL tables", {
  fit <- make_fit()
  fit$sets <- list(cs = NULL, purity = NULL, cs_index = NULL)
  res <- extractResults(fit)

  expect_null(res$cs_df)
  expect_null(res$variant_df)
})
