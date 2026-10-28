# test-susie_posterior.R
# - Used in susie_run.R

library(testthat)

fit <- list(
  alpha = rbind(c(0.5, 0.5), c(0.2, 0.8)),
  mu    = rbind(c(2, 4), c(1, 1)),
  mu2   = rbind(c(5, 17), c(2, 2)),
  V = c(1, 1),
  X_column_scale_factors = c(1, 2))

test_that("posterior mean and SD sum over retained effects", {
  # mean: colSums(alpha * mu) / scale = c(1 + 0.2, 2 + 0.8) / c(1, 2)
  expect_equal(susie_get_posterior_mean_custom(fit), c(1.2, 1.4))
  # var: colSums(alpha * mu2 - (alpha * mu)^2)
  #   = c(2.5 - 1 + 0.4 - 0.04, 8.5 - 4 + 1.6 - 0.64) = c(1.86, 5.46)
  expect_equal(susie_get_posterior_sd_custom(fit), sqrt(c(1.86, 5.46)) / c(1, 2))
})

test_that("a single retained effect works (the reason for the custom versions)", {
  fit1 <- fit
  fit1$V <- c(1, 0)
  expect_equal(susie_get_posterior_mean_custom(fit1), c(1, 2) / c(1, 2))
  expect_equal(susie_get_posterior_sd_custom(fit1), sqrt(c(2.5 - 1, 8.5 - 4)) / c(1, 2))
})

test_that("no retained effects gives zeros", {
  fit0 <- fit
  fit0$V <- c(0, 0)
  expect_equal(susie_get_posterior_mean_custom(fit0), c(0, 0))
  expect_equal(susie_get_posterior_sd_custom(fit0), c(0, 0))
})
