#--------------------------------------------------------------------------------------
#
#    test-repl_sc_filter_jang_artifacts.R
#      - Used in repl_sc_pi1_enrichment.R and repl_sc_beta_correlation.R
#        (function in repl_sc_functions.R)
#
#--------------------------------------------------------------------------------------

library(testthat)
suppressMessages({
  library(dplyr)
  library(tibble)
})

filter_quiet <- function(...) suppressMessages(filter_jang_artifacts(...))

toy <- tibble(variant_id = paste0("rs", 1:6),
              Fixed_P    = c(1,   1,   0.2, 1,  NA,  0.5),
              Random_P   = c(1,   0.3, 1,   NA, 1,   0.5))

test_that("only rows with Fixed_P == 1 and Random_P == 1 are removed", {
  out <- filter_quiet(toy)
  expect_false("rs1" %in% out$variant_id)
  expect_true(all(c("rs2", "rs3", "rs6") %in% out$variant_id))
})

test_that("rows where the test is NA are dropped, as dplyr::filter() does", {
  out <- filter_quiet(toy)
  expected <- dplyr::filter(toy, !(Fixed_P == 1 & Random_P == 1))
  expect_identical(out, expected)
})

test_that("data.table input gives the same rows as data.table's own subsetting", {
  skip_if_not_installed("data.table")
  dt <- data.table::as.data.table(toy)
  out <- filter_quiet(dt)
  expected <- dt[!(Fixed_P == 1 & Random_P == 1)]
  expect_identical(out, expected)
})
