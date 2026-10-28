# test-susie_splitIntoChunks.R
# - Used in susie_run.R

library(testthat)

test_that("chunks partition the gene list, remainder in the last chunk", {
  chunks <- lapply(1:3, function(i) which(splitIntoChunks(i, 3, 10)))
  expect_equal(chunks, list(1:3, 4:6, 7:10))
  expect_equal(sort(unlist(chunks)), 1:10)
})

test_that("a single chunk selects everything", {
  expect_true(all(splitIntoChunks(1, 1, 7)))
})

test_that("splitIntoBatches gives consecutive batch IDs", {
  expect_equal(splitIntoBatches(7, 3), c(1, 1, 1, 2, 2, 2, 3))
})
