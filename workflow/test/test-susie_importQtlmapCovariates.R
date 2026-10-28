# test-susie_importQtlmapCovariates.R
# - Used in susie_run.R

library(testthat)

test_that("covariates are transposed to samples x covariates and constant ones dropped", {
  cov_file <- tempfile(fileext = ".txt")
  writeLines(c("SampleID\tS1\tS2\tS3",
               "genPC1\t0.1\t0.2\t0.3",
               "const\t1\t1\t1",
               "expPC1\t-1\t0\t2"), cov_file)

  cov <- importQtlmapCovariates(cov_file)

  expect_equal(dim(cov), c(3, 2))
  expect_equal(rownames(cov), c("S1", "S2", "S3"))
  expect_equal(colnames(cov), c("genPC1", "expPC1"))
  expect_equal(unname(cov["S3", ]), c(0.3, 2))
})
