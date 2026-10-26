#--------------------------------------------------------------------------------------
#
#    test-to_plink_pheno.R
#
#--------------------------------------------------------------------------------------
#
# - Used in twas_prep_expr_for_FUSION.R and twas_prep_covar_for_FUSION.R
#
#--------------------------------------------------------------------------------------

library(testthat)
library(dplyr)
library(tibble)

## Toy data  --------------------------------------------------------------------------
# 2 genes x 3 samples; sample columns deliberately not in alphabetical order
expr_toy <- tibble(
  TargetID = c("ENSG01", "ENSG02"),
  S3 = c(0.3, -0.3),
  S1 = c(0.1, -0.1),
  S2 = c(0.2, -0.2)
)

# 2 covariates x 2 samples, as read by twas_prep_covar_for_FUSION.R
covar_toy <- tibble(
  covariate = c("PC1", "sex"),
  S1 = c(0.5, 1),
  S2 = c(-0.5, 0)
)

## Tests  -----------------------------------------------------------------------------
test_that("expression table is transposed to FID, IID, one column per gene", {
  out <- to_plink_pheno(expr_toy, id_col = "TargetID")
  expect_equal(colnames(out), c("FID", "IID", "ENSG01", "ENSG02"))
  expect_equal(nrow(out), 3)
  expect_true(all(out$FID == 0))
})

test_that("sample order follows the input column order", {
  out <- to_plink_pheno(expr_toy, id_col = "TargetID")
  expect_equal(out$IID, c("S3", "S1", "S2"))
})

test_that("values land in the right sample x gene cell", {
  out <- to_plink_pheno(expr_toy, id_col = "TargetID")
  expect_equal(out$ENSG01[out$IID == "S1"], 0.1)
  expect_equal(out$ENSG02[out$IID == "S3"], -0.3)
})

test_that("works with a different ID column name (covariates)", {
  out <- to_plink_pheno(covar_toy, id_col = "covariate")
  expect_equal(colnames(out), c("FID", "IID", "PC1", "sex"))
  expect_equal(out$sex, c(1, 0))
  expect_equal(out$PC1, c(0.5, -0.5))
})

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
