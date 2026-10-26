#--------------------------------------------------------------------------------------
#
#    test-make_ld_metatable.R
#      - Used in ctwas_run.R (function in ctwas_functions.R)
#
#--------------------------------------------------------------------------------------

library(testthat)

region_info <- data.frame(
  chrom = c(1L, 22L),
  start = c(10583L, 16000000L),
  stop  = c(1892607L, 17000000L)
)

meta_quiet <- function(...) suppressMessages(make_ld_metatable(...))

test_that("one row per region, with the original columns kept", {
  out <- meta_quiet(region_info, "ld_mats", prefix = "1000G_EUR_Phase3_hg38")
  expect_equal(nrow(out), 2)
  expect_equal(out[, c("chrom", "start", "stop")], region_info)
})

test_that("LD and SNP file paths match the names written by create_ld_matrices", {
  out <- meta_quiet(region_info, "ld_mats", prefix = "1000G_EUR_Phase3_hg38")
  expect_equal(out$LD_file,
               c("ld_mats/1000G_EUR_Phase3_hg38_chr1.R_snp.10583_1892607.RDS",
                 "ld_mats/1000G_EUR_Phase3_hg38_chr22.R_snp.16000000_17000000.RDS"))
  expect_equal(out$SNP_file,
               c("ld_mats/1000G_EUR_Phase3_hg38_chr1.R_snp.10583_1892607.Rvar",
                 "ld_mats/1000G_EUR_Phase3_hg38_chr22.R_snp.16000000_17000000.Rvar"))
})

test_that("the prefix argument is used in the file names", {
  out <- meta_quiet(region_info[1, ], "ld_mats", prefix = "other_ref")
  expect_equal(out$LD_file, "ld_mats/other_ref_chr1.R_snp.10583_1892607.RDS")
})
