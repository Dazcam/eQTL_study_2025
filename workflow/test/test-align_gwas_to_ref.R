#--------------------------------------------------------------------------------------
#
#    test-align_gwas_to_ref.R
#      - Used in ctwas_run.R (function in ctwas_functions.R)
#
#--------------------------------------------------------------------------------------

library(testthat)
suppressMessages(library(tidyverse))

# Reference .bim alleles (A1 = ALT, A2 = REF)
ref_bim <- tribble(
  ~Chr, ~SNP,  ~CM, ~BP, ~A1, ~A2,
  1,    "rs1", 0,   100, "A", "G",
  1,    "rs2", 0,   200, "C", "T",
  1,    "rs3", 0,   300, "A", "C",
  1,    "rs4", 0,   400, "G", "T"
)

# rs1 aligned, rs2 and rs4 flipped, rs3 different alleles, rs5 not in reference
toy_gwas <- tribble(
  ~SNP,  ~CHR, ~A1, ~A2, ~Z,
  "rs1", "1",  "A", "G",  2.0,
  "rs2", "1",  "T", "C",  3.0,
  "rs3", "1",  "G", "T",  1.5,
  "rs4", "1",  "T", "G", -4.0,
  "rs5", "1",  "A", "C",  1.0
)

# Run without the function's log output (messages and the printed summary table)
align_quiet <- function(gwas, bim) {
  out <- NULL
  invisible(capture.output(out <- suppressMessages(align_gwas_to_ref(gwas, bim))))
  out
}

test_that("output has the cTWAS z_snp columns", {
  out <- align_quiet(toy_gwas, ref_bim)
  expect_equal(names(out), c("id", "A1", "A2", "z"))
})

test_that("aligned SNPs are unchanged", {
  out <- align_quiet(toy_gwas, ref_bim) |> filter(id == "rs1")
  expect_equal(c(out$A1, out$A2, out$z), c("A", "G", "2"))
})

test_that("flipped SNPs take the reference alleles and the opposite sign of Z", {
  out <- align_quiet(toy_gwas, ref_bim) |> filter(id %in% c("rs2", "rs4")) |> arrange(id)
  expect_equal(out$A1, c("C", "G"))
  expect_equal(out$A2, c("T", "T"))
  expect_equal(out$z,  c(-3.0, 4.0))
})

test_that("no SNP ends up with A1 equal to A2 (allele swap regression)", {
  out <- align_quiet(toy_gwas, ref_bim)
  expect_false(any(out$A1 == out$A2))
})

test_that("SNPs with different alleles are kept unchanged", {
  out <- align_quiet(toy_gwas, ref_bim) |> filter(id == "rs3")
  expect_equal(c(out$A1, out$A2, out$z), c("G", "T", "1.5"))
})

test_that("SNPs absent from the reference are kept unchanged", {
  out <- align_quiet(toy_gwas, ref_bim) |> filter(id == "rs5")
  expect_equal(c(out$A1, out$A2, out$z), c("A", "C", "1"))
})

test_that("SNPs with a missing or blank value are dropped", {
  gwas_na <- toy_gwas |>
    mutate(Z  = if_else(SNP == "rs1", NA_real_, Z),
           A2 = if_else(SNP == "rs3", "", A2))
  out <- align_quiet(gwas_na, ref_bim)
  expect_equal(sort(out$id), c("rs2", "rs4", "rs5"))
})

test_that("the final SNP count message matches the number of rows returned", {
  invisible(capture.output(
    expect_message(align_gwas_to_ref(toy_gwas, ref_bim), "Final number of SNPs in GWAS: 5\\s*$")
  ))
})
