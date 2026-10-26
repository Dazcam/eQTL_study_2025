#--------------------------------------------------------------------------------------
#
#    test-align_gwas_to_ref.R
#      - Used in ctwas_run.R (function in ctwas_functions.R)
#
#--------------------------------------------------------------------------------------

library(testthat)
suppressMessages({
  library(dplyr)
  library(tibble)
  library(stringr)
})

# Reference .bim alleles (A1 = ALT, A2 = REF)
ref_bim <- tribble(
  ~Chr, ~SNP,  ~CM, ~BP, ~A1, ~A2,
  1,    "rs1", 0,   100, "A", "G",
  1,    "rs2", 0,   200, "C", "T",
  1,    "rs3", 0,   300, "A", "C",
  1,    "rs4", 0,   400, "G", "T",
  1,    "rs6", 0,   600, "A", "T",
  1,    "rs7", 0,   700, "C", "G"
)

# rs1 aligned, rs2 and rs4 flipped, rs3 different alleles, rs5 not in reference,
# rs6 (A/T) and rs7 (G/C) strand-ambiguous
toy_gwas <- tribble(
  ~SNP,  ~CHR, ~A1, ~A2, ~Z,
  "rs1", "1",  "A", "G",  2.0,
  "rs2", "1",  "T", "C",  3.0,
  "rs3", "1",  "G", "T",  1.5,
  "rs4", "1",  "T", "G", -4.0,
  "rs5", "1",  "A", "C",  1.0,
  "rs6", "1",  "A", "T",  0.5,
  "rs7", "1",  "G", "C", -0.5
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

test_that("only aligned and flipped, unambiguous SNPs are kept", {
  out <- align_quiet(toy_gwas, ref_bim)
  expect_equal(sort(out$id), c("rs1", "rs2", "rs4"))
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

test_that("SNPs with different alleles are removed", {
  expect_false("rs3" %in% align_quiet(toy_gwas, ref_bim)$id)
})

test_that("SNPs absent from the reference are removed", {
  expect_false("rs5" %in% align_quiet(toy_gwas, ref_bim)$id)
})

test_that("strand-ambiguous SNPs are removed, whether aligned or flipped", {
  out <- align_quiet(toy_gwas, ref_bim)
  expect_false(any(c("rs6", "rs7") %in% out$id))
})

test_that("SNPs with a missing or blank value are dropped", {
  gwas_na <- toy_gwas |>
    mutate(Z  = if_else(SNP == "rs1", NA_real_, Z),
           A2 = if_else(SNP == "rs2", "", A2))
  out <- align_quiet(gwas_na, ref_bim)
  expect_equal(out$id, "rs4")
})

test_that("the log reports removals and the final SNP count", {
  invisible(capture.output({
    expect_message(align_gwas_to_ref(toy_gwas, ref_bim),
                   "Removing 2 SNPs not in the reference or with different alleles")
    expect_message(align_gwas_to_ref(toy_gwas, ref_bim),
                   "Removing 2 strand-ambiguous SNPs")
    expect_message(align_gwas_to_ref(toy_gwas, ref_bim),
                   "Final number of SNPs in GWAS: 3\\s*$")
  }))
})
