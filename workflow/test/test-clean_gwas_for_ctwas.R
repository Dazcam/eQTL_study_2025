#--------------------------------------------------------------------------------------
#
#    test-clean_gwas_for_ctwas.R
#      - Used in ctwas_run.R (function in ctwas_functions.R)
#
#--------------------------------------------------------------------------------------

library(testthat)
suppressMessages({
  library(dplyr)
  library(tibble)
  library(stringr)
})

# Run without the function's log output (messages and the printed chromosome table)
clean_quiet <- function(gwas) {
  out <- NULL
  invisible(capture.output(out <- suppressMessages(clean_gwas_for_ctwas(gwas))))
  out
}

clean_gwas <- tribble(
  ~SNP,  ~CHR, ~A1, ~A2, ~Z,
  "rs1", "1",  "A", "G", 1.0,
  "rs2", "2",  "C", "T", 2.0
)

test_that("a clean GWAS is returned unchanged", {
  expect_equal(clean_quiet(clean_gwas), clean_gwas)
})

test_that("non-rsID SNPs are removed", {
  gwas <- bind_rows(clean_gwas,
                    tibble(SNP = "1:12345:A:G", CHR = "1", A1 = "A", A2 = "G", Z = 3.0))
  expect_equal(clean_quiet(gwas)$SNP, c("rs1", "rs2"))
})

test_that("duplicate rsIDs keep only the first occurrence", {
  gwas <- bind_rows(clean_gwas,
                    tibble(SNP = "rs1", CHR = "1", A1 = "A", A2 = "C", Z = 9.0))
  out <- clean_quiet(gwas)
  expect_equal(out$SNP, c("rs1", "rs2"))
  expect_equal(out$Z[out$SNP == "rs1"], 1.0)
})

test_that("SNPs on alt contigs are removed", {
  gwas <- bind_rows(clean_gwas,
                    tibble(SNP = "rs3", CHR = "7_KI270803v1_alt", A1 = "A", A2 = "G", Z = 4.0))
  expect_equal(clean_quiet(gwas)$SNP, c("rs1", "rs2"))
})

test_that("the log reports how many SNPs were removed for each reason", {
  gwas <- bind_rows(
    clean_gwas,
    tibble(SNP = c("1:12345:A:G", "rs1", "rs3"),
           CHR = c("1", "1", "7_KI270803v1_alt"),
           A1 = "A", A2 = "G", Z = 0)
  )
  invisible(capture.output({
    expect_message(clean_gwas_for_ctwas(gwas), "Found 1 non-rsID SNPs")
    expect_message(clean_gwas_for_ctwas(gwas), "Found 1 duplicate rsIDs")
    expect_message(clean_gwas_for_ctwas(gwas), "Count of SNPs with 'alt' chromosomes: 1")
  }))
})
