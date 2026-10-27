#--------------------------------------------------------------------------------------
#
#    test-smr_align_freq_to_ref.R
#      - Used in smr_format_gwas.R and smr_create_query.R (function in smr_functions.R)
#
#--------------------------------------------------------------------------------------

library(testthat)
suppressMessages({
  library(dplyr)
  library(tibble)
})

# Run without the function's log output (messages and the printed status table)
align_quiet <- function(...) {
  out <- NULL
  invisible(capture.output(out <- suppressMessages(align_freq_to_ref(...))))
  out
}

# Reference .frq as read in the scripts (all character except MAF)
frq <- tribble(
  ~CHR, ~SNP,  ~A1, ~A2, ~MAF, ~NCHROBS,
  "1",  "rs1", "A", "G", 0.1,  1006L,
  "1",  "rs2", "C", "T", 0.2,  1006L,
  "1",  "rs3", "G", "A", 0.3,  1006L,
  "2",  "rs4", "T", "C", 0.4,  1006L
)

# rs1 aligned, rs2 flipped, rs3 different alleles, rs4 on another chromosome in the
# GWAS (no CHR + SNP match), rs5 absent from the reference
gwas <- tribble(
  ~SNP,  ~CHR, ~A1, ~A2, ~b,
  "rs1", "1",  "A", "G", 0.1,
  "rs2", "1",  "T", "C", -0.2,
  "rs3", "1",  "C", "T", 0.3,
  "rs4", "1",  "T", "C", 0.4,
  "rs5", "1",  "A", "G", 0.5
)

test_that("status and frequency follow the A1 of the input table (GWAS mode)", {
  out <- align_quiet(gwas, frq, by = c("CHR", "SNP"), freq_col = "freq")
  expect_equal(out$allele_status,
               c("aligned", "flipped", "mismatched", "mismatched", "mismatched"))
  expect_equal(out$freq, c(0.1, 0.8, NA, NA, NA))
})

test_that("alleles and effects are not changed", {
  out <- align_quiet(gwas, frq, by = c("CHR", "SNP"), freq_col = "freq")
  expect_equal(out$A1, gwas$A1)
  expect_equal(out$A2, gwas$A2)
  expect_equal(out$b, gwas$b)
})

test_that("reference columns are dropped and row order is kept", {
  out <- align_quiet(gwas, frq, by = c("CHR", "SNP"), freq_col = "freq")
  expect_equal(names(out), c(names(gwas), "allele_status", "freq"))
  expect_equal(out$SNP, gwas$SNP)
})

test_that("keep_unmatched keeps the existing frequency for unmatched SNPs (eQTL mode)", {
  eqtl <- gwas |> select(SNP, A1, A2, b) |> mutate(Freq = c(0.15, 0.75, 0.35, 0.45, 0.55))
  out <- align_quiet(eqtl, frq, by = "SNP", freq_col = "Freq", keep_unmatched = TRUE)
  # rs4 matches on SNP alone here, so it is aligned
  expect_equal(out$allele_status,
               c("aligned", "flipped", "mismatched", "aligned", "mismatched"))
  expect_equal(out$Freq, c(0.1, 0.8, 0.35, 0.4, 0.55))
})

test_that("the log reports the status counts under the label", {
  invisible(capture.output(
    expect_message(align_freq_to_ref(gwas, frq, by = c("CHR", "SNP"), freq_col = "freq",
                                     label = "GWAS"),
                   "=== GWAS frequency alignment summary ===")
  ))
})

