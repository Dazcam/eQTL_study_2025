#--------------------------------------------------------------------------------------
#
#    test-repl_sc_make_snp_gr.R
#      - Used in repl_sc_atac_enrichment.R (function in repl_sc_functions.R)
#
#--------------------------------------------------------------------------------------

library(testthat)
suppressMessages(library(tibble))

coords <- tibble(refsnp_id   = c("rs1", "rs2", "rs2", "rs3"),
                 chr_name    = c("chr1", "chr1", "chr2", "chr5"),
                 chrom_start = c(100, 200, 300, 400),
                 chrom_end   = c(100, 200, 300, 400))

test_that("only requested SNPs are kept, one range per SNP (first record)", {
  skip_if_not_installed("GenomicRanges")
  gr <- make_snp_gr(c("rs2", "rs3", "rs9"), coords)
  expect_equal(length(gr), 2)
  expect_equal(gr$refsnp_id, c("rs2", "rs3"))
  expect_equal(as.character(GenomicRanges::seqnames(gr)), c("chr1", "chr5"))
  expect_equal(GenomicRanges::start(gr), c(200, 400))
})
