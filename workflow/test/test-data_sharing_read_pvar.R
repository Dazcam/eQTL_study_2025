# test-data_sharing_read_pvar.R
# - Used in data_sharing_eQTL_perm.R, data_sharing_eQTL_nominal.R

pvar_file <- tempfile(fileext = ".pvar")
writeLines(c("##fileformat=VCFv4.2",
             "#CHROM\tPOS\tID\tREF\tALT\tINFO",
             "1\t100\trs1\tA\tG\t.",
             "1\t200\t.\tC\tT\t.",
             "1\t300\trs2\tG\tA\t.",
             "1\t301\trs2\tG\tC\t.",
             "2\t400\trs3\tT\tC\t."), pvar_file)

test_that("dot IDs dropped, first row kept per duplicated ID", {
  suppressMessages(res <- read_pvar(pvar_file))
  expect_equal(names(res$pvar), c("CHROM", "POS", "ID", "REF", "ALT"))
  expect_equal(res$pvar$ID, c("rs1", "rs2", "rs3"))
  expect_equal(res$pvar$ALT, c("G", "A", "C"))
})

test_that("duplicated SNP IDs are reported with their counts", {
  suppressMessages(res <- read_pvar(pvar_file))
  expect_equal(res$dup_snps$ID, "rs2")
  expect_equal(res$dup_snps$n, 2L)
})
