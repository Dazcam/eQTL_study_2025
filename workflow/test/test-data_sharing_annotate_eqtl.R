# test-data_sharing_annotate_eqtl.R
# - Used in data_sharing_eQTL_perm.R, data_sharing_eQTL_nominal.R

pvar <- tibble::tibble(CHROM = c("1", "chr2", "1"), POS = c(100, 400, 300),
                       ID = c("rs1", "rs3", "rs2"), REF = c("A", "T", "G"),
                       ALT = c("G", "C", "A"))
gene_lookup <- tibble::tibble(ensembl_gene_id = c("ENSG01", "ENSG02"),
                              external_gene_name = c("GENEA", "GENEB"))
eqtl <- tibble::tibble(ensembl_id   = c("ENSG01", "ENSG01", "ENSG02", "ENSG03", "ENSG02", "ENSG02"),
                       SNP          = c("rs1", ".", "rs3", "rs2", "rs9", "rs3"),
                       af           = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.3),
                       pval_nominal = c(1e-5, 0.1, 0.01, 0.2, 0.3, 0.01),
                       slope        = c(0.5, 0.1, -0.2, 0.3, 0.1, -0.2),
                       slope_se     = c(0.1, 0.1, 0.1, 0.1, 0.1, 0.1),
                       qval         = c(0.01, 0.5, 0.05, 0.6, 0.7, 0.05))

run <- function(...) suppressMessages(annotate_eqtl(...))

test_that("dot IDs and SNPs absent from the pvar are removed, in input order", {
  res <- run(eqtl, pvar, gene_lookup)
  expect_equal(res$tbl$SNP, c("rs1", "rs3", "rs2", "rs3"))
  expect_equal(res$counts$rows_in, 6L)
  expect_equal(res$counts$dot_id, 1L)
  expect_equal(res$counts$no_pvar_match, 1L)
  expect_equal(res$counts$rows_out, 4L)
})

test_that("output columns, chr prefix removed, AF renamed", {
  res <- run(eqtl, pvar, gene_lookup)
  expect_equal(names(res$tbl), c("ensembl_id", "symbol", "CHROM", "SNP", "POS", "REF",
                                 "ALT", "AF", "slope", "slope_se", "pval_nominal"))
  expect_equal(res$tbl$CHROM, c("1", "2", "1", "2"))
  expect_equal(res$tbl$AF, c(0.1, 0.3, 0.4, 0.3))
})

test_that("missing gene symbols and duplicate gene-SNP pairs are kept and counted", {
  res <- run(eqtl, pvar, gene_lookup)
  expect_equal(res$tbl$symbol, c("GENEA", "GENEB", NA, "GENEB"))
  expect_equal(res$counts$no_gene_symbol, 1L)
  expect_equal(res$counts$duplicate_gene_snp, 1L)
})

test_that("extra columns are appended after pval_nominal", {
  res <- run(eqtl, pvar, gene_lookup, extra_cols = "qval")
  expect_equal(tail(names(res$tbl), 2), c("pval_nominal", "qval"))
  expect_equal(res$tbl$qval, c(0.01, 0.05, 0.6, 0.05))
})

test_that("an empty table gives an empty output with zero counts", {
  res <- run(eqtl[0, ], pvar, gene_lookup)
  expect_equal(nrow(res$tbl), 0L)
  expect_equal(unlist(res$counts), c(rows_in = 0, dot_id = 0, no_pvar_match = 0,
                                     no_gene_symbol = 0, duplicate_gene_snp = 0, rows_out = 0))
})
