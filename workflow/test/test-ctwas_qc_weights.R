#--------------------------------------------------------------------------------------
#
#    test-qc_ctwas_weights.R
#      - Used in ctwas_run.R when run_weight_qc is on (function in ctwas_functions.R)
#
#--------------------------------------------------------------------------------------

library(testthat)
suppressMessages(library(tibble))

z_snp <- tibble(id = c("rs1", "rs2", "rs3"), z = c(2, -1, 0.5))

# A valid gene: two SNPs, weights as a 1-column matrix, LD as a 2 x 2 matrix
good_gene <- list(
  wgt   = matrix(c(0.5, 0.2), ncol = 1, dimnames = list(c("rs1", "rs2"), NULL)),
  R_wgt = matrix(c(1, 0.3, 0.3, 1), nrow = 2)
)

# Run without the per-gene progress lines
qc_quiet <- function(weights) {
  out <- NULL
  invisible(capture.output(out <- qc_ctwas_weights(weights, z_snp)))
  out
}

reason_for <- function(gene) unname(qc_quiet(list(g = gene))$reason)

test_that("returns one row per gene with gene, status and reason", {
  out <- qc_quiet(list(g1 = good_gene, g2 = good_gene))
  expect_equal(names(out), c("gene", "status", "reason"))
  expect_equal(out$gene, c("g1", "g2"))
})

test_that("a valid gene passes", {
  out <- qc_quiet(list(g = good_gene))
  expect_equal(unname(out$status), "PASS")
  expect_true(is.na(out$reason))
})

test_that("weights that are not a matrix fail", {
  expect_equal(reason_for(modifyList(good_gene, list(wgt = c(0.5, 0.2)))), "wgt_not_matrix")
})

test_that("missing LD fails", {
  gene <- good_gene
  gene["R_wgt"] <- list(NULL)
  expect_equal(reason_for(gene), "R_wgt_NULL")
})

test_that("LD that is not square fails", {
  expect_equal(reason_for(modifyList(good_gene, list(R_wgt = matrix(1, 2, 3)))), "R_wgt_not_square")
})

test_that("weights and LD with different numbers of SNPs fail", {
  expect_equal(reason_for(modifyList(good_gene, list(R_wgt = diag(3)))), "dim_mismatch")
})

test_that("a weight SNP absent from the GWAS fails", {
  gene <- good_gene
  rownames(gene$wgt) <- c("rs1", "rs99")
  expect_equal(reason_for(gene), "SNP_missing_in_z")
})

test_that("SNPs in perfect LD fail", {
  expect_equal(reason_for(modifyList(good_gene, list(R_wgt = matrix(1, 2, 2)))), "perfect_LD")
})

test_that("a gene whose z-score can't be computed fails", {
  # All-zero weights give 0 / 0 = NaN
  gene <- modifyList(good_gene, list(wgt = matrix(0, 2, 1, dimnames = list(c("rs1", "rs2"), NULL))))
  expect_match(reason_for(gene), "^compute_error")
})
