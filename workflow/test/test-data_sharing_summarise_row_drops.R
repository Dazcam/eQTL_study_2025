# test-data_sharing_summarise_row_drops.R
# - Used in data_sharing_eQTL_perm.R, data_sharing_eQTL_nominal.R

drops <- tibble::tibble(cell_type = c("Glu-A", "NPC"), rows_in = c(10L, 5L),
                        dot_id = c(1L, 0L), no_pvar_match = c(2L, 1L),
                        no_gene_symbol = c(3L, 0L), duplicate_gene_snp = c(0L, 1L),
                        rows_out = c(7L, 4L))

test_that("a TOTAL row is appended with column sums", {
  suppressMessages(res <- summarise_row_drops(drops))
  expect_equal(res$cell_type, c("Glu-A", "NPC", "TOTAL"))
  expect_equal(res$rows_in, c(10L, 5L, 15L))
  expect_equal(res$no_pvar_match, c(2L, 1L, 3L))
  expect_equal(res$rows_out, c(7L, 4L, 11L))
  expect_equal(names(res), names(drops))
})

test_that("rows not accounted for are noted in the log", {
  bad <- drops
  bad$rows_out[1] <- 6L
  expect_message(summarise_row_drops(bad), "1 rows unaccounted")
})
