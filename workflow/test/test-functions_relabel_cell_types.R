# test-functions_relabel_cell_types.R
# - Used in 16 tables and 17 figures scripts

labels <- list("Glu-UL" = "Glu-A", "Glu-DL" = "Glu-B", "-Q4-" = "-")

test_that("L1 and L2 cell types are relabelled", {
  expect_equal(relabel_cell_types(c("Glu-UL", "Glu-DL", "Glu-UL-0", "Glu-DL-2", "GABA", "NPC-1"), labels),
               c("Glu-A", "Glu-B", "Glu-A-0", "Glu-B-2", "GABA", "NPC-1"))
})

test_that("pseudotime bins lose -Q4-", {
  expect_equal(relabel_cell_types(c("NPC-to-Glu-UL-Q4-Bin1", "NPC-to-Glu-DL-Q4-Bin4"), labels),
               c("NPC-to-Glu-A-Bin1", "NPC-to-Glu-B-Bin4"))
})

test_that("strings containing cell types and repeated names are relabelled", {
  expect_equal(relabel_cell_types(c("Glu-UL_vs_GABA", "Glu-UL_vs_Glu-DL"), labels),
               c("Glu-A_vs_GABA", "Glu-A_vs_Glu-B"))
})

test_that("other labels, NA and factors are handled", {
  expect_equal(relabel_cell_types(c("Ext", NA, "Endo-Peri"), labels), c("Ext", NA, "Endo-Peri"))
  expect_equal(relabel_cell_types(factor(c("Glu-UL", "MG")), labels), c("Glu-A", "MG"))
  expect_equal(relabel_cell_types(character(0), labels), character(0))
})

test_that("a named character vector works like a list", {
  expect_equal(relabel_cell_types("Glu-DL-1", c("Glu-DL" = "Glu-B")), "Glu-B-1")
})

test_that("unnamed labels stop with an error", {
  expect_error(relabel_cell_types("Glu-UL", c("Glu-A")), "named")
})
