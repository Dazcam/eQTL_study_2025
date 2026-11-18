# test-functions_cell_type_levels.R
# - Used in 16 tables and 17 figures scripts

cell_types <- c("Glu-UL", "Glu-DL", "GABA", "NPC", "OPC", "MG", "Endo-Peri",
                "Glu-UL-0", "Glu-UL-1", "Glu-DL-0", "GABA-0", "NPC-0", "NPC-1",
                "NPC-to-Glu-UL-Q4-Bin1", "NPC-to-Glu-UL-Q4-Bin2")
labels <- list("Glu-UL" = "Glu-A", "Glu-DL" = "Glu-B", "-Q4-" = "-")

test_that("family layout puts each L1 type before its subclusters, then pseudotime", {
  x <- rev(cell_types)
  expect_equal(cell_type_levels(x, cell_types),
               c("Glu-UL", "Glu-UL-0", "Glu-UL-1", "Glu-DL", "Glu-DL-0", "GABA", "GABA-0",
                 "NPC", "NPC-0", "NPC-1", "OPC", "MG", "Endo-Peri",
                 "NPC-to-Glu-UL-Q4-Bin1", "NPC-to-Glu-UL-Q4-Bin2"))
})

test_that("level layout gives L1, then L2, then pseudotime", {
  x <- rev(cell_types)
  expect_equal(cell_type_levels(x, cell_types, layout = "level"),
               c("Glu-UL", "Glu-DL", "GABA", "NPC", "OPC", "MG", "Endo-Peri",
                 "Glu-UL-0", "Glu-UL-1", "Glu-DL-0", "GABA-0", "NPC-0", "NPC-1",
                 "NPC-to-Glu-UL-Q4-Bin1", "NPC-to-Glu-UL-Q4-Bin2"))
})

test_that("manuscript labels match after relabelling, and only cell types in x are returned", {
  expect_equal(cell_type_levels(c("NPC-to-Glu-A-Bin2", "GABA", "Glu-A-1", "Glu-B", "Glu-A", "GABA"),
                                cell_types, labels),
               c("Glu-A", "Glu-A-1", "Glu-B", "GABA", "NPC-to-Glu-A-Bin2"))
  expect_equal(cell_type_levels(factor(c("MG", "OPC")), cell_types), c("OPC", "MG"))
})

test_that("Endo-Peri is not treated as a subcluster", {
  expect_equal(cell_type_levels(c("Endo-Peri", "MG"), cell_types, layout = "level"), c("MG", "Endo-Peri"))
})

test_that("a cell type not in config stops with an error", {
  expect_error(cell_type_levels(c("Glu-A", "NPC-3"), cell_types, labels), "NPC-3")
  expect_error(cell_type_levels("Glu-A", cell_types), "Glu-A")
})
