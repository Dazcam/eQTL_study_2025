#--------------------------------------------------------------------------------------
#
#    test-eqtl_sharing_get_parent_cell_type.R
#      - Used in eqtl_sharing_extract_universe.R and eqtl_sharing_effect_sizes.R
#        (function in eqtl_sharing_functions.R)
#
#--------------------------------------------------------------------------------------

library(testthat)

cell_types_L1 <- c("Glu-UL", "Glu-DL", "GABA", "NPC", "MG", "OPC", "Endo-Peri")

test_that("L1 cell types return themselves, including hyphenated names", {
  expect_equal(get_parent_cell_type("Glu-UL", cell_types_L1), "Glu-UL")
  expect_equal(get_parent_cell_type("Endo-Peri", cell_types_L1), "Endo-Peri")
})

test_that("subclusters return their L1 parent", {
  expect_equal(get_parent_cell_type("Glu-UL-0", cell_types_L1), "Glu-UL")
  expect_equal(get_parent_cell_type("GABA-2", cell_types_L1), "GABA")
  expect_equal(get_parent_cell_type("NPC-1", cell_types_L1), "NPC")
})

test_that("works with sapply over a vector, as in eqtl_sharing_effect_sizes.R", {
  expect_equal(unname(sapply(c("GABA-0", "MG", "Glu-DL-2"), get_parent_cell_type,
                             cell_types_L1 = cell_types_L1)),
               c("GABA", "MG", "Glu-DL"))
})
