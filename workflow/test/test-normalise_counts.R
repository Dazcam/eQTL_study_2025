#--------------------------------------------------------------------------------------
#
#    test-normalise_counts.R
#      - Used in tensorqtl_prep_input_files.R
#
#--------------------------------------------------------------------------------------

library(testthat)

test_that("bryois removes genes with mean CPM < 1", {
  cnts <- matrix(
    c(0.5, 0.4, 0.5, 0.4, 0.5,   # g1: mean ~0.46 -> removed
      2, 2, 2, 2, 2,             # g2: mean 2 -> kept
      3, 3, 3, 3, 3,             # g3
      4, 4, 4, 4, 4,             # g4
      5, 5, 5, 5, 5,             # g5
      6, 6, 6, 6, 6),            # g6
    nrow = 6, byrow = TRUE,
    dimnames = list(paste0("g", 1:6), paste0("s", 1:5))
  )
  out <- normalise_counts(cnts, "bryois", cell_type = "test_celltype")
  expect_equal(rownames(out), c("g2", "g3", "g4", "g5", "g6"))
})

test_that("bryois logs how many genes were removed and for which cell type", {
  cnts <- matrix(
    c(0.5, 0.4, 0.5, 0.4, 0.5,
      2, 2, 2, 2, 2,
      3, 3, 3, 3, 3,
      4, 4, 4, 4, 4,
      5, 5, 5, 5, 5,
      6, 6, 6, 6, 6),
    nrow = 6, byrow = TRUE,
    dimnames = list(paste0("g", 1:6), paste0("s", 1:5))
  )
  expect_message(
    normalise_counts(cnts, "bryois", cell_type = "RG"),
    "Removed 1 genes with mean CPM < 1 for RG"
  )
})

test_that("fujita filters genes below log2(CPM+1) = 2.0 in every sample", {
  cnts <- matrix(
    c(0.1, 0.1, 0.1, 0.1, 0.1,   # g1: low in all samples -> removed
      10, 10, 10, 10, 10,        # g2..g6: above threshold -> kept
      11, 11, 11, 11, 11,
      12, 12, 12, 12, 12,
      13, 13, 13, 13, 13,
      14, 14, 14, 14, 14),
    nrow = 6, byrow = TRUE,
    dimnames = list(paste0("g", 1:6), paste0("s", 1:5))
  )
  out <- normalise_counts(cnts, "fujita", cell_type = "test_celltype")
  expect_false("g1" %in% rownames(out))
  expect_equal(rownames(out), paste0("g", 2:6))
})

test_that("quantile normalisation keeps fujita-style filtering and returns same gene set size", {
  cnts <- matrix(
    c(10, 12, 11, 13, 9,
      8,  9,  10, 11, 7,
      15, 14, 16, 13, 15,
      20, 21, 19, 22, 20,
      6,  7,  8,  6,  7),
    nrow = 5, byrow = TRUE,
    dimnames = list(paste0("g", 1:5), paste0("s", 1:5))
  )
  out <- normalise_counts(cnts, "quantile", cell_type = "test_celltype")
  expect_equal(nrow(out), 5)
})

test_that("xue drops genes with >=90% zeros and z-scores the remainder", {
  zero_heavy <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 5)   # 9/10 zero (90%) -> removed
  kept_pattern <- c(4, 5, 6, 5, 4, 5, 6, 5, 4, 5)  # no zeros -> kept
  cnts <- rbind(
    g1 = zero_heavy,
    g2 = kept_pattern, g3 = kept_pattern,
    g4 = kept_pattern, g5 = kept_pattern, g6 = kept_pattern
  )
  colnames(cnts) <- paste0("s", 1:10)

  out <- normalise_counts(cnts, "xue", cell_type = "test_celltype")
  expect_false("g1" %in% rownames(out))
  expect_equal(rownames(out), paste0("g", 2:6))
  # z-scored rows should have ~zero mean
  expect_equal(unname(rowMeans(out)), rep(0, nrow(out)), tolerance = 1e-8)
})

test_that("an unrecognised normalisation method raises an informative error", {
  cnts <- matrix(1:4, nrow = 2, dimnames = list(c("g1", "g2"), c("s1", "s2")))
  expect_error(normalise_counts(cnts, "made_up_method"), "Unknown norm_method")
})
