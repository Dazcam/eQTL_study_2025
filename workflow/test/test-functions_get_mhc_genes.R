#--------------------------------------------------------------------------------------
#
#    test-functions_get_mhc_genes.R
#      - Used in smr_create_query.R and ctwas_run.R (function in functions.R)
#
#--------------------------------------------------------------------------------------

library(testthat)
suppressMessages({
  library(tibble)
})

# Broad xMHC window (hg38), as in config mhc_region
mhc_chr   <- "6"
mhc_start <- 25000000
mhc_end   <- 34000000

lookup <- tribble(
  ~ensembl_gene_id, ~external_gene_name, ~chromosome_name, ~start_position, ~end_position,
  "ENSG01",         "INSIDE",            "6",              26457904,        26476720,   # BTN2A1-like, inside
  "ENSG02",         "LEFT_EDGE",         "6",              24990000,        25010000,   # spans start
  "ENSG03",         "RIGHT_EDGE",        "6",              33990000,        34100000,   # spans end
  "ENSG04",         "BEFORE",            "6",              24000000,        24999999,   # ends 1 bp before
  "ENSG05",         "AFTER",             "6",              34000001,        34100000,   # starts 1 bp after
  "ENSG06",         "OTHER_CHR",         "5",              30000000,        30010000,   # right range, wrong chr
  "ENSG07",         "CHR_PREFIX",        "chr6",           30000000,        30010000,   # 'chr' prefix
  "ENSG08",         "REVERSED",          "6",              30010000,        30000000,   # start > end
  "ENSG09",         "TOUCH_START",       "6",              24000000,        25000000,   # ends exactly at start
  "ENSG01",         "INSIDE",            "6",              26457904,        26476720,   # duplicate row
  "ENSG10",         "DUP_ONE_IN",        "6",              10000000,        10010000,   # duplicate ID: one row out,
  "ENSG10",         "DUP_ONE_IN",        "6",              30000000,        30010000,   #   one row in
  "ENSG11",         "MISSING_POS",       "6",              NA,              NA          # no coordinates
)

test_that("genes overlapping the window are returned once each", {
  out <- get_mhc_genes(lookup, mhc_chr, mhc_start, mhc_end)
  expect_setequal(out, c("ENSG01", "ENSG02", "ENSG03", "ENSG07", "ENSG08", "ENSG09", "ENSG10"))
  expect_equal(anyDuplicated(out), 0)
})

test_that("genes outside the window, on another chromosome or without positions are excluded", {
  out <- get_mhc_genes(lookup, mhc_chr, mhc_start, mhc_end)
  expect_false(any(c("ENSG04", "ENSG05", "ENSG06", "ENSG11") %in% out))
})

test_that("the region chromosome may be given with a 'chr' prefix or as a number", {
  expect_setequal(get_mhc_genes(lookup, "chr6", mhc_start, mhc_end),
                  get_mhc_genes(lookup, mhc_chr, mhc_start, mhc_end))
  expect_setequal(get_mhc_genes(lookup, 6, mhc_start, mhc_end),
                  get_mhc_genes(lookup, mhc_chr, mhc_start, mhc_end))
})

test_that("character positions (as read with all-character col_types) work", {
  chr_lookup <- lookup
  chr_lookup$start_position <- as.character(chr_lookup$start_position)
  chr_lookup$end_position   <- as.character(chr_lookup$end_position)
  expect_setequal(get_mhc_genes(chr_lookup, mhc_chr, mhc_start, mhc_end),
                  get_mhc_genes(lookup, mhc_chr, mhc_start, mhc_end))
})

test_that("custom column names can be used", {
  renamed <- setNames(lookup, c("id", "symbol", "chr", "start", "end"))
  out <- get_mhc_genes(renamed, mhc_chr, mhc_start, mhc_end,
                       id_col = "id", chr_col = "chr", start_col = "start", end_col = "end")
  expect_setequal(out, get_mhc_genes(lookup, mhc_chr, mhc_start, mhc_end))
})

test_that("a missing column gives a clear error", {
  expect_error(get_mhc_genes(lookup[, -5], mhc_chr, mhc_start, mhc_end),
               "missing column\\(s\\): end_position")
})
