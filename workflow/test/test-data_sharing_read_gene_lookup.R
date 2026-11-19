# test-data_sharing_read_gene_lookup.R
# - Used in data_sharing_eQTL_perm.R, data_sharing_eQTL_nominal.R, data_sharing_susie_cred_sets.R

lookup_file <- tempfile(fileext = ".tsv")
writeLines(c("ensembl_gene_id\texternal_gene_name\tchromosome_name",
             "ENSG01\tGENEA\t1",
             "ENSG01\tGENEA\t1",
             "ENSG02\tGENEB\t2",
             "ENSG02\tGENEB2\t2",
             "ENSG03\t\t3"), lookup_file)

test_that("one row per Ensembl ID, first symbol kept, blank symbol is NA", {
  suppressMessages(res <- read_gene_lookup(lookup_file))
  expect_equal(names(res), c("ensembl_gene_id", "external_gene_name"))
  expect_equal(res$ensembl_gene_id, c("ENSG01", "ENSG02", "ENSG03"))
  expect_equal(res$external_gene_name, c("GENEA", "GENEB", NA))
})

test_that("Ensembl IDs with more than one symbol are counted in the log", {
  expect_message(read_gene_lookup(lookup_file), "Gene lookup: 1 Ensembl IDs")
})
