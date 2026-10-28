#--------------------------------------------------------------------------------------
#
#    SuSiE fine-mapping of significant eGenes
#
#--------------------------------------------------------------------------------------
#
# Pipeline:  08SUSIE | Rule: run_susie (one job per cell type and gene batch)
#            Upstream:   prep_susie_input, prep_susie_gene_meta, vcf_to_dosage,
#                        05TENSORQTL (permutation output)
#            Downstream: merge_susie -> sort_susie
#
# Purpose:   - Select genes with tensorQTL qval < fdr that have gene metadata
#            - Keep this job's chunk of the gene list
#            - For each gene: residualise expression and cis genotypes (TSS +/-
#              cisdistance) on the covariates and fit SuSiE (L = 10)
#            - Write credible set and variant results
#
# Inputs:    expression_matrix  phenotype_id + one column per sample
#            phenotype_meta     phenotype_id, group_id, chromosome, phenotype_pos (TSS)
#            sample_meta        one-column sample list (sample_id)
#            phenotype_list     tensorQTL permutation output (phenotype_id, qval)
#            covariates         covariates in rows, samples in columns
#            genotype_matrix    tabix-indexed dosage file (CHROM, POS, REF, ALT, samples)
#
# Outputs:   out_cred_hp  variants in high-purity credible sets (min.abs.corr >= 0.5)
#            out_cred     one row per credible set, all purities
#            out_snp      all variants with PIP, posterior and per-effect results
#
# Notes:     Run with Rscript from a shell rule; messages go to the rule's log.
#            SuSiE estimates credible set purity from a random sample of 100 variants
#            when a set has more than 100, so the seed is fixed for reproducibility.
#            Adapted from the eQTL Catalogue qtlmap run_susie.R.
#
#--------------------------------------------------------------------------------------

message("\n\nRunning SuSiE fine-mapping ...")

##  Load packages, functions and variables  -------------------------------------------
suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
})

script_dir <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)))
source(file.path(script_dir, "susie_functions.R"))

option_list <- list(
  make_option("--expression_matrix", type = "character", help = "Expression matrix (phenotype_id + samples)"),
  make_option("--phenotype_meta", type = "character", help = "Gene metadata (TSS positions)"),
  make_option("--sample_meta", type = "character", help = "Sample list (sample_id)"),
  make_option("--phenotype_list", type = "character", help = "tensorQTL permutation output"),
  make_option("--fdr", type = "double", help = "qval threshold for genes to fine-map"),
  make_option("--covariates", type = "character", help = "Covariates file"),
  make_option("--genotype_matrix", type = "character", help = "Tabix-indexed dosage file"),
  make_option("--cisdistance", type = "integer", help = "Window either side of the TSS (bp)"),
  make_option("--chunk", type = "character", help = "Chunk to run, as 'index n_chunks'"),
  make_option("--out_cred_hp", type = "character", help = "Output: high-purity CS variants"),
  make_option("--out_cred", type = "character", help = "Output: credible sets"),
  make_option("--out_snp", type = "character", help = "Output: all variants"))
opt <- parse_args(OptionParser(option_list = option_list))

seed <- 1

# Print what each variable is set to
variables <- c(expression_matrix = opt$expression_matrix,
               phenotype_meta    = opt$phenotype_meta,
               sample_meta       = opt$sample_meta,
               phenotype_list    = opt$phenotype_list,
               fdr               = opt$fdr,
               covariates        = opt$covariates,
               genotype_matrix   = opt$genotype_matrix,
               cisdistance       = opt$cisdistance,
               chunk             = opt$chunk,
               out_cred_hp       = opt$out_cred_hp,
               out_cred          = opt$out_cred,
               out_snp           = opt$out_snp,
               seed              = seed)
message("\nVariables")
message("============================")
message(paste(sprintf("%-18s  %s", names(variables), variables), collapse = "\n"))
message("============================\n")

##  Load inputs  ------------------------------------------------------------------------
message("Loading expression, sample list, gene metadata and covariates ...")
expression_tbl  <- readr::read_tsv(opt$expression_matrix, show_col_types = FALSE)
sample_metadata <- utils::read.csv(opt$sample_meta, sep = "\t", stringsAsFactors = FALSE)
phenotype_meta  <- utils::read.csv(opt$phenotype_meta, sep = "\t", stringsAsFactors = FALSE)
covariates_matrix <- importQtlmapCovariates(opt$covariates)

# Gene x sample expression matrix, samples in sample-list order
samples <- intersect(sample_metadata$sample_id, colnames(expression_tbl))
expr_mat <- as.matrix(as.data.frame(expression_tbl)[, samples])
rownames(expr_mat) <- expression_tbl$phenotype_id
message(nrow(expr_mat), " genes and ", length(samples), " samples in expression matrix.")

##  Select genes  -----------------------------------------------------------------------
message("Selecting genes at FDR < ", opt$fdr, " ...")
phenotype_table <- readr::read_tsv(opt$phenotype_list, show_col_types = FALSE)
filtered_list <- dplyr::filter(phenotype_table, qval < opt$fdr)
phenotype_list <- dplyr::semi_join(phenotype_meta, filtered_list,
                                   by = c("group_id" = "phenotype_id")) |>
  dplyr::filter(phenotype_id %in% expression_tbl$phenotype_id)
message(nrow(filtered_list), " genes pass FDR; ", nrow(phenotype_list),
        " have gene metadata and expression.")

# This job's chunk of the gene list
chunk_vector <- strsplit(opt$chunk, split = " ") |> unlist() |> as.numeric()
chunk_id <- chunk_vector[1]
n_chunks <- chunk_vector[2]
selected_chunk <- splitIntoChunks(chunk_id, n_chunks, length(phenotype_list$phenotype_id))
selected_phenotypes <- phenotype_list$phenotype_id[selected_chunk] %>% setNames(as.list(.), .)
message(length(selected_phenotypes), " genes in chunk ", chunk_id, " of ", n_chunks, ".")

##  Fine-map  ---------------------------------------------------------------------------
set.seed(seed)
results <- purrr::map(selected_phenotypes, ~finemapPhenotype(., expr_mat, phenotype_meta,
                                                             opt$genotype_matrix,
                                                             covariates_matrix,
                                                             opt$cisdistance))

##  Collect results  --------------------------------------------------------------------
cis_distance <- opt$cisdistance
region_df <- dplyr::transmute(phenotype_list, phenotype_id,
                              finemapped_region = paste0("chr", chromosome, ":",
                                                         phenotype_pos - cis_distance, "-",
                                                         phenotype_pos + cis_distance))

res <- purrr::map(results, extractResults) |>
  purrr::transpose()

# All variants
variant_df <- purrr::map_df(res$variant_df, identity, .id = "phenotype_id")
if (nrow(variant_df) > 0) {
  variant_df <- variant_df |>
    dplyr::left_join(region_df, by = "phenotype_id") |>
    tidyr::separate(variant_id, c("chr", "pos", "ref", "alt"), sep = "_", remove = FALSE) |>
    dplyr::mutate(chr = stringr::str_remove_all(chr, "chr")) |>
    dplyr::mutate(cs_index = cs_id) |>
    dplyr::mutate(cs_id = paste(phenotype_id, cs_index, sep = "_"))
}

# Credible sets
cs_df <- purrr::map_df(res$cs_df, identity, .id = "phenotype_id")

if (nrow(cs_df) > 0) {
  cs_df <- dplyr::left_join(cs_df, region_df, by = "phenotype_id") |>
    dplyr::mutate(cs_index = cs_id) |>
    dplyr::mutate(cs_id = paste(phenotype_id, cs_index, sep = "_")) |>
    dplyr::select(phenotype_id, cs_id, cs_index, finemapped_region, cs_log10bf,
                  cs_avg_r2, cs_min_r2, cs_size, low_purity)

  # Variants in high-purity credible sets
  in_cs_variant_df <- dplyr::filter(variant_df, !is.na(cs_index) & !low_purity) |>
    dplyr::select(phenotype_id, variant_id, chr, pos, ref, alt, cs_id, cs_index,
                  finemapped_region, pip, z, cs_min_r2, cs_avg_r2, cs_size,
                  posterior_mean, posterior_sd, cs_log10bf)
} else {
  # Empty tables with the same columns
  in_cs_variant_df <- dplyr::tibble(
    phenotype_id = character(), variant_id = character(), chr = numeric(),
    pos = numeric(), ref = numeric(), alt = numeric(), cs_id = numeric(),
    cs_index = numeric(), finemapped_region = numeric(), pip = numeric(),
    z = numeric(), cs_min_r2 = numeric(), cs_avg_r2 = numeric(), cs_size = numeric(),
    posterior_mean = numeric(), posterior_sd = numeric(), cs_log10bf = numeric())

  cs_df <- dplyr::tibble(
    phenotype_id = numeric(), cs_id = numeric(), cs_index = numeric(),
    finemapped_region = numeric(), cs_log10bf = numeric(), cs_avg_r2 = numeric(),
    cs_min_r2 = numeric(), cs_size = numeric(), low_purity = numeric())
}

if (nrow(variant_df) > 0) {
  variant_df <- dplyr::select(variant_df, phenotype_id, variant_id, chr, pos, ref, alt,
                              cs_id, cs_index, low_purity, finemapped_region, pip, z,
                              posterior_mean, posterior_sd, alpha1:sd10)
}

##  Write outputs  ----------------------------------------------------------------------
message("Writing ", nrow(cs_df), " credible sets (", nrow(in_cs_variant_df),
        " high-purity CS variants) ...")
write.table(in_cs_variant_df, opt$out_cred_hp, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
write.table(cs_df, opt$out_cred, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
write.table(variant_df, opt$out_snp, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

message("Done.")

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
