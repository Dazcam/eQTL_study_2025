#--------------------------------------------------------------------------------------
#
#    pseudotime_prep_tensorqtl.R
#
#    Prepare input files for tensorQTL from pseudotime quantile pseudobulk data.
#
#    Adapted from tensorqtl_prep_input_files.R for the pseudotime trajectory
#    analysis. Applies quantile normalisation only (no sensitivity analysis).
#
#    Pipeline:
#      1. Load pseudobulk counts for a single trajectory/quantile/bin
#      2. Filter to genotype-matched samples
#      3. TMM normalise then quantile normalise
#      4. Filter lowly expressed genes (log2 CPM < 2.0 in all samples)
#      5. Run PCA to generate expression PCs as covariates
#      6. Combine with genotype PCs, sex, and PCW covariates
#      7. Map gene symbols to Ensembl IDs
#      8. Write BED file and covariate matrix
#
#    Usage: called via Snakemake rule pseudotime_prep_tensorqtl
#
#--------------------------------------------------------------------------------------

## Set up logging  -------------------------------------------------------------------

if (exists("snakemake")) {
  log_smk <- function() {
    if (exists("snakemake") && length(snakemake@log) != 0) {
      log <- file(snakemake@log[[1]], open = "wt")
      sink(log, append = TRUE)
      sink(log, append = TRUE, type = "message")
    }
  }
  log_smk()
}

message("\n\nPrepping pseudotime quantile input files for tensorQTL ...")

## Load packages  --------------------------------------------------------------------

library(edgeR)
library(tidyverse)
library(limma)

## Parameters  ------------------------------------------------------------------------

cov_file     <- snakemake@input[["cov_file"]]
sex_file     <- snakemake@input[["sex_file"]]
gene_lookup  <- snakemake@input[["gene_lookup"]]
exp_out      <- snakemake@output[["exp_out"]]
cov_out      <- snakemake@output[["cov_out"]]
pseudoblk_dir <- snakemake@params[["pseudoblk_dir"]]
out_dir      <- snakemake@params[["out_dir"]]
trajectory   <- snakemake@wildcards[["trajectory"]]
quantile_n   <- snakemake@wildcards[["quantile_n"]]
bin          <- snakemake@wildcards[["bin"]]

# Log variables
message("\nVariables")
message("============================")
tibble(
  variable = c("cov_file", "sex_file", "gene_lookup", "exp_out", "cov_out",
               "pseudoblk_dir", "out_dir", "trajectory", "quantile_n", "bin"),
  value    = c(cov_file, sex_file, gene_lookup, exp_out, cov_out,
               pseudoblk_dir, out_dir, trajectory, quantile_n, bin)) |>
  knitr::kable(format = "simple", align = "l") |>
  print()
message("============================\n")

## Load covariate data  ---------------------------------------------------------------

# Genotype PCs — fixed at 4 for pseudotime analysis
genotype_cov_tbl <- read_delim(cov_file, delim = ' ',
                               col_names = c('sample', 'sample2',
                                             paste0('genPC', seq(1, 10, 1)))) |>
  dplyr::select(-sample2)

# Sex and PCW metadata
meta_data <- read_tsv(sex_file) |>
  dplyr::select(sample, sex_code, PCW) |>
  filter(sex_code != 0) |>
  filter(PCW != 'unknown')

cov_tbl <- meta_data |>
  left_join(genotype_cov_tbl)

message("Sample counts")
message("============================")
message("Samples in genotype covariate file: ", nrow(genotype_cov_tbl))
message("Samples in metadata file:           ", nrow(meta_data))
message("Samples in merged covariate table:  ", nrow(cov_tbl))
message("NAs in covariate table:             ", anyNA(cov_tbl))
message("============================\n")

## Load gene lookup  ------------------------------------------------------------------

gene_lookup_tbl <- read_tsv(gene_lookup)

## Load pseudobulk counts  ------------------------------------------------------------

pseudobulk_file <- file.path(
  pseudoblk_dir,
  paste0("Q", quantile_n),
  paste0("Q", quantile_n, "_bin", bin, "_", trajectory, "_pseudobulk.csv")
)

message("Loading pseudobulk counts from: ", pseudobulk_file)
pseudblk_cnts <- read.csv(pseudobulk_file, row.names = 1)
message("Pseudobulk dimensions (samples x genes): ",
        paste0(dim(pseudblk_cnts)[1], " x ", dim(pseudblk_cnts)[2]))
print(pseudblk_cnts[1:5, 1:5])

# Transpose to genes as rows, samples as columns (TensorQTL format)
message("\nTransposing to genes x samples ...")
pseudblk_cnts <- t(pseudblk_cnts)
message("Dimensions after transpose: ",
        paste0(dim(pseudblk_cnts)[1], " x ", dim(pseudblk_cnts)[2]))

## Filter to genotype-matched samples  ------------------------------------------------

available_samples <- intersect(cov_tbl$sample, colnames(pseudblk_cnts))
if (length(available_samples) == 0) {
  stop("No overlapping samples between genotypes and expression for ",
       trajectory, " Q", quantile_n, " bin ", bin,
       ". Check sample name formats.")
}

message("\nSample overlap:")
message("  Samples in pseudobulk:  ", ncol(pseudblk_cnts))
message("  Samples in covariates:  ", nrow(cov_tbl))
message("  Overlapping samples:    ", length(available_samples))
message("  Samples in genotypes but not expression: ",
        paste(setdiff(cov_tbl$sample, colnames(pseudblk_cnts)), collapse = " "))

# Filter zero-sum genes
if (sum(rowSums(pseudblk_cnts) == 0) != 0) {
  n_removed <- sum(rowSums(pseudblk_cnts) == 0)
  message("\nRemoved ", n_removed, " genes with zero counts across all samples")
  pseudblk_cnts <- pseudblk_cnts[rowSums(pseudblk_cnts != 0) > 0, , drop = FALSE]
}

# Reorder samples to match available subset
pseudblk_cnts <- pseudblk_cnts[, match(available_samples, colnames(pseudblk_cnts)),
                                drop = FALSE]

## Normalisation — TMM then quantile  -------------------------------------------------

# TMM normalisation
message("\nTMM normalising counts ...")
dge              <- DGEList(counts = pseudblk_cnts)
dge              <- calcNormFactors(dge, method = "TMM")
normalised_cnts  <- edgeR::cpm(dge, normalized.lib.sizes = TRUE)
message("Dimensions after TMM: ",
        paste0(dim(normalised_cnts)[1], " x ", dim(normalised_cnts)[2]))

# Log2 transform
normalised_cnts <- log2(normalised_cnts + 1)

# Filter lowly expressed genes: log2 CPM < 2.0 in all samples
low_expr_genes <- rowSums(normalised_cnts >= 2.0) == 0
if (sum(low_expr_genes) > 0) {
  message("Removed ", sum(low_expr_genes),
          " genes with log2 CPM < 2.0 in all samples")
  normalised_cnts <- normalised_cnts[!low_expr_genes, , drop = FALSE]
}

# Quantile normalisation
message("\nApplying quantile normalisation ...")
normalised_cnts <- normalizeQuantiles(normalised_cnts)
message("Dimensions after quantile normalisation: ",
        paste0(dim(normalised_cnts)[1], " x ", dim(normalised_cnts)[2]))
print(normalised_cnts[1:5, 1:5])

## PCA on expression values  ----------------------------------------------------------

message("\nRunning PCA on normalised expression ...")
pca          <- prcomp(t(normalised_cnts), scale. = TRUE)
exp_pc_scores <- as.data.frame(pca$x[, 1:50]) |>
  as_tibble(rownames = 'sample') |>
  rename_with(~ paste0("exp", .), .cols = starts_with("PC"))

# Log variance explained
exp_var_explained <- summary(pca)$importance[2, 1:20]
message("\nVariance explained by first 20 expression PCs:")
print(round(exp_var_explained, 4))

## Build covariate matrix  ------------------------------------------------------------

cov_full_tbl <- cov_tbl |>
  inner_join(exp_pc_scores, by = "sample") |>
  mutate(sample = factor(sample, levels = available_samples)) |>
  arrange(match(sample, available_samples)) |>
  dplyr::rename(id = sample) |>
  mutate(PCW = as.numeric(PCW))

# Correlation between genotype and expression PCs
cor_matrix <- cor(cov_full_tbl[-1], use = "pairwise.complete.obs")
message("\nCorrelation between genotype PCs and first 5 expression PCs:")
print(round(cor_matrix[paste0("genPC", 1:4), paste0("expPC", 1:5)], 3))

# Transpose to covariates as rows, samples as columns
cov_matrix          <- t(as.matrix(cov_full_tbl[, -1]))
colnames(cov_matrix) <- cov_full_tbl$id
rownames(cov_matrix) <- colnames(cov_full_tbl)[-1]

message("\nWriting covariate matrix to: ", cov_out)
dir.create(dirname(cov_out), recursive = TRUE, showWarnings = FALSE)
write.table(cov_matrix,
            file      = cov_out,
            sep       = "\t",
            quote     = FALSE,
            col.names = NA,
            row.names = TRUE,
            append    = FALSE)
message("Covariate matrix dimensions: ",
        paste0(nrow(cov_matrix), " covariates x ", ncol(cov_matrix), " samples"))

## Prepare gene expression BED file  --------------------------------------------------

message("\nHarmonising gene IDs to Ensembl ...")
message("Genes in normalised counts: ", nrow(normalised_cnts))

# Genes already in Ensembl format
pseudblk_ensembl_cnts <- normalised_cnts |>
  as_tibble(rownames = 'genes') |>
  filter(str_detect(genes, '^ENSG')) |>
  inner_join(gene_lookup_tbl, by = join_by(genes == ensembl_gene_id)) |>
  dplyr::select(-hgnc_symbol, -hgnc_id, -entrezgene_id, -external_gene_name) |>
  drop_na() |>
  mutate(
    TSS       = if_else(strand == 1, start_position, end_position),
    cis_start = TSS,
    cis_end   = TSS + 1
  ) |>
  distinct(genes, .keep_all = TRUE) |>
  dplyr::select(Chr = chromosome_name, start = cis_start, end = cis_end,
                TargetID = genes, any_of(available_samples))
message("Genes already in Ensembl format: ", nrow(pseudblk_ensembl_cnts))

# Genes in symbol format — map to Ensembl
pseudblk_symbol_cnts <- normalised_cnts |>
  as_tibble(rownames = 'genes') |>
  left_join(gene_lookup_tbl, by = join_by(genes == external_gene_name)) |>
  relocate(genes, ensembl_gene_id)

message("Genes mapped from symbol to Ensembl: ", nrow(pseudblk_symbol_cnts))
message("NAs in symbol mapping: ",
        sum(is.na(pseudblk_symbol_cnts$ensembl_gene_id)))
message("Symbol duplicates: ",
        length(unique(pseudblk_symbol_cnts$genes[duplicated(pseudblk_symbol_cnts$genes)])))
message("Ensembl duplicates: ",
        length(unique(pseudblk_symbol_cnts$ensembl_gene_id[
          duplicated(pseudblk_symbol_cnts$ensembl_gene_id) &
          !is.na(pseudblk_symbol_cnts$ensembl_gene_id)])))

# Drop zero-expression genes
pseudblk_symbol_cnts <- pseudblk_symbol_cnts |>
  filter(rowSums(across(all_of(available_samples))) > 0)

# Final BED file
message("\nWriting final normalised BED file to: ", exp_out)
dir.create(dirname(exp_out), recursive = TRUE, showWarnings = FALSE)

final_bed <- pseudblk_symbol_cnts |>
  drop_na() |>
  mutate(
    TSS       = if_else(strand == 1, start_position, end_position),
    cis_start = TSS,
    cis_end   = TSS + 1
  ) |>
  distinct(ensembl_gene_id, .keep_all = TRUE) |>
  dplyr::select(Chr = chromosome_name, start = cis_start, end = cis_end,
                TargetID = ensembl_gene_id, any_of(available_samples)) |>
  rbind(pseudblk_ensembl_cnts) |>
  distinct(TargetID, .keep_all = TRUE) |>
  arrange(Chr, as.numeric(start), as.numeric(end)) |>
  filter(Chr %in% seq(1, 22, 1)) |>
  dplyr::rename('#Chr' = Chr)

write_tsv(final_bed, exp_out)

message("Final BED dimensions: ",
        paste0(nrow(final_bed), " genes x ", ncol(final_bed) - 4, " samples"))

message("\nDone.")

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
