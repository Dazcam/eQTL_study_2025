#--------------------------------------------------------------------------------------
#
#    Prep input files for tensorQTL
#
#--------------------------------------------------------------------------------------

# Harmonises the pseudobulk and covariate data for each cell type into a TensorQTL 
# - TMM normalises pseudobulk data
# - Applies various additional normalisation options for eQTL sensitivity testing
# - Runs PCA to calculate expression PCs to add as covariates
# - Standardises gene IDs to Ensembl encoding
# - Output: normalised pseudobulk GeX counts & combined covariate file per cell-type

## Info  ------------------------------------------------------------------------------

## Set up logging for Snakemake
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

message("\n\nPrepping input files for tensorQTL input ...")

##  Load Packages, functions and variables  -------------------------------------------
# Install and load required libraries
library(edgeR)
library(tidyverse)
library(limma)  # For quantile normalization
library(sva)    # For ComBat

# Input and output paths
cov_file <- snakemake@input[["cov_file"]]
sex_file <- snakemake@input[["sex_file"]]
gene_lookup <- snakemake@input[["gene_lookup"]]
cov_out <- snakemake@output[["cov_out"]]
exp_out <- snakemake@output[["exp_out"]]
pseudoblk_dir <- snakemake@params[["pseudoblk_dir"]]
report_dir <- snakemake@params[["report_dir"]]
out_dir <- snakemake@params[["out_dir"]]
cell_type <- snakemake@wildcards[["cell_type"]]
norm_method <- snakemake@params[["norm_method"]] 
batch_var <- snakemake@params[["batch_var"]]

# For testing
# cov_file <- "../results/04GENOTYPES-POST/covariates/pca.eigenvec"
# sex_file <- "reports/02SCANPY/02_metadata_per_sample.tsv"
# gene_lookup <- "../resources/sheets/gene_lookup_hg38.tsv"
# cell_type <- "RG"
# cov_out <- paste0("../05TENSORQTL/prep_input/", cell_type, "_base_covariates.txt")
# exp_out <- paste0("../05TENSORQTL/prep_input/", cell_type, "_tmm.bed")
# pseudoblk_dir <- "../results/03SCANPY/pseudobulk/"
# out_dir <- "../results/05TENSORQTL/prep_input/"
# report_dir <- "../workflow/reports/05TENSORQTL/"
# norm_method <- c('fujita')
# batch_var <- c('quantile')

# Make a tibble showing what each variable is set to
message("\nVariables")
cat("============================")
tibble(
  variable = c("cov_file", "sex_file", "gene_lookup", "cov_out", "exp_out", 
               "pseudoblk_dir", "report_dir", "out_dir", "cell_type", "norm_method", "batch_var"),
  value    = c(cov_file, sex_file, gene_lookup, cov_out, exp_out, 
               pseudoblk_dir, report_dir, out_dir, cell_type, norm_method, batch_var)) |> 
  knitr::kable(format = "simple", align = "l") |>
  print()
message("\n============================\n")

## Load data  -------------------------------------------------------------------------
# Genotype covariate file generate in calc_genotype_pcs rule
genotype_cov_tbl <- read_delim(cov_file, delim = ' ', 
                               col_names = c('sample', 'sample2', paste0('genPC', seq(1,10,1)))) |>
  dplyr::select(-sample2)

# Sex and PCW sample info generated in scanpy_clustering script
meta_data <- read_tsv(sex_file) |>
  dplyr::select(sample, sex_code, PCW) |>
  filter(sex_code != 0) |>
  filter(PCW != 'unknown') 

cov_tbl <- meta_data |>
  left_join(genotype_cov_tbl)

# Load Biomart gene lookup tables
# Run once: hg38_lookup <- get_biomart_gene_lookup('hg38'); see workflow/scripts/functions.R
gene_lookup <- read_tsv(gene_lookup)

#Print message with formatted tables
message("Sample Counts")
message("============================\n")

message("Sample cnt in genotype covariate file: ", nrow(genotype_cov_tbl))
message("Sample cnt in metadata file: ", nrow(meta_data))

message("Counts by PCW:\n")
meta_data %>%
  group_by(PCW) %>%
  dplyr::count() %>%
  knitr::kable(format = "simple", align = "c") %>%
  print()

message("\nCounts by Sex:")
meta_data %>%
  group_by(sex_code) %>%
  dplyr::count() %>%
  knitr::kable(format = "simple", align = "c") %>%
  print()

message("\nNAs in cov table:", meta_data |> anyNA())
message("Sample in mrgd covariate tbl: ", nrow(cov_tbl))

message("\n============================\n")

## Generate cell-specific TMM normalised counts and covariates for QTL analyses  ------
report_list <- list()
#for (cell_type in cell_types) { # Only need loop for local testing

message('\n\nCreating normalised counts and covariate files for: ', cell_type, ' with norm method: ', norm_method)
report_tibble <- tibble(cell_type = cell_type, norm_method = norm_method)  

# Load counts: don't use dplyr as we want to keep rownames
message('\nLoading counts ... \n')
pseudblk_cnts <- read.csv(paste0(pseudoblk_dir, cell_type, "_pseudobulk.csv"), row.names = 1)
print(pseudblk_cnts[1:5, 1:5])

# Transpose to genes as rows, samples as columns (TensorQTL format)
message('\nTransposing counts ... \n')
pseudblk_cnts <- t(pseudblk_cnts)
message('\nDimensions of ', cell_type, ' pseudblk cnts counts: ', 
        paste0(dim(pseudblk_cnts)[1], ' x ', dim(pseudblk_cnts)[2]))
print(pseudblk_cnts[1:5, 1:5])

# Subset sample_lst to available samples
available_samples <- intersect(cov_tbl$sample, colnames(pseudblk_cnts))
if (length(available_samples) == 0) {
  stop('No overlapping samples between genotypes and expression for ', cell_type, 
       '. Check sample name formats.')
}
message('\nSamples in genotypes but not in expression for ', cell_type, ': ', 
        paste(setdiff(cov_tbl$sample, colnames(pseudblk_cnts)), collapse = ' '))

# Filter zero-sum genes (rows) and reorder samples (columns)
if (sum(rowSums(pseudblk_cnts) == 0) != 0) {
  message('Removed ', sum(rowSums(pseudblk_cnts) == 0), ' genes not expressed in cell type')
  pseudblk_cnts <- pseudblk_cnts[rowSums(pseudblk_cnts != 0) > 0, , drop = FALSE]
}

# Reorder samples to match available subset of sample_lst
pseudblk_cnts <- pseudblk_cnts[, match(available_samples, colnames(pseudblk_cnts)), drop = FALSE]

# TMM normalization (base for all methods)
# Create a DGEList object
message('\nTMM normalising counts ... \n')
dge <- DGEList(counts = pseudblk_cnts)
dge <- calcNormFactors(dge, method = "TMM") # Do we need two normalisation factors here?
normalised_cnts <- edgeR::cpm(dge, normalized.lib.sizes = TRUE)
print(normalised_cnts[1:5, 1:5])
message('\nDimensions of normalised counts: ', 
        paste0(dim(normalised_cnts)[1], ' x ', dim(normalised_cnts)[2]))

# Branch for different normalization methods and filters
message('\nApplying normalization method: ', norm_method, ' ...\n')
if (norm_method == "bryois") {
  # Bryois: Filter genes with mean CPM < 1; PMID:35915177
  mean_cpm <- rowMeans(normalised_cnts)
  low_expr_genes <- mean_cpm < 1
  if (sum(low_expr_genes) > 0) {
    message('Removed ', sum(low_expr_genes), ' genes with mean CPM < 1 for ', cell_type)
    normalised_cnts <- normalised_cnts[!low_expr_genes, , drop = FALSE]
  }
} else if (norm_method == "fujita") {
  # Fujita: log2(CPM +1), filter genes with log2 CPM < 2.0 in all samples; PMID:38514782
  normalised_cnts <- log2(normalised_cnts + 1)  # Add 1 to avoid log(0)
  print(normalised_cnts[1:5, 1:5])
  low_expr_genes <- rowSums(normalised_cnts >= 2.0) == 0  # Genes with log2 CPM < 2.0 in all samples
  if (sum(low_expr_genes) > 0) {
    message('Removed ', sum(low_expr_genes), ' genes with log2 CPM < 2.0 in all samples for ', cell_type)
    normalised_cnts <- normalised_cnts[!low_expr_genes, , drop = FALSE]
  }
} else if (norm_method == "quantile") {
  # Quantile normalization across samples with Bryois filter
  normalised_cnts <- normalizeQuantiles(normalised_cnts)
  mean_cpm <- rowMeans(normalised_cnts)
  low_expr_genes <- mean_cpm < 1
  if (sum(low_expr_genes) > 0) {
    message('Removed ', sum(low_expr_genes), ' genes with mean CPM < 1 post-quantile for ', cell_type)
    normalised_cnts <- normalised_cnts[!low_expr_genes, , drop = FALSE]
  }
} else if (norm_method == "combat") {
  # ComBat batch correction (using specified batch_var, e.g., PCW)
  batch <- cov_tbl[[batch_var]][match(available_samples, cov_tbl$sample)]  # Align to samples
  if (any(is.na(batch))) stop("NA in batch variable; check cov_tbl")
  normalised_cnts <- ComBat(dat = normalised_cnts, batch = batch, mod = NULL, par.prior = TRUE)
  # Optionally apply Bryois filter post-ComBat
  mean_cpm <- rowMeans(normalised_cnts)
  low_expr_genes <- mean_cpm < 1
  if (sum(low_expr_genes) > 0) {
    message('Removed ', sum(low_expr_genes), ' genes with mean CPM < 1 post-ComBat for ', cell_type)
    normalised_cnts <- normalised_cnts[!low_expr_genes, , drop = FALSE]
  }
} else if (norm_method == "Xue_opt11") {
  # Xue Option #11: QC (≥90% zeros), log1p, z-score (STD); PMID:36823676
  zero_prop <- rowMeans(normalised_cnts == 0)
  low_expr_genes <- zero_prop >= 0.9
  if (sum(low_expr_genes) > 0) {
    message('Removed ', sum(low_expr_genes), ' genes with >=90% zeros for ', cell_type)
    normalised_cnts <- normalised_cnts[!low_expr_genes, , drop = FALSE]
  }
  normalised_cnts <- log(normalised_cnts + 1)  # log1p
  normalised_cnts <- t(scale(t(normalised_cnts)))  # Z-score per gene (mean 0, SD 1)
} else {
  stop("Unknown norm_method: ", norm_method)
}

message('\nDimensions after method-specific filtering: ',
        paste0(dim(normalised_cnts)[1], ' x ', dim(normalised_cnts)[2]))

# PCA Analysis
message('\nRunning PCA on expression values ...\n')
pca <- prcomp(t(normalised_cnts), scale. = TRUE)
exp_pc_scores <- as.data.frame(pca$x[, 1:50]) |> 
  as_tibble(rownames = 'sample') |>
  rename_with(~ paste0("exp", .), .cols = starts_with("PC"))

# Variance explained by expression PCs
exp_var_explained <- summary(pca)$importance[2, 1:50]  # Proportion of variance
report_tibble$exp_var_explained <- list(exp_var_explained)

# Combine exp pcs with genotype covariates
cov_full_tbl <- cov_tbl |> 
  inner_join(exp_pc_scores, by = "sample") |> 
  mutate(sample = factor(sample, levels = available_samples)) |> 
  arrange(match(sample, available_samples)) |> 
  dplyr::rename(id = sample) |> 
  mutate(PCW = as.numeric(PCW)) 

# Correlation between genotype and expression PCs
cor_matrix <- cor(cov_full_tbl[-1], use = "pairwise.complete.obs")
report_tibble$cor_matrix <- list(cor_matrix)

# Transpose and set row names
cov_matrix <- t(as.matrix(cov_full_tbl[, -1]))  # Exclude 'id' column, transpose
colnames(cov_matrix) <- cov_full_tbl$id         # Set sample IDs as column names
rownames(cov_matrix) <- colnames(cov_full_tbl)[-1]  # Set covariate names as row names
cov_matrix[1:5, 1:5]

# Write to file without column name for the index
message('\nWriting ', cell_type, ' covariate matrix for ', norm_method, ' ...')
write.table(cov_matrix, 
            file = paste0(out_dir, cell_type, "_", norm_method, "_base_covariates.txt"), 
            sep = "\t", 
            quote = FALSE, 
            col.names = NA,    # Write sample IDs as column names
            row.names = TRUE,    # Write covariate names as row names
            append = FALSE)


## Prepare gene expression data -----

# Gene expression data
# Add Ensembl IDs using BiomaRt Gene lookup; some genes are symbol, so are emsembl ID
message('\nHarmonising gene expression gene IDs to Ensembl ... ')
message('Rows in pseudobulk before merging with gene lookup: ', nrow(normalised_cnts))
pseudblk_ensembl_cnts <- normalised_cnts |>
  as_tibble(rownames = 'genes') |>
  filter(str_detect(genes, '^ENSG')) |>
  inner_join(gene_lookup, by = join_by(genes == ensembl_gene_id)) |>
  dplyr::select(-hgnc_symbol, -hgnc_id, -entrezgene_id, -external_gene_name) |> # Need to rm these col or all cols NA
  drop_na() |>
  mutate(
    TSS = if_else(strand == 1, start_position, end_position),  # Determine TSS based on strand
    cis_start = TSS,  # Prevent negative coordinates
    cis_end = TSS + 1 
  ) |>
  distinct(genes, .keep_all = TRUE) |> # Keep first occurrence dirty for now
  dplyr::select(Chr = chromosome_name, start = cis_start, end = cis_end, TargetID = genes, any_of(available_samples)) 
message('Rows annotated to Ensembl ID rather than symbol: ', nrow(pseudblk_ensembl_cnts))

pseudblk_cnts <- normalised_cnts |>
  as_tibble(rownames = 'genes') |>
  left_join(gene_lookup, by = join_by(genes == external_gene_name)) |>
  relocate(genes, ensembl_gene_id) 

message('Counts tbl dimensions after first merge on gene id (Genes x [Samples + 8 annotation cols]): ', 
        paste0(dim(pseudblk_cnts)[1], ' x ', dim(pseudblk_cnts)[2]))
report_tibble$initial_dims <- paste0(dim(pseudblk_cnts)[1], "x", dim(pseudblk_cnts)[2])

# Handle NAs and duplicates: scope here to try to retain more genes
pseudblk_cnts_nas <- pseudblk_cnts |>
  dplyr::filter(if_any(everything(), is.na)) 
message('Counts tbl gene ID NAs (note some : ', dim(pseudblk_cnts_nas)[1])
report_tibble$NAs <- dim(pseudblk_cnts_nas)[1]

pseudblk_cnts_symbol_dups <- pseudblk_cnts |>
  drop_na() |>
  group_by(genes) %>%
  dplyr::filter(n() > 1) %>%
  ungroup()
message('Counts tbl gene symbol dups: ', length(pseudblk_cnts_symbol_dups$genes |> unique()))
report_tibble$symbol_dups <- length(pseudblk_cnts_symbol_dups$genes |> unique())

pseudblk_cnts_ens_dups <- pseudblk_cnts |>
  drop_na() |>
  group_by(ensembl_gene_id) %>%
  dplyr::filter(n() > 1) %>%
  ungroup()
message('Counts tbl gene ensembl dups: ', length(pseudblk_cnts_ens_dups$genes |> unique()))
report_tibble$ensembl_dups <- length(pseudblk_cnts_ens_dups$genes |> unique())

# Double check: Drop genes with zero expression across all samples
pseudblk_cnts <- pseudblk_cnts %>%
  filter(rowSums(across(all_of(available_samples))) > 0)

# Final BED file
message('Writing final normalised gene expression matrix for ', norm_method, ' ... ')
pseudblk_cnts <- pseudblk_cnts |>
  drop_na() |>
  mutate(
    TSS = if_else(strand == 1, start_position, end_position),  # Determine TSS based on strand
    cis_start = TSS,  # Prevent negative coordinates
    cis_end = TSS + 1 
  ) |>
  distinct(ensembl_gene_id, .keep_all = TRUE) |> # Keep first occurrence dirty for now
  dplyr::select(Chr = chromosome_name, start = cis_start, end = cis_end, TargetID = ensembl_gene_id, any_of(available_samples)) |>
  rbind(pseudblk_ensembl_cnts) |> # Salvage genes marked as ENSEMBL ID, rather than symbol initially
  distinct(TargetID, .keep_all = TRUE) |>
  arrange(Chr, as.numeric(start), as.numeric(end)) |>
  dplyr::filter(Chr %in% seq(1,22,1)) |>
  dplyr::rename('#Chr' = Chr) |> # Required or tabix chokes at tensorQTL step
  write_tsv(paste0(out_dir, cell_type, '_', norm_method, '.bed'))

message('Final counts tbl dims: ', paste0(dim(pseudblk_cnts)[1], ' x ', dim(pseudblk_cnts)[2]), '\n')
report_tibble$final_dims <- paste0(dim(pseudblk_cnts)[1], 'x', dim(pseudblk_cnts)[2])
report_list[[paste(cell_type, norm_method, sep="_")]] <- report_tibble  # Updated key for multi-norm

# Combine all reports into a single tibble
final_report <- bind_rows(report_list) |>
  print()

write_rds(final_report, paste0(report_dir, 'gene_exp_matrix_report_', cell_type, '_', norm_method, '.rds'))


#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------