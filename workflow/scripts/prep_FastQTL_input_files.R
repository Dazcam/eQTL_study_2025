#--------------------------------------------------------------------------------------
#
#    Prep FastQTL input files - draft
#
#--------------------------------------------------------------------------------------

## Info  ------------------------------------------------------------------------------

#  Prep input bed and covarite files for FastQTL 
#  Running locally for now

# Need to add loop for all cell types

##  Load Packages, functions and variables  -------------------------------------------
# Install and load required libraries
library(edgeR)
library(variancePartition)
library(tidyverse)
#library(peer) # Old package: non-standard installation

# Load data
data_dir <- "~/Desktop/eQTL_study_2025/results/03SCANPY/pseudobulk/"
cov_dir <- "~/Desktop/eQTL_study_2025/resources/test_genotypes/"
lookup_dir <- "~/Desktop/eQTL_study_2025/results/02PARSE/combine_plate1/all-sample/DGE_filtered/"


pseudobulk_counts <- read.csv(paste0(data_dir, "cell_type_1_pseudobulk.csv"), row.names = 1)  # rows=genes, columns=samples
pseudobulk_counts[1:5, 1:5]

cov_orig_tbl <- read_tsv(paste0(cov_dir, 'filtered_covariates.txt'), col_types = 'ccd') 
cov_tbl <- cov_orig_tbl |> 
  dplyr::select(Sample, Sex, PCW, PC1, PC2, PC3) |>
  rename_with(~ paste0("gen", .), .cols = matches("^PC\\d+$"))

gene_lookup <- read_csv(paste0(lookup_dir, 'all_genes.csv'), col_types = 'ccc')

# Need to order samples in input files in same order
overlap_lst <- read_tsv(paste0(cov_dir, 'overlapping_samples.txt'), col_types = 'c', col_names = 'Sample') 

window <- 250000

# TMM normalization
# Create a DGEList object
dge <- DGEList(counts = pseudobulk_counts)
dge <- calcNormFactors(dge, method = "TMM") # Do we need two normalisation factors here?
normalized_counts <- cpm(dge, normalized.lib.sizes = TRUE)
normalized_counts[1:5, 1:5]

# PCA Analysis
pca <- prcomp(normalized_counts, scale. = TRUE)
pc_scores <- as.data.frame(pca$x[, 1:10]) |> # Extract first 10 PCs
  as_tibble(rownames = 'Sample') |>
  rename_with(~ paste0("exp", .), .cols = starts_with("PC"))
  

# PEER Analysis
# num_factors <- 10
# model <- PEER()
# PEER_setPhenoMean(model, as.matrix(t(normalized_counts)))
# PEER_setNk(model, num_factors)
# PEER_update(model)
# 
# peer_factors <- as.data.frame(PEER_getX(model)[, 1:num_factors])
# colnames(peer_factors) <- paste0("PEER", 1:num_factors)
# peer_factors$Sample <- rownames(peer_factors)

# Merge with covariates
cov_full_tbl <- cov_tbl |> 
  inner_join(pc_scores, by = "Sample") |> 
  #inner_join(peer_factors, by = "Sample") |>
  mutate(Sample = factor(Sample, levels = overlap_lst$Sample)) |>
  mutate(Sex = ifelse(Sex == 'Male', 1, 2)) |>
  arrange(Sample) |>
  rename('id' = 'Sample') |>
  write_tsv(paste0(data_dir, "fastqtl_covariates.txt"))

transposed_cov <- cov_full_tbl |>
  as.data.frame() |>
  t() 
rownames(transposed_cov)
write.table(transposed_cov, paste0(data_dir, "fastqtl_covariates_t.txt"), quote = F, sep = '\t', col.names = F)

cov_full_tbl %>%
  column_to_rownames("Sample") %>%  # Move 'Sample' to row names
  as_tibble() %>%                   # Convert to tibble (preserves row names)
  t() %>%                           # Transpose the data
  as_tibble(rownames = "Sample") %>% # Convert back to tibble with row names as 'Sample'
  mutate(Sample = c("", Sample[-1]))

# Load gene lookup tables
hg38_lookup <- get_biomart_gene_lookup('hg38')
hg19_lookup <- get_biomart_gene_lookup('hg19')

# Add hg38 Ensembl IDs using Parse gene lookup
counts_tbl <- t(normalized_counts) |>
  as_tibble(rownames = 'genes') |>
  left_join(gene_lookup, by = join_by(genes == gene_name)) |>
  relocate(genes, gene_id, genome) 
  
# 363 gene not in Parse lookup table
# How were these annotated?? By Scanpy??
# Rm for now
counts_tbl_nas <- counts_tbl |>
  filter(if_any(everything(), is.na)) 

counts_tbl_symbol_dups <- counts_tbl |>
  drop_na() |>
  left_join(hg38_lookup, by = join_by(gene_id == ensembl_gene_id)) |>
  group_by(genes) %>%
  filter(n() > 1) %>%
  ungroup()

counts_tbl_ens_dups <- counts_tbl |>
  drop_na() |>
  left_join(hg38_lookup, by = join_by(gene_id == ensembl_gene_id)) |>
  group_by(gene_id) %>%
  filter(n() > 1) %>%
  ungroup()

counts_tbl <- counts_tbl |>
  drop_na() |>
  left_join(hg38_lookup, by = join_by(gene_id == ensembl_gene_id)) |>
  mutate(
    TSS = if_else(strand == 1, start_position, end_position),  # Determine TSS based on strand
    cis_start = TSS,  # Prevent negative coordinates
    cis_end = TSS + 1 
  ) |>
  distinct(gene_id, .keep_all = TRUE) |> # Keep first occurance dirty for now
  dplyr::select(Chr = chromosome_name, start = cis_start, end = cis_end, TargetID = gene_id, all_of(overlap_lst[[1]])) |>
  arrange(Chr, as.numeric(start), as.numeric(end)) |>
  filter(Chr %in% seq(1,22,1)) |>
  rename('#Chr' = Chr) |> # Required or tabix chokes at fastQTL step
  write_tsv(paste0(data_dir, 'fastqtl_cell_type_1_tmm.bed'))

# counts_tbl |>
#   filter(Chr == 1) |>
#   write_tsv(paste0(data_dir, 'fastqtl_cell_type_1_tmm_22.bed'))
  
# FastQTL choked on test data for any chr not chr1-22
counts_tbl |>
  group_by(`#Chr`) |>
  count() |>
  print(n = Inf)


read_tsv(paste0(data_dir, 'fastqtl_cell_type_1_tmm.bed')) |>
  anyNA()
#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
  

