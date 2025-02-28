#--------------------------------------------------------------------------------------
#
#    Prep QTL input files - draft
#
#--------------------------------------------------------------------------------------

## Info  ------------------------------------------------------------------------------

#  Prep input bed and covarite files for FastQTL 
#  Running locally for now

# TODO:

#- Add a loop for all genes
#- Do we need TMM and CPM?
#- Add proper covariates when full data available
#- Try and salvage NA genes across symbol and ens IDs for each cell type
#- Need a comprehensive strategy to deal with duplicates (dirty solution for now)
#- Explore Peer and different PC num for covariates (variance partition)
#- Add sublibrary, plate 
#- Some cell types do not fully overlap with overlap list changed final counts tbl 
#> to any_of instead of all_of to avoid error

##  Load Packages, functions and variables  -------------------------------------------
# Install and load required libraries
library(edgeR)
library(variancePartition)
library(tidyverse)
#library(peer) # Old package: non-standard installation

## Load data  -------------------------------------------------------------------------
data_dir <- "~/Desktop/eQTL_study_2025/results/03SCANPY/pseudobulk/"
cov_dir <- "~/Desktop/eQTL_study_2025/resources/test_genotypes/"
lookup_dir <- "~/Desktop/eQTL_study_2025/results/02PARSE/combine_plate1/all-sample/DGE_filtered/"
gene_lookup <- read_csv(paste0(lookup_dir, 'all_genes.csv'), col_types = 'ccc')
cell_types <- c('RG', 'ExN-1', 'InN-1', 'ExN-2', 'InN-2', 'Endo-Peri', 'MG', 'Mig-N')

# Load biomart gene lookup tables
hg38_lookup <- get_biomart_gene_lookup('hg38')
hg19_lookup <- get_biomart_gene_lookup('hg19')

# Need to order samples in input files in same order
overlap_lst <- read_tsv(paste0(cov_dir, 'overlapping_samples.txt'), col_types = 'c', col_names = 'Sample') 

## Covariates
cov_orig_tbl <- read_tsv(paste0(cov_dir, 'filtered_covariates.txt'), col_types = 'ccd') 
cov_tbl <- cov_orig_tbl |> 
  dplyr::select(Sample, Sex, PCW, PC1, PC2, PC3) |>
  rename_with(~ paste0("gen", .), .cols = matches("^PC\\d+$"))

## Generate cell-specific TMM normalised counts and covariates for QTL analyses  ------

#for (cell_type in cell_types) {

report_list <- list()

for (cell_type in cell_types) {
  
  message('Creating TMM normalised counts and covariate files for: ', cell_type)
  report_tibble <- tibble(cell_type = cell_type)
  
  # Counts: don;t use dplyr as we want to keep rownames
  message('Loading counts ... \n')
  pseudobulk_counts <- read.csv(paste0(data_dir, cell_type, "_pseudobulk.csv"), row.names = 1) # rows=genes, columns=samples
  
  message('Check for samples not in overlap list: ', 
          paste(setdiff(overlap_lst[[1]], rownames(pseudobulk_counts)), collapse = ' '))
  
  if (sum(colSums(pseudobulk_counts) == 0) != 0) {
    message('Removed ', sum(colSums(pseudobulk_counts) == 0), ' genes not expressed in cell type')
    pseudobulk_counts <- pseudobulk_counts <- pseudobulk_counts[, colSums(pseudobulk_counts != 0) > 0]}
  print(pseudobulk_counts[1:5, 1:5])

  # TMM normalization
  # Create a DGEList object
  message('\nTMM normalising counts ... \n')
  dge <- DGEList(counts = pseudobulk_counts)
  dge <- calcNormFactors(dge, method = "TMM") # Do we need two normalisation factors here?
  normalized_counts <- cpm(dge, normalized.lib.sizes = TRUE)
  print(normalized_counts[1:5, 1:5])
  
  # PCA Analysis
  message('\nRunning PCA on expression values ... ')
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
  message('Creating covariate files ... ')
  cov_full_tbl <- cov_tbl |> 
    inner_join(pc_scores, by = "Sample") |> 
    #inner_join(peer_factors, by = "Sample") |>
    mutate(Sample = factor(Sample, levels = overlap_lst$Sample)) |>
    mutate(Sex = ifelse(Sex == 'Male', 1, 2)) |>
    arrange(Sample) |>
    rename('id' = 'Sample') |>
    write_tsv(paste0(data_dir, cell_type, "_covariates.txt"))
  
  transposed_cov <- cov_full_tbl |>
    as.data.frame() |>
    t() 
  rownames(transposed_cov)
  write.table(transposed_cov, paste0(data_dir, cell_type, "_covariates_t.txt"), quote = F, sep = '\t', col.names = F)
  
  # cov_full_tbl %>%
  #   column_to_rownames("Sample") %>%  # Move 'Sample' to row names
  #   as_tibble() %>%                   # Convert to tibble (preserves row names)
  #   t() %>%                           # Transpose the data
  #   as_tibble(rownames = "Sample") %>% # Convert back to tibble with row names as 'Sample'
  #   mutate(Sample = c("", Sample[-1]))
  
  # Add hg38 Ensembl IDs using Parse gene lookup
  counts_tbl <- t(normalized_counts) |>
    as_tibble(rownames = 'genes') |>
    left_join(gene_lookup, by = join_by(genes == gene_name)) |>
    relocate(genes, gene_id, genome) 
  message('Counts tbl dimensions (Genes x [Samples + 3 annotation cols]): ', 
          paste0(dim(counts_tbl)[1], ' x ', dim(counts_tbl)[2]))
  report_tibble$initial_dims <- paste0(dim(pseudobulk_counts)[1], "x", dim(pseudobulk_counts)[2])
    
  # Some genes are not in Parse lookup table why?
  # How were these annotated?? By Scanpy??
  # Rm for now
  counts_tbl_nas <- counts_tbl |>
    filter(if_any(everything(), is.na)) 
  message('Counts tbl gene ID NAs: ', dim(counts_tbl_nas)[1])
  report_tibble$NAs <- dim(counts_tbl_nas)[1]
  
  counts_tbl_symbol_dups <- counts_tbl |>
    drop_na() |>
    left_join(hg38_lookup, by = join_by(gene_id == ensembl_gene_id)) |>
    group_by(genes) %>%
    filter(n() > 1) %>%
    ungroup()
  message('Counts tbl gene symbol dups: ', length(counts_tbl_symbol_dups$genes |> unique()))
  report_tibble$symbol_dups <- length(counts_tbl_symbol_dups$genes |> unique())
  
  counts_tbl_ens_dups <- counts_tbl |>
    drop_na() |>
    left_join(hg38_lookup, by = join_by(gene_id == ensembl_gene_id)) |>
    group_by(gene_id) %>%
    filter(n() > 1) %>%
    ungroup()
  message('Counts tbl gene ensembl dups: ', length(counts_tbl_ens_dups$genes |> unique()))
  report_tibble$ensembl_dups <- length(counts_tbl_ens_dups$genes |> unique())
  
  counts_tbl <- counts_tbl |>
    drop_na() |>
    left_join(hg38_lookup, by = join_by(gene_id == ensembl_gene_id)) |>
    mutate(
      TSS = if_else(strand == 1, start_position, end_position),  # Determine TSS based on strand
      cis_start = TSS,  # Prevent negative coordinates
      cis_end = TSS + 1 
    ) |>
    distinct(gene_id, .keep_all = TRUE) |> # Keep first occurrence dirty for now
    dplyr::select(Chr = chromosome_name, start = cis_start, end = cis_end, TargetID = gene_id, any_of(overlap_lst[[1]])) |>
    arrange(Chr, as.numeric(start), as.numeric(end)) |>
    filter(Chr %in% seq(1,22,1)) |>
    rename('#Chr' = Chr) |> # Required or tabix chokes at fastQTL step
    write_tsv(paste0(data_dir, cell_type, '_tmm.bed'))
  
  message('Final counts tbl dims: ', paste0(dim(counts_tbl)[1], ' x ', dim(counts_tbl)[2]), '\n')
  report_tibble$final_dims <- paste0(dim(counts_tbl)[1], 'x', dim(counts_tbl)[2])
  
  report_list[[cell_type]] <- report_tibble
  
}

# Combine all reports into a single tibble
final_report <- bind_rows(report_list) |>
  print()

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
  

