#--------------------------------------------------------------------------------------
#
#    Prep QTL input files - draft
#
#--------------------------------------------------------------------------------------

## Info  ------------------------------------------------------------------------------

#  Prep input bed and covarite files for FastQTL 
#  Running locally for now

# TODO:

#- Try and salvage NA genes across symbol and ens IDs for each cell type
#- Need a comprehensive strategy to deal with duplicates (dirty solution for now)
#- Explore Peer and different PC num for covariates (variance partition)
#- Add sublibrary, plate?
#- Some cell types do not fully overlap with overlap list changed final counts tbl 
#- to any_of instead of all_of to avoid error
#- Duplicate sample 14993 and 14993_FC, keeping the latter (need to remove this before imputation)

##  Load Packages, functions and variables  -------------------------------------------
# Install and load required libraries
library(edgeR)
library(variancePartition)
library(tidyverse)
library(readxl)
library(ggcorrplot)
#library(peer) # Old package: non-standard installation


## Load data  -------------------------------------------------------------------------
data_dir <- "~/Desktop/eQTL_study_2025/results/03SCANPY/pseudobulk/"
cov_dir <- "~/Desktop/eQTL_study_2025/results/04GENOTYPES/TOPMED/covariates/"
resources_dir <- "~/Desktop/eQTL_study_2025/resources/sheets/"
scanpy_dir <- "~/Desktop/eQTL_study_2025/results/03SCANPY/"
script_dir <- "~/Desktop/eQTL_study_2025/workflow/scripts/"
list_dir <- "~/Desktop/eQTL_study_2025/results/04GENOTYPES/TOPMED/filtered/"
lookup_dir <- "~/Desktop/eQTL_study_2025/results/02PARSE/combine_plate1/all-sample/DGE_filtered/"
geno_dir <- "~/Desktop/GenotypeQCtoHRC/eqtlhg38tm/"
gene_lookup <- read_csv(paste0(lookup_dir, 'all_genes.csv'), col_types = 'ccc')
cell_types <- c('RG', 'ExN-1', 'InN-1', 'ExN-2', 'InN-2', 'Endo-Peri', 'MG', 'Mig-N')
source(paste0(script_dir, 'functions.R'))

# Load biomart gene lookup tables
hg38_lookup <- get_biomart_gene_lookup('hg38')
hg19_lookup <- get_biomart_gene_lookup('hg19')

# Need to order samples in input files in same order
sample_lst <- read_tsv(paste0(list_dir, 'geno_sample_list.txt'), col_types = 'c', col_names = 'sample') 

## Covariates
genotype_cov_tbl <- read_delim(paste0(cov_dir, 'pca.eigenvec'), delim = ' ', 
                      col_names = c('sample', 'sample2', paste0('genPC', seq(1,10,1)))) |>
  dplyr::select(-sample2) 
  
## Metadata
cov_tbl <- read_excel(paste0(resources_dir, 'Fetal single cell eQTL_genotypes_hg19.xlsx'), sheet = 'Sheet2') |>
  dplyr::rename(sample = ...1) |>
  mutate(sample = str_replace_all(sample, " \\(A\\)", "_A")) |>
  mutate(sample = str_replace_all(sample, " \\(a\\)", "_a")) |>
  dplyr::select(sample, PCW, Sex) |>
  left_join(genotype_cov_tbl) |>
  print(n = Inf)

cov_tbl |>
  group_by(PCW) |>
  dplyr::count()
  
cov_tbl |>
  group_by(Sex) |>
  dplyr::count()

message('Samples with NA across genotypes: ')
cov_tbl |>
  filter(if_any(starts_with("genPC"), ~ is.na(.)))

cov_tbl <- cov_tbl |>
  filter(if_any(starts_with("genPC"), ~ !is.na(.)))

## Sex assignments from expression data
sex_assign_exp_tbl <- read_csv(paste0(scanpy_dir, 'sex_assignments.csv')) |>
  dplyr::rename(sample = Sample_ID) |>
  mutate(Sex_exp = case_when(
    Predicted_Sex == 'Female' ~ "F",
    Predicted_Sex == 'Male' ~ "M",
    Predicted_Sex == 'Ambiguous' ~ "A"))

sex_assign_geno_tbl <- read.table(
  paste0(geno_dir, 'eQTL_genotypes_hg19.qc3.sexcheck'), header = TRUE) |>
  as_tibble() |>
  mutate(sex_geno = case_when(
    F <= 0.2 ~ 'F',
    F >= 0.8 ~ 'M'
  ))

# Compare
sex_assign_final_tbl <- cov_tbl |>
  dplyr::select(sample, meta_sex = Sex) |>
  left_join(sex_assign_exp_tbl |> dplyr::select(sample, Sex_exp)) |>
  left_join(sex_assign_geno_tbl |> dplyr::select(sample = FID, sex_geno)) |>
  mutate(
    # Count "F" and "M" occurrences across the three columns
    f_count = rowSums(across(c(meta_sex, Sex_exp, sex_geno), ~ . == "F"), na.rm = TRUE),
    m_count = rowSums(across(c(meta_sex, Sex_exp, sex_geno), ~ . == "M"), na.rm = TRUE),
    # Determine sex_final based on counts
    sex_final = case_when(
      f_count == 3 | (f_count == 2 & m_count <= 1) ~ "F",  # All F or 2F + 0/1M/NA
      m_count == 3 | (m_count == 2 & f_count <= 1) ~ "M",  # All M or 2M + 0/1F/NA
      TRUE ~ "A"                                            # Otherwise ambiguous
    ),
    # Assign numeric codes: F=2, M=1, A=0
    sex_code = case_when(
      sex_final == "F" ~ 2,
      sex_final == "M" ~ 1,
      sex_final == "A" ~ 0
    )
  ) 

cov_tbl <- cov_tbl |>
  left_join(sex_assign_final_tbl |> dplyr::select(sample, sex_code)) |>
  dplyr::select(-Sex) |>
  relocate(sample, PCW, Sex = sex_code)

## Generate cell-specific TMM normalised counts and covariates for QTL analyses  ------

#for (cell_type in cell_types) {

report_list <- list()

for (cell_type in cell_types) {
  
  message('Creating TMM normalised counts and covariate files for: ', cell_type)
  report_tibble <- tibble(cell_type = cell_type)
  
  # Counts: don't use dplyr as we want to keep rownames
  message('Loading counts ... \n')
  pseudobulk_counts <- read.csv(paste0(data_dir, cell_type, "_pseudobulk.csv"), row.names = 1) # rows=genes, columns=samples
  
  message('Check for samples not in overlap list: \n', 
          paste(setdiff(sample_lst[[1]], rownames(pseudobulk_counts)), collapse = ' '))
  
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
  exp_pc_scores <- as.data.frame(pca$x[, 1:30]) |> 
    as_tibble(rownames = 'sample') |>
    rename_with(~ paste0("exp", .), .cols = starts_with("PC"))
  
  # Variance explained by expression PCs
  exp_var_explained <- summary(pca)$importance[2, 1:30]  # Proportion of variance
  report_tibble$exp_var_explained <- list(exp_var_explained)
    
  # Combine with covariates
  cov_full_tbl <- cov_tbl |> 
    left_join(exp_pc_scores, by = "sample") |>
    mutate(PCW = str_replace(PCW, '1st tri\\s*\\?', '0')) |>
    mutate(PCW = replace(PCW, is.na(PCW), '0')) |>
    mutate(PCW = as.numeric(PCW))
  
  # Correlation between genotype and expression PCs
  cor_matrix <- cor(cov_full_tbl[-1], use = "pairwise.complete.obs")
  
  # Store results
  report_tibble$cor_matrix <- list(cor_matrix)
  report_list[[cell_type]] <- report_tibble

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
    inner_join(exp_pc_scores, by = "sample") |> 
    #inner_join(peer_factors, by = "Sample") |>
    mutate(sample = factor(sample, levels = sample_lst$sample)) |>
    rename('id' = 'sample') |>
    dplyr::select(id, PCW, Sex, genPC1:genPC3, expPC1:expPC10) |>
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
  # Potentially blanks and duplicates here as more rows after join!!
  pseudobulk_counts <- t(normalized_counts) |>
    as_tibble(rownames = 'genes') |>
    left_join(gene_lookup, by = join_by(genes == gene_name)) |>
    relocate(genes, gene_id, genome) 
  message('Counts tbl dimensions (Genes x [Samples + 3 annotation cols]): ', 
          paste0(dim(pseudobulk_counts)[1], ' x ', dim(pseudobulk_counts)[2]))
  report_tibble$initial_dims <- paste0(dim(pseudobulk_counts)[1], "x", dim(pseudobulk_counts)[2])
    
  # Some genes are not in Parse lookup table why?
  # How were these annotated?? By Scanpy??
  # Rm for now
  pseudobulk_counts_nas <- pseudobulk_counts |>
    filter(if_any(everything(), is.na)) 
  message('Counts tbl gene ID NAs: ', dim(pseudobulk_counts_nas)[1])
  report_tibble$NAs <- dim(pseudobulk_counts_nas)[1]
  
  pseudobulk_counts_symbol_dups <- pseudobulk_counts |>
    drop_na() |>
    left_join(hg38_lookup, by = join_by(gene_id == ensembl_gene_id)) |>
    group_by(genes) %>%
    filter(n() > 1) %>%
    ungroup()
  message('Counts tbl gene symbol dups: ', length(pseudobulk_counts_symbol_dups$genes |> unique()))
  report_tibble$symbol_dups <- length(pseudobulk_counts_symbol_dups$genes |> unique())
  
  pseudobulk_counts_ens_dups <- pseudobulk_counts |>
    drop_na() |>
    left_join(hg38_lookup, by = join_by(gene_id == ensembl_gene_id)) |>
    group_by(gene_id) %>%
    filter(n() > 1) %>%
    ungroup()
  message('Counts tbl gene ensembl dups: ', length(pseudobulk_counts_ens_dups$genes |> unique()))
  report_tibble$ensembl_dups <- length(pseudobulk_counts_ens_dups$genes |> unique())
  
  pseudobulk_counts <- pseudobulk_counts |>
    drop_na() |>
    left_join(hg38_lookup, by = join_by(gene_id == ensembl_gene_id)) |>
    mutate(
      TSS = if_else(strand == 1, start_position, end_position),  # Determine TSS based on strand
      cis_start = TSS,  # Prevent negative coordinates
      cis_end = TSS + 1 
    ) |>
    distinct(gene_id, .keep_all = TRUE) |> # Keep first occurrence dirty for now
    dplyr::select(Chr = chromosome_name, start = cis_start, end = cis_end, TargetID = gene_id, any_of(sample_lst[[1]])) |>
    arrange(Chr, as.numeric(start), as.numeric(end)) |>
    filter(Chr %in% seq(1,22,1)) |>
    rename('#Chr' = Chr) |># Required or tabix chokes at fastQTL step
    write_tsv(paste0(data_dir, cell_type, '_tmm.bed'))
  
  message('Final counts tbl dims: ', paste0(dim(pseudobulk_counts)[1], ' x ', dim(pseudobulk_counts)[2]), '\n')
  report_tibble$final_dims <- paste0(dim(pseudobulk_counts)[1], 'x', dim(pseudobulk_counts)[2])
  
  report_list[[cell_type]] <- report_tibble
  
}

# Combine all reports into a single tibble
final_report <- bind_rows(report_list) |>
  print()

# Plot
var_explained_plt_lst <- list()
cor_plt_lst <- list()

for (cell_type in names(report_list)) {
  
  exp_var_explained <- report_list[[cell_type]]$exp_var_explained[[1]]
  
  # Variance Explained Plot for Expression PCs
  var_df <- tibble(
    PC = 1:30,
    Variance = exp_var_explained,
    Cumulative = cumsum(exp_var_explained)
  )
  p_var <- ggplot(var_df, aes(x = PC)) +
    geom_line(aes(y = Variance), color = "blue") +
    geom_point(aes(y = Variance), color = "blue") +
    geom_line(aes(y = Cumulative), color = "red") +
    labs(title = paste("Expression PC Variance Explained -", cell_type),
         y = "Proportion of Variance",
         x = "Principal Component") +
    theme_minimal()
  print(p_var)
  var_explained_plt_lst[[cell_type]] <- p_var
  #ggsave(paste0(data_dir, cell_type, "_exp_pc_variance.png"), p_var, width = 8, height = 6)
  
  # Correlation Heatmap
  p_cor <- ggcorrplot::ggcorrplot(cor_matrix, 
                      hc.order = F, 
                      type = "lower", 
                      lab = TRUE, 
                      lab_size = 2,
                      title = paste("Correlation of PCs -", cell_type)) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  print(p_cor)
  cor_plt_lst[[cell_type]] <- p_cor
# ggsave(paste0(data_dir, cell_type, "_pc_correlation.png"), p_cor, width = 10, height = 10)

}

cowplot::plot_grid(plotlist = var_explained_plt_lst)
cowplot::plot_grid(plotlist = cor_plt_lst)
#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
  


