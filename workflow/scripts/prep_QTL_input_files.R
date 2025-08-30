#--------------------------------------------------------------------------------------
#
#    Prep input files for tensorQTL
#
#--------------------------------------------------------------------------------------

# As output wildcards need to be standardised in smk, I run this twice each instance
# tracking either the gene exp {cell_type}_tmm.bed files or the covaraite
# {cell_type}_genPC_{geno_pc}_expPC_{exp_pc}_tmm.bed files. Could separate the scripts
# but I'd like to run the checks together

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

# Input and output paths
# cov_file <- snakemake@input[["cov_file"]]
# sex_file <- snakemake@input[["sex_file"]]
# gene_lookup <- snakemake@input[["gene_lookup"]]
# cov_out <- snakemake@input[["cov_out"]]
# exp_out <- snakemake@input[["exp_out"]]
# pseudoblk_dir <- snakemake@input[["pseudoblk_dir"]]
# report_dir <- snakemake@input[["report_dir"]]
# out_dir <- snakemake@input[["out_dir"]]

cov_file <- "~/Desktop/eQTL_study_2025/results/04GENOTYPES/TOPMED/covariates/pca.eigenvec"
sex_file <- "~/Desktop/eQTL_study_2025/results/03SCANPY/sex_assign_final_tbl.tsv"
cell_types <- c('RG', 'InN', 'ExN-2', 'OPC', 'Endo-Peri', 'MG')
gene_lookup <- "~/Desktop/eQTL_study_2025/resources/sheets/gene_lookup_hg38.tsv"
geno_pcs <- 3
exp_pcs <- c(10, 20, 30, 40, 50)
pseudobulk_dir <- "~/Desktop/eQTL_study_2025/results/03SCANPY/pseudobulk/"
out_dir <- "~/Desktop/eQTL_study_2025/results/05TENSORQTL/"
report_dir <- "~/Desktop/eQTL_study_2025/workflow/reports/05TENSORQTL/"

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
message("Sample Counts\n")
message("============================\n")

message("Sample cnt in genotype covariate file: ", nrow(genotype_cov_tbl))
message("Sample cnt in metadata file: ", nrow(meta_data))

message("Counts by PCW:\n")
meta_data %>%
  group_by(PCW) %>%
  dplyr::count() %>%
  knitr::kable(format = "simple", align = "c") %>%
  print()

message("\nCounts by Sex:\n")
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
for (cell_type in cell_types) {
  
  message('\n\nCreating TMM normalised counts and covariate files for: ', cell_type)
  report_tibble <- tibble(cell_type = cell_type)
  
  # Load counts: don't use dplyr as we want to keep rownames
  message('\nLoading counts ... \n')
  pseudblk_cnts <- read.csv(paste0(pseudblk_dir, cell_type, "_pseudobulk.csv"), row.names = 1)
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

  # TMM normalization
  # Create a DGEList object
  message('\nTMM normalising counts ... \n')
  dge <- DGEList(counts = pseudblk_cnts)
  dge <- calcNormFactors(dge, method = "TMM") # Do we need two normalisation factors here?
  normalised_cnts <- edgeR::cpm(dge, normalized.lib.sizes = TRUE)
  print(normalised_cnts[1:5, 1:5])
  message('\nDimensions of normalised counts: ', 
          paste0(dim(normalised_cnts)[1], ' x ', dim(normalised_cnts)[2]))
  
  
  # Filter genes with log2 CPM < 2.0 in all samples - Fugita
  # message('\nFiltering genes with log2 CPM < 2.0 in all samples ... \n')
  # log2_cpm <- log2(normalised_cnts + 1)  # Add 1 to avoid log(0)
  # print(log2_cpm[1:5, 1:5])
  # low_expr_genes <- rowSums(log2_cpm >= 2.0) == 0  # Genes with log2 CPM < 2.0 in all samples
  # if (sum(low_expr_genes) > 0) {
  #   message('Removed ', sum(low_expr_genes), ' genes with log2 CPM < 2.0 in all samples for ', cell_type)
  #   normalised_cnts <- normalised_cnts[!low_expr_genes, , drop = FALSE]
  # }
  # message('\nDimensions after log2 CPM filtering: ', 
  #         paste0(dim(normalised_cnts)[1], ' x ', dim(normalised_cnts)[2]))
  
  #Filter genes with mean CPM < 1 - Bryois - more stringent
  message('\nFiltering genes with mean CPM < 1 ... \n')
  mean_cpm <- rowMeans(normalised_cnts)
  low_expr_genes <- mean_cpm < 1
  if (sum(low_expr_genes) > 0) {
    message('Removed ', sum(low_expr_genes), ' genes with mean CPM < 1 for ', cell_type)
    normalised_cnts <- normalised_cnts[!low_expr_genes, , drop = FALSE]
  }
  message('\nDimensions after mean CPM filtering: ',
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
  message('Writing ', cell_type, ' covariate matrix ...')
  write.table(cov_matrix, 
              file = paste0(out_dir, cell_type, "_covariates.txt"), 
              sep = "\t", 
              quote = FALSE, 
              col.names = NA,    # Write sample IDs as column names
              row.names = TRUE,    # Write covariate names as row names
              append = FALSE)


  ## Prepare gene expression data -----
  
  # Gene expression data
  # Add Ensembl IDs using BiomaRt Gene lookup; some genes are symbol, so are emsembl ID
  message('Harmonising gene expression gene IDs to Ensembl ... ')
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
  message('Writing final tmm normalised gene expression matrix ... ')
  pseudblk_cnts <- pseudblk_cnts |>
    drop_na() |>
    mutate(
      TSS = if_else(strand == 1, start_position, end_position),  # Determine TSS based on strand
      cis_start = TSS,  # Prevent negative coordinates
      cis_end = TSS + 1 
    ) |>
    distinct(ensembl_gene_id, .keep_all = TRUE) |> # Keep first occurrence dirty for now
    dplyr::select(Chr = chromosome_name, start = cis_start, end = cis_end, TargetID = ensembl_gene_id, any_of(available_samples)) |>
    rbind(pseudblk_ensembl_cnts) |> # Salvage genes marked as ENSEMBL ID, rather than symbol initally
    distinct(TargetID, .keep_all = TRUE) |>
    arrange(Chr, as.numeric(start), as.numeric(end)) |>
    dplyr::filter(Chr %in% seq(1,22,1)) |>
    dplyr::rename('#Chr' = Chr) |> # Required or tabix chokes at fastQTL step
    write_tsv(paste0(out_dir, cell_type, '_tmm.bed'))
  
  message('Final counts tbl dims: ', paste0(dim(pseudblk_cnts)[1], ' x ', dim(pseudblk_cnts)[2]), '\n')
  report_tibble$final_dims <- paste0(dim(pseudblk_cnts)[1], 'x', dim(pseudblk_cnts)[2])
  report_list[[cell_type]] <- report_tibble
  
}

# Combine all reports into a single tibble
final_report <- bind_rows(report_list) |>
  print()

write_rds(final_report, paste0(report_dir, 'gene_exp_matrix_report.rds'))


#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
  


