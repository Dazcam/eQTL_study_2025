#--------------------------------------------------------------------------------------
#
#    Create corr tbl for my sig. eQTL and Fugita nominal eQTL (single cell pairs)
#
#--------------------------------------------------------------------------------------


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

message("\n\nGenerating data for eQTL beta correlation analysis (single cell pairs) ...")

# -------------------------------------------------------------------------------------

out_file <- snakemake@output[[1]]
out_dir <- dirname(out_file)
in_dir <- snakemake@params[['in_dir']] 
fugita_dir <- snakemake@params[['fugita_dir']] 
#'../resources/public_datasets/fugita_2024/'
#'../results/05TENSORQTL/tensorqtl_perm/'

# Load required libraries
library(tidyverse)
library(data.table)

# Define variables
gen_PCs <- 4
norm_method <- 'quantile' 

expPC_map_all <- c(
  "Glu-UL" = 50,
  "Glu-DL" = 40,
  "GABA"   = 30
)

celltype_map <- tibble(
  my_cell     = c("Glu-UL", "Glu-DL", "GABA"),
  fugita_cell = c("Exc",    "Exc",    "Inh")
)


paired_results <- list()

for (i in seq_len(nrow(celltype_map))) {
  
  my_ct <- celltype_map$my_cell[i]
  fu_ct <- celltype_map$fugita_cell[i]
  expPC <- expPC_map_all[[my_ct]]
  
  message('\n\nGenerating beta cor table for ', my_ct, ' and ', fu_ct,'\n')
  
  message('Loading sig.eQTL data for', my_ct,'\n')
  log_dir <- paste0(in_dir, my_ct, '_', norm_method, '_genPC_', gen_PCs, '_expPC_', expPC, '/')
  file_path <- paste0(log_dir, my_ct, '_', norm_method, '_perm.cis_qtl.txt.gz')
  
  # Read the file
  df_sig <- read_delim(file_path, delim = '\t', col_names = TRUE) %>%
    filter(qval < 0.05) %>%
    select(snp = variant_id, gene = phenotype_id, beta = slope) %>%
    filter(str_detect(snp, "^rs"), str_detect(gene, "^ENSG"))
  message('\nSig. eQTL in ', my_ct, ': ', nrow(df_sig), '\n')
  
  write_tsv(df_sig, file.path(out_dir, paste0(my_ct, "_beta_cor_single_df_sig.tsv")))
  
  # Pool: for dup SNP-gene, keep largest abs(beta)
  pooled_my <- df_sig %>%
    group_by(snp, gene) %>%
    slice_max(abs(beta), n = 1, with_ties = FALSE) %>%  # If ties in abs, picks first
    ungroup() |>
    # filter(str_detect(snp, '^rs')) |> 
    # filter(str_detect(gene, '^ENSG')) |>
    mutate(key = paste(snp, gene, sep = '_')) |>
    select(key, beta_my = beta)
  message('Pooled sig. eQTL after dup rm : ', nrow(pooled_my))
  message('Any NAs in pooled sig. eQTL? ', anyNA(pooled_my))
  message('\nPooled sig. eQTL tbl:\n')
  print(pooled_my)
  
  write_tsv(pooled_my, file.path(out_dir, paste0(my_ct, "_beta_cor_single_pooled_my.tsv")))

  
  message('\nLoading Fugita data ...\n')
  fugita_path <- paste0(fugita_dir, 'celltype-eqtl-sumstats.', fu_ct,'.tsv')
  
  # Read the file
  fugita_df <- read_delim(fugita_path, delim = '\t', col_names = TRUE, na = c("", "NA", "na", "-")) |>
    select(snp = snps, gene = gene_id, beta, celltype) |>
    mutate(beta = as.numeric(beta)) |>
    filter(str_detect(snp, '^rs'),
           str_detect(gene, '^ENSG'),
           !is.na(beta))
  message('\nTotal SNP-gene pairs in Fugita ', fu_ct, ': ' , nrow(fugita_df), '\n')

  setDT(fugita_df)
  
  pooled_fugita <- fugita_df[
    , .SD[which.max(abs(beta))],
    by = .(snp, gene)
  ][
    , .(
      key = paste(snp, gene, sep = "_"),
      beta_fugita = beta
    )
  ] %>%
    as_tibble()
  
  rm(fugita_df)
  gc()

  message('\nFugita SNP-gene pairs after dup rm: ', nrow(pooled_fugita))
  message('Any NAs in pooled Fugita data? ', anyNA(pooled_fugita))
  message('\nPooled Fugita tbl:\n')
  print(pooled_fugita)

  # Join with your pooled data for paired betas
  paired_betas <- pooled_my %>%
    inner_join(pooled_fugita, by = c('key')) %>%
    select(key, beta_my, beta_fugita)
  paired_results[[my_ct]] <- paired_betas
  message('Number of overlapping SNP-gene pairs: ', nrow(paired_betas))
  message('Any NAs in paired betas tbl? ', anyNA(paired_betas))
  message('Paired betas tbl:\n')
  print(paired_betas)
  
  message('Writing tbl to:', file.path(out_dir, paste0(my_ct, "_beta_cor_single_tbl.tsv")))
  write_tsv(paired_betas, file.path(out_dir, paste0(my_ct, "_beta_cor_single_tbl.tsv")))
                                    
  
}

file.create(out_file)
message('All Done.', out_file)

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
