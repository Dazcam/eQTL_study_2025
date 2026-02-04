#--------------------------------------------------------------------------------------
#
#    Create corr tbl for my sig. SMR eQTL and Fugita nominal eQTL (single cell pairs)
#
#--------------------------------------------------------------------------------------


## Info  ------------------------------------------------------------------------------

# Set up logging for Snakemake
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
p_smr <- 0.05
p_heidi <- 0.01
disorder <- 'scz'

celltype_map <- tibble(
  my_cell     = c("Glu-UL", "Glu-DL", "GABA"),
  fugita_cell = c("Exc",    "Exc",    "Inh")
)

paired_results <- list()

for (i in seq_len(nrow(celltype_map))) {
  
  my_ct <- celltype_map$my_cell[i]
  fu_ct <- celltype_map$fugita_cell[i]
  
  message('\n\nGenerating beta cor table for ', my_ct, ' and ', fu_ct,'\n')
  
  message('Loading sig.eQTL data for', my_ct,'\n')
  file_path <- file.path(in_dir, my_ct, paste0(my_ct, "_", disorder, ".smr"))
  
  # Read the file
  smr <- suppressMessages(read_tsv(file_path))
  n_probes <- length(unique(smr$probeID))
  smr_sig <- smr |> filter(p_SMR < (p_smr / n_probes))
  smr_sig_heidi <- smr_sig |> filter(p_HEIDI > p_heidi)
  
  message('n_probes: ', n_probes)
  message('Sig SNPs at p_SMR: ', smr_sig)
  message('Sig SNPs at p_SMR and non-sig HEIDI: ', smr_sig_heidi)
  
  df_sig <- smr_sig_heidi %>%
    select(snp = topSNP, gene = Gene, beta = b_eQTL) %>%
    filter(str_detect(snp, "^rs"), str_detect(gene, "^ENSG"))
  message('\nSig. eQTL in ', my_ct, ': ', nrow(df_sig), '\n')
  
  # Pool: for dup SNP-gene, keep largest abs(beta)
  pooled_my <- df_sig %>%
    # group_by(snp, gene) %>%
    # slice_max(abs(beta), n = 1, with_ties = FALSE) %>%  # If ties in abs, picks first
    # ungroup() |>
    # filter(str_detect(snp, '^rs')) |> 
    # filter(str_detect(gene, '^ENSG')) |>
    mutate(key = paste(snp, gene, sep = '_')) |>
    select(key, beta_my = beta)
  message('Pooled sig. eQTL after dup rm : ', nrow(pooled_my))
  message('Any NAs in pooled sig. eQTL? ', anyNA(pooled_my))
  message('\nPooled sig. eQTL tbl:\n')
  print(pooled_my)

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
