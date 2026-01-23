#--------------------------------------------------------------------------------------
#
#    Create data for my sig. eQTL and Fugita nominal eQTL Beta correlation
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

message("\n\nGenerating data for eQTL beta correlation analysis ...")

# -------------------------------------------------------------------------------------

out_file <- snakemake@output[[1]]
in_dir <- snakemake@params[['in_dir']] 
fugita_dir <- snakemake@params[['fugita_dir']] 
#'../resources/public_datasets/fugita_2024/'
#'../results/05TENSORQTL/tensorqtl_perm/'

# Load required libraries
library(tidyverse)
library(data.table)

# Define variables
cell_types <- c(
  "Glu-UL", "Glu-UL-0", "Glu-UL-1", "Glu-UL-2",
  "Glu-DL", "Glu-DL-0", "Glu-DL-1", "Glu-DL-2", 
  "NPC", "NPC-0", "NPC-1", "NPC-2", 
  "GABA", "GABA-0", "GABA-1", "GABA-2",
  "Endo-Peri",
  "OPC",
  "MG"
)

# Exp PC map
expPC_map_all <- c(
  "Glu-UL"     = 50,
  "Glu-DL"     = 40,
  "GABA"       = 30,
  "NPC"        = 30,
  "MG"         = 30,
  "OPC"        = 30,
  "Endo-Peri"  = 30,
  "Glu-UL-0"     = 30,
  "Glu-UL-1"     = 50,
  "Glu-UL-2"     = 40,
  "Glu-DL-0"     = 40,
  "Glu-DL-1"     = 40,
  "Glu-DL-2"     = 30,
  "GABA-0"       = 30,
  "GABA-1"       = 40,
  "GABA-2"       = 40,
  "NPC-0"        = 50,
  "NPC-1"        = 30,
  "NPC-2"        = 20
)


exp_PCs <- 40
gen_PCs <- 4
norm_method <- 'quantile' 



fugita_cell_types <- c('Ast', 'Exc', 'Inh', 'Oli', 'OPC', 'Mic', 'End')

# Step 1: Pool unique sig. SNP-gene pairs from eQTL data, max(abs(beta)) for dups
my_data_list <- list()

message('Loading sig.eQTL data ...\n')
for (cell_type in cell_types) {
  
  expPC <- expPC_map_all[[ct]]
  
  log_dir <- paste0(in_dir, cell_type, '_', norm_method, '_genPC_', gen_PCs, '_expPC_', exp_PCs, '/')
  file_path <- paste0(log_dir, cell_type, '_', norm_method, '_perm.cis_qtl.txt.gz')
  
  # Read the file
  df_sig <- read_delim(file_path, delim = '\t', col_names = TRUE) %>%
    filter(qval < 0.05) %>%
    select(snp = variant_id, gene = phenotype_id, beta = slope)
  message('\nSig. eQTL in ', cell_type, ':', nrow(df_sig), '\n')
  
  # Append to list
  my_data_list[[cell_type]] <- df_sig
  rm(df_sig)
}

# Combine all significant data
all_my_sig <- bind_rows(my_data_list)
message('Pooled Sig. eQTL before dup rm: ', nrow(all_my_sig))

# Pool: for dup SNP-gene, keep largest abs(beta)
pooled_my <- all_my_sig %>%
  group_by(snp, gene) %>%
  slice_max(abs(beta), n = 1, with_ties = FALSE) %>%  # If ties in abs, picks first
  ungroup() |>
  filter(str_detect(snp, '^rs')) |> 
  filter(str_detect(gene, '^ENSG')) |>
  mutate(key = paste(snp, gene, sep = '_')) |>
  select(key, beta_my = beta)
message('Pooled sig. eQTL after dup rm : ', nrow(pooled_my))
message('Any NAs in pooled sig. eQTL? ', anyNA(pooled_my))
message('\nPooled sig. eQTL tbl:\n')
pooled_my

message('\nLoading Fugita data ...\n')
fugita_lst <- list()
for (cell_type in fugita_cell_types) {

  file_path <- paste0(fugita_dir, 'celltype-eqtl-sumstats.', cell_type,'.tsv')
  
  # Read the file
  fugita_df <- read_delim(file_path, delim = '\t', col_names = TRUE, na = c("", "NA", "na", "-")) |>
    select(snp = snps, gene = gene_id, beta, celltype) |>
    mutate(beta = as.numeric(beta)) |>
    filter(str_detect(snp, '^rs')) |> 
    filter(str_detect(gene, '^ENSG')) |>
    filter(!is.na(beta))
  message('\nTotal SNP-gene pairs in Fugita ', cell_type, ': ' , nrow(fugita_df), '\n')
  
  # Append to list
  fugita_lst[[cell_type]] <- fugita_df
  rm(fugita_df)
}

# Combine all Fugita data
all_fugita_df <- bind_rows(fugita_lst)
message('Combined Fugita SNP-gene pairs before dup rm: ', nrow(all_fugita_df))
message('Are Fugita betas numeric after coercion? ', is.numeric(all_fugita_df$beta))


# For duplicates in Fugita, keep max abs(beta)
# pooled_fugita <- all_fugita_df %>%
#   group_by(snp, gene) %>%
#   slice_max(abs(beta), n = 1, with_ties = FALSE) %>%
#   ungroup() %>%
#   mutate(key = paste(snp, gene, sep = '_')) |>
#   select(key, beta_fugita = beta)

# Need to use DT as computation takes ages otherwise
setDT(all_fugita_df)

# De-duplicate: keep the row with max(abs(beta)) for each (snp, gene) group
# If ties in abs(beta), keeps the first occurrence (same as your with_ties = FALSE)
pooled_fugita_dt <- all_fugita_df[
  , .SD[which.max(abs(beta))], 
  by = .(snp, gene)
]

# Create the key and select only needed columns
pooled_fugita_dt <- pooled_fugita_dt[
  , .(key = paste(snp, gene, sep = '_'),
      beta_fugita = beta)                 # ← rename here
]

# Revert to tibble
pooled_fugita <- as_tibble(pooled_fugita_dt)

# Clean up
rm(all_fugita_df, pooled_fugita_dt)
gc() 


message('\nFugita SNP-gene pairs after dup rm: ', nrow(pooled_fugita))
message('Any NAs in pooled Fugita data? ', anyNA(pooled_fugita))
message('\nPooled Fugita tbl:\n')
pooled_fugita

# Join with your pooled data for paired betas
paired_betas <- pooled_my %>%
  GABAer_join(pooled_fugita, by = c('key')) %>%
  select(key, beta_my, beta_fugita)
message('Number of overlapping SNP-gene pairs: ', nrow(paired_betas))
message('Paired betas tbl:\n')
paired_betas

message('Writing tbl to:', out_file)
write_tsv(paired_betas, out_file)

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
