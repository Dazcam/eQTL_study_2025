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

# Load required libraries
library(tidyverse)
library(data.table)

# Define variables
cell_types <- c(
  "ExN-UL", "ExN-UL-0", "ExN-UL-1", "ExN-UL-2",
  "ExN-DL", "ExN-DL-0", "ExN-DL-1", "ExN-DL-2", "ExN-DL-3",
  "RG", "RG-0", "RG-1", "RG-2", "RG-3",
  "InN", "InN-0", "InN-1",
  "Endo-Peri",
  "OPC",
  "MG"
)

exp_PCs <- 40
gen_PCs <- 4
norm_method <- 'quantile' 
in_dir <- '../results/05TENSORQTL/tensorqtl_perm/'

fugita_dir <- '../resources/public_datasets/fugita_2024/'
fugita_cell_types <- c('Ast', 'Exc', 'Inh', 'Oli', 'OPC', 'Mic', 'End')

# Step 1: Pool unique sig. SNP-gene pairs from eQTL data, max(abs(beta)) for dups
my_data_list <- list()

message('Loading sig.eQTL data ...\n')
for (cell_type in cell_types) {
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
message('Are Fugita betas numeric: ')


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
pooled_fugita_dt[
  , key := paste(snp, gene, sep = '_')
][
  , .(key, beta = beta_fugita)
]

# Revert to tibble
pooled_fugita <- as_tibble(pooled_fugita_dt)

# Clean up
rm(all_fugita_df, pooled_fugita_dt)
gc() 


message('Fugita SNP-gene pairs after dup rm: ', nrow(pooled_fugita))
message('Any NAs in pooled Fugita data? ', anyNA(pooled_fugita))
message('Pooled Fugita tbl:\n')
pooled_fugita

# Join with your pooled data for paired betas
paired_betas <- pooled_my %>%
  inner_join(pooled_fugita, by = c('key')) %>%
  select(key, beta_my, beta_fugita)
message('Number of overlapping SNP-gene pairs: ', nrow(paired_betas))
message('Paired betas tbl:\n')
paired_betas

message('Writing tbl to:', out_file)
write_tsv(paired_betas, out_file)

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------

# # Run correlation analysis (Pearson)
# cor_result <- cor.test(paired_betas$beta_my, paired_betas$beta_fugita, method = 'pearson')
# print(cor_result) 
# 
# cor_r <- round(cor_result$estimate, 3)
# cor_p <- ifelse(cor_result$p.value < 1e-3, "< 0.001", round(cor_result$p.value, 3))
# cor_label <- paste0("r = ", cor_r, ", p = ", cor_p)
# 
# # Generate publication-ready plot
# p <- ggplot(paired_betas, aes(x = beta_fugita, y = beta_my)) +
#   geom_point(alpha = 0.6, color = 'steelblue', size = 1.5) +
#   geom_smooth(method = 'lm', color = 'red', se = TRUE, linewidth = 1) +
#   annotate("text", x = Inf, y = Inf, label = cor_label, hjust = 1.02, vjust = 1.02,
#            size = 4, fontface = "bold") +
#   labs(
#     x = 'Beta (Query)',
#     y = 'Beta (Reference)'
#   ) +
#   theme_minimal() +  # Publication-ready theme from ggpubr
#   theme(
#     plot.title = element_text(hjust = 0.5, face = 'bold'),
#     axis.title = element_text(face = 'bold')
#   )
# 
# # Display plot
# print(p)
# 
# # Save plot to file (e.g., for publication)
# ggsave('eqtl_beta_correlation_plot.png', p, width = 8, height = 6, dpi = 300, bg = 'white')
# ggsave('eqtl_beta_correlation_plot.pdf', p, width = 8, height = 6)  