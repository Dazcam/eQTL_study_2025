#--------------------------------------------------------------------------------------
#
#    Generate eQTL supplementry plts for manuscript
#
#--------------------------------------------------------------------------------------

# 1: Genotype PC elbow plt
# 2: eQTL Expression QC plt L1

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

message("\n\nGenerating supplementry plot for the manuscript ...")

# -------------------------------------------------------------------------------------
library(tidyverse)
library(cowplot)
library(ggrepel)


# ----- Set variables
geno_dir <- snakemake@params[['geno_dir']]
expr_dir <- snakemake@params[['expr_dir']] 
out_file <- snakemake@output[[1]]
out_dir <- dirname(out_file)

message('Smk variables set to:')
message('  geno_dir: ', geno_dir)
message('  expr_dir: ', expr_dir)
message('  out_dir: ', out_dir)
message('  out_file (tracked by snakemake): ', out_file)

cell_types <- c("Glu-UL", "Glu-DL", "GABA", "NPC", "MG", "OPC", "Endo-Peri")
cell_types_L2 <- c("Glu-UL-0", "Glu-UL-1", "Glu-UL-2",
                   "Glu-DL-0", "Glu-DL-1", "Glu-DL-2", 
                   "NPC-0", "NPC-1", "NPC-2", 
                   "GABA-0", "GABA-1", "GABA-2")
exp_PCs <- c(10, 20, 30, 40, 50)
gen_PCs <- 4
norm_methods <- 'quantile'

custom_palette <- c(
  'Glu-UL' = '#4363d8',
  'Glu-DL' = '#00B6EB',
  'NPC' = '#FF5959',
  'GABA' = '#3CBB75FF',
  'Endo-Peri' = '#B200ED',
  'MG' = '#F58231',
  'OPC' = '#FDE725FF'
)

custom_palette_L2 <- c(
  "NPC-0" = '#FF5959',
  "NPC-1" = '#FF5959',
  "NPC-2" = '#FF5959',
  "Glu-UL-0" = '#4363d8',
  "Glu-UL-1" = '#4363d8',
  "Glu-UL-2" = '#4363d8',
  "Glu-DL-0" = '#00B6EB',
  "Glu-DL-1" = '#00B6EB',
  "Glu-DL-2" = '#00B6EB',
  "GABA-0" = '#3CBB75FF',
  "GABA-1" = '#3CBB75FF',
  "GABA-2" = '#3CBB75FF'
)

base_theme <- theme(
  axis.text.y = element_text(size = 10),
  axis.text.x = element_text(size = 10),
  panel.grid.major.x = element_blank(),
  panel.grid.minor = element_blank(),
  panel.background   = element_rect(fill = "white", colour = NA),
  plot.background    = element_rect(fill = "white", colour = NA),
  strip.text = element_text(face = "bold"),
  legend.position = "none",
  plot.title = element_text(hjust = 0.5, face = "bold"),
  panel.spacing.x = unit(2, "lines"),
)

# ----- Load data
message('Loading eigenval file ...')
eigenval <- read_table(paste0(geno_dir, "pca.eigenval"), col_names = "eigenval")


# ----- Geno PC Elbow plt
# Calc proportion and cumulative variance explained
message('Generating Elbow plot ...')
variance <- eigenval %>%
  mutate(PC = row_number(),
         variance_explained = eigenval / sum(eigenval) * 100,
         cumulative = cumsum(variance_explained))

# Elbow plot
geno_plt <- ggplot(variance, aes(x = PC, y = variance_explained)) +
  geom_line(color = "steelblue", linewidth = 1) +
  geom_point(color = "steelblue", size = 2) +
  theme_minimal(base_size = 14) +
  theme(panel.grid.major.x = element_blank()) +
  labs(
    x = "Principal Component",
    y = "Variance Explained (%)",
    title = element_blank()
  ) +
  scale_x_continuous(
    breaks = seq(0, 10, by = 2),
    labels = seq(0, 10, by = 2)
  )

message('Saving Elbow plot ...')
ggsave(
  filename = paste0(out_dir, "/genotype_elbow_plt.pdf"),
  plot = geno_plt,
  width = 5,
  height = 5,
  units = "in",
  device = "pdf",
  dpi = 300
)

# ----- eQTL expression PC plt
message('Loading eQTL expression data ...')
summarise_logs <- function(cell_types, in_dir, expPCs = exp_PCs, 
                           genPC = gen_PCs, norm_method = norm_methods) {
  # Generate all combinations of cell types and expPC values
  params_grid <- expand.grid(cell_type = cell_types, expPC = expPCs, norm_method = norm_methods,
                             stringsAsFactors = FALSE)
  
  # Function to process a single log file
  process_log <- function(cell_type, expPC, norm_method) {
    message('Reading eQTL log file for: ', cell_type, ', ', norm_method , 
            ', with genPC_', genPC, ' expPC_', expPC)
    
    # Construct the log file path
    log_dir <- paste0(in_dir, cell_type, '_', norm_method, '_genPC_', genPC, '_expPC_', expPC, '/')
    log_file <- paste0(log_dir, cell_type, '_', norm_method, '_perm.tensorQTL.cis.log')
    
    # Check if file exists
    if (!file.exists(log_file)) {
      message('Log file not found: ', log_file)
      return(NULL)
    }
    
    # Read and process the log file
    lines <- readLines(log_file, warn = FALSE)
    
    gene_line <- str_subset(lines, 'phenotypes tested:')
    sample_line <- str_subset(lines, 'samples')
    covariate_line <- tail(str_subset(lines, 'covariates'), 1) # need to fix this
    variants_line <- str_subset(lines, 'variants')
    qtl_line <- str_subset(lines, 'QTL phenotypes @ FDR 0.05:')
    window_line <- str_subset(lines, 'cis-window:')
    run_line <- str_subset(lines, 'Time elapsed:')
    beta_line <- str_subset(lines, 'Correlation between Beta-approximated and empirical p-values:')
    min_p_line <- str_subset(lines, 'min p-value threshold @ FDR 0.05:')
    prop_sig_line <- str_subset(lines, 'Proportion of significant phenotypes')
    
    # Extract numbers
    num_genes <- as.integer(str_extract(gene_line, '\\d+$'))
    num_samples <- as.integer(str_extract(sample_line, '\\d+'))
    covariates <- as.integer(str_extract(covariate_line, '\\d+$'))
    variants <- as.integer(str_extract(variants_line, '\\d+'))
    num_qtls <- as.integer(str_extract(qtl_line, '\\d+$'))
    cis_window <- as.integer(str_extract(window_line, '\\d+'))
    run_time <- as.integer(str_extract(run_line, '\\d+'))
    beta_cor <- as.numeric(str_extract(beta_line, '[0-9\\.eE-]+$'))
    min_p <- as.numeric(str_extract(min_p_line, '[0-9\\.eE-]+$'))
    prop_sig <- as.integer(str_extract(prop_sig_line, '\\d+$'))
    
    message('Collating sumstats for: ', cell_type, ' with genPC_', genPC, ' expPC_', expPC, ' ...')
    
    # Create tibble with all relevant fields
    tibble(
      cell_type = cell_type,
      genPC = genPC,
      expPC = expPC,
      norm_methods = norm_method,
      genes_tested = num_genes,
      samples = num_samples,
      covariates = covariates,
      variants = variants,
      qtl_fdr_0.05 = num_qtls,
      cis_window_kb = cis_window,
      beta_cor = beta_cor,
      min_p = formatC(min_p, format = "e", digits = 2),
      prop_sig = prop_sig,
      run_time_min = run_time
    )
  }
  
  # Apply the processing function to all combinations
  summary_list <- pmap(params_grid, function(cell_type, expPC, norm_method) {
    process_log(cell_type, expPC, norm_method)
  })
  
  # Identify which elements are non-null
  valid_idx <- !sapply(summary_list, is.null)
  table(valid_idx)  # how many TRUE vs FALSE
  length(summary_list)
  sum(valid_idx)
  
  # Keep only valid entries
  summary_list <- summary_list[valid_idx]
  
  # Name the list elements for clarity
  names(summary_list) <- paste0(
    params_grid$cell_type, '_', 
    params_grid$norm_method, '_genPC_', genPC, 
    '_expPC_', params_grid$expPC
  )[valid_idx]
  
  # Return the list of summary tibbles
  return(summary_list)
}

summary_L1_list <- summarise_logs(cell_types, expr_dir)
summary_L1_tbl <- bind_rows(summary_L1_list)

summary_L2_list <- summarise_logs(cell_types_L2, expr_dir)
summary_L2_tbl <- bind_rows(summary_L2_list)

# Plot
message('Generating expression L1 PC plot ...')
expr_L1_plt <- ggplot(summary_L1_tbl %>% filter(norm_methods == 'quantile'),
       aes(x = expPC, y = qtl_fdr_0.05, color = cell_type, group = cell_type)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  scale_color_manual(values = custom_palette) +
  labs(
    x = "Expression PCs",
    y = "eGene count (FDR < 0.05)",
    color = "Level 1 Cell Type"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.major.x = element_blank(),
    plot.title = element_blank(),
    legend.position = "top",
    legend.title = element_blank(),
    plot.margin = margin(10, 10, 10, 10)
  ) + ylim(0, 2100)

# Plot level 2
message('Generating expression L2 PC plot ...')
expr_L2_plt <- ggplot(summary_L2_tbl %>% filter(norm_methods == 'quantile'),
                   aes(x = expPC, y = qtl_fdr_0.05, color = cell_type, group = cell_type)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  scale_color_manual(values = custom_palette_L2) +
  labs(
    x = "Expression PCs",
    y = "eGene count (FDR < 0.05)",
    color = "Level 2 Cell Type"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.major.x = element_blank(),
    plot.title = element_blank(),
    legend.position = "top",
    legend.title = element_blank(),
    plot.margin = margin(10, 10, 10, 10)
  ) + ylim(0, 2100)

# Combine plots
final_plt <- plot_grid(expr_L1_plt, expr_L2_plt, labels = 'AUTO', 
                       ncol = 2, label_size = 24)


message('Saving expression PC plot ...')
ggsave(
  filename = out_file,
  plot = final_plt,
  width = 10,
  height = 5,
  units = "in",
  device = "pdf",
  dpi = 300
)

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------