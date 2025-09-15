#--------------------------------------------------------------------------------------
#
#    Plot TSS densities and munge eQTL log files
#
#--------------------------------------------------------------------------------------

# Cell types processed together here rather than independent jobs
# Most fo this has moved to TensorQTL Rmd reprt file - consider removing script

## Load libraries and variables -------------------------------------------------------
library(tidyverse)

if (exists("snakemake")) {
  log_smk <- function() {
    if (exists("snakemake") & length(snakemake@log) != 0) {
      log <- file(snakemake@log[1][[1]], open = "wt")
      sink(log, append = TRUE)
      sink(log, append = TRUE, type = "message")
    }
  }
  log_smk()
  message('Generating TSS enrichment plots and tensorQTL summsary tables remotely ...')
  root_dir <- snakemake@params[["root_dir"]]
  cell_types <- snakemake@params[["cell_types"]]  # <- Get cell types from config
} else {
  root_dir <- '~/Desktop/eQTL_study_2025/'
  cell_types <- c('RG', 'ExN-1', 'InN', 'ExN-2', 'Endo-Peri', 'MG', 'OPC')  # fallback default
}

tensorqtl_dir <- paste0(root_dir, 'results/04TENSORQTL/')
fig_dir <- paste0(root_dir, 'results/figs/')
table_dir <- paste0(root_dir, 'results/tables/')

## Create TSS density plot  -----------------------------------------------------------
plt_list <- list()
for (cell_type in cell_types) {
  
  message('Reading eQTL data for: ', cell_type)
  eqtl_tbl <- read_tsv(paste0(tensorqtl_dir, cell_type, '/', cell_type, '_perm.cis_qtl.txt.gz')) %>%
    dplyr::filter(!is.na(start_distance)) %>%
    mutate(distance_kb = start_distance / 1000) %>%
    dplyr::filter(qval < 0.05)
  
  message('eQTL detected at FDR < 0.05: ', nrow(eqtl_tbl))
  
  if (nrow(eqtl_tbl) > 0) {
    
    message('Plotting ...')
    
    eqtl_count <- nrow(eqtl_tbl)
    
    eqtl_density_plt <- eqtl_tbl  %>%
      ggplot(aes(x = distance_kb)) +
      geom_density() +
      geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
      annotate("text", x = 100, y = Inf, 
               label = paste("eQTLs:", eqtl_count), 
               hjust = 0, vjust = 1.5, size = 4, fontface = "bold") +
      scale_x_continuous(limits = c(-600, 600), 
                         breaks = seq(-600, 600, by = 200),
                         labels = seq(-600, 600, by = 200)) +
      labs(x = "Distance to TSS (kb)",
           y = "Density",
           title = cell_type) +
      theme_minimal() +
      theme(
        plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        axis.text = element_text(color = "black"),
        axis.title = element_text(face = "bold")
      )
    
    plt_list[[cell_type]] <- eqtl_density_plt
    
  } else {
    
    message('Skipped plotting.')
    
  }
  
}

message('Combining TSS plots ...')
combined_plot <- cowplot::plot_grid(plotlist = plt_list)
message('Plotting done.\n\n')


summary_tbl <- map_dfr(cell_types, function(cell_type) {
  
  message('Reading eQTL log file for: ', cell_type)
  log_file <- paste0(tensorqtl_dir, cell_type, '/', cell_type, '_perm.tensorQTL.cis.log')
  
  lines <- readLines(log_file, warn = FALSE)
  
  gene_line <- str_subset(lines, 'phenotypes tested:')
  sample_line <- str_subset(lines, 'samples')
  covariate_line <- tail(str_subset(lines, 'covariates'), 1)
  covariate_line <- as.integer(str_extract(covariate_line, '\\d+'))
  variants_line <- str_subset(lines, 'variants')
  qtl_line  <- str_subset(lines, 'QTL phenotypes @ FDR 0.05:')
  window_line <- str_subset(lines, 'cis-window:') 
  run_line <- str_subset(lines, 'Time elapsed:') 
  beta_line <- str_subset(lines, 'Correlation between Beta-approximated and empirical p-values:') 
  min_p_line <- str_subset(lines, 'min p-value threshold @ FDR 0.05:')
  prop_sig_line <- str_subset(lines, 'Proportion of significant phenotypes')
  
  # Extract number after colon (and optional whitespace) at end of line
  num_genes <- as.integer(str_extract(gene_line, '\\d+$'))
  num_samples <- as.integer(str_extract(sample_line, '\\d+'))
  covariates <- as.integer(str_extract(covariate_line, '\\d+$'))
  variants <- as.integer(str_extract(variants_line, '\\d+'))
  num_qtls  <- as.integer(str_extract(qtl_line, '\\d+$'))
  cis_window  <- as.integer(str_extract(window_line, '\\d+'))
  run_time  <- as.integer(str_extract(run_line, '\\d+'))
  beta_cor  <- as.numeric(str_extract(beta_line, '[0-9\\.eE-]+$'))
  min_p  <- as.numeric(str_extract(min_p_line, '[0-9\\.eE-]+$'))
  prop_sig <- as.integer(str_extract(prop_sig_line, '\\d+$'))
  
  message('Collating sumstats for: ', cell_type, ' ...')
  tibble(
    cell_type = cell_type,
    genes_tested = num_genes,
    samples = num_samples,
    covariates = covariates,
    variants = variants,
    qtl_fdr_0.05 = num_qtls,
    cis_window_kb = cis_window,
    beta_cor = beta_cor,
    min_p = min_p,
    prop_sig = prop_sig,
    run_time_min = run_time
  ) |>
    mutate(min_p = formatC(min_p, format = "e", digits = 2)) 
})

# Save files
message('Saving TSS plot and sumstats table ...')
ggsave(filename = paste0(fig_dir, "tss_density_plot_combined.png"), 
       plot = combined_plot, width = 12, height = 8, dpi = 300)
write_tsv(summary_tbl, paste0(table_dir, "tensorqtl_summary_table.tsv"))
message('All done.')

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
