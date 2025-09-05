#--------------------------------------------------------------------------------------
#
#    pi1 enrichments for sig. eQTL with bulk data
#   
#--------------------------------------------------------------------------------------

# Questions: What is the lower threshold for overlaps to calc pi1?
# MG and OPC currently failing as too few SNPs to judge pi1

# pi1 process
# For each cell type:
#   1. Extract sig. eQTL @ FDR < 0.05
#   2. Extract pvals in ref eQTL data for SNPs ID'ed in step 1
#   3. Compute pi1 stat using the pvals from step 2 (i.e. ref data)

# -------------------------------------------------------------------------------------

library(tidyverse)
library(readr)
library(qvalue)
library(ggplot2)
library(cowplot)

# Define color palette
custom_palette <- c(
  "RG" = "#FF5959",
  "ExN-UL" = "#00B6EB",
  "InN" = "#3CBB75FF",
  "ExN-DL" = "#CEE5FD",
  "Endo-Peri" = "#B200ED",
  "MG" = "#F58231",
  "OPC" = "#FDE725FF"
)

# Define directories and cell types
cell_types <- c("ExN-UL", "ExN-DL", "InN", "RG", "MG", "OPC", "Endo-Peri")
exp_pc <- 40
gen_pc <- 4
in_dir <- "~/Desktop/eQTL_study_2025/workflow/reports/05TENSORQTL/logs/"
public_dir <- "~/Desktop/eQTL_study_2025/resources/public_datasets/"
tensorqtl_dir <- paste0(in_dir, 'results/05TENSORQTL/')
obrien_dir <- paste0(public_dir, 'obrien_2018/')

# Load full eQTL results from O'Brien and Bryois (placeholder paths)
# Replace with actual file paths once obtained
obrien_full <- read_tsv(paste0(obrien_dir, 'all_eqtls_gene.txt.gz')) %>%
  dplyr::select(variant_id, phenotype_id = gene_id, pval = pval_nominal, slope_ref = slope)
# bryois_full <- read_tsv(paste0(root_dir, 'resources/public_datasets/bryois_2022/full_eqtl_results.tsv')) %>%
#   dplyr::select(variant_id, phenotype_id, pval)

compute_pi1 <- function(my_eqtl, full_eqtl, ref_name, min_overlap = 50) {
  message("Computing pi1 statistic for ", ref_name)
  
  # Merge and extract p-values
  merged <- left_join(my_eqtl, full_eqtl, by = c("variant_id", "phenotype_id"))
  valid_merged <- merged[!is.na(merged$pval) & is.finite(merged$pval) & merged$pval >= 0 & merged$pval <= 1, ]
  valid_pvals <- valid_merged$pval
  overlap_count <- length(valid_pvals)
  
  message("Number of overlapping eQTLs: ", overlap_count)
  
  # Check if too few p-values
  if (overlap_count < min_overlap) {
    message("Too few overlapping eQTLs (<50) for pi1 calculation.")
    return(list(pi1 = NA, prop_replicating = NA, prop_same_direction = NA, overlap_count = overlap_count, pvals = valid_pvals))
  }
  
  # Attempt qvalue calculation
  qobj <- tryCatch(
    qvalue(p = valid_pvals, pi0.method = "smoother"),
    error = function(e) {
      message("Smoother method failed, switching to bootstrap")
      tryCatch(
        qvalue(p = valid_pvals, pi0.method = "bootstrap"),
        error = function(e2) {
          message("Bootstrap also failed: ", e2$message)
          NULL
        }
      )
    }
  )
  
  if (!is.null(qobj)) {
    pi1 <- 1 - qobj$pi0
    prop_replicating <- mean(qobj$qvalues < 0.05)
    rep_idx <- which(qobj$qvalues < 0.05)
    replicating <- valid_merged[rep_idx, ]
    prop_same_direction <- if (nrow(replicating) > 0) {
      mean(sign(replicating$slope_my) == sign(replicating$slope_ref), na.rm = TRUE)
    } else {
      NA
    }
    message("pi1 calculated: ", pi1)
  } else {
    pi1 <- NA
    prop_replicating <- NA
    prop_same_direction <- NA
  }
  
  return(list(pi1 = pi1, prop_replicating = prop_replicating, prop_same_direction = prop_same_direction, overlap_count = overlap_count, pvals = valid_pvals))
}

# Results storage
enrichment_results <- tibble()
overlapping_eqtls <- tibble()
plot_list <- list()

# Loop for O'Brien enrichment (pi1)
message("\n--- Starting O'Brien Enrichment Analysis ---")
for (cell_type in cell_types) {
  message("\nProcessing cell type: ", cell_type)
  
  my_eqtl <- read_tsv(paste0(in_dir, cell_type, '_genPC_', gen_pc, '_expPC_', exp_pc, '/', cell_type, '_perm.cis_qtl.txt.gz')) %>%
    filter(qval < 0.05) %>%
    dplyr::select(variant_id, phenotype_id, slope_my = slope)
  
  if (nrow(my_eqtl) == 0) {
    message("No significant eQTLs found for ", cell_type)
    next
  }
  
  pi1_result <- compute_pi1(my_eqtl, obrien_full, "O'Brien")
  
  enrichment_results <- bind_rows(enrichment_results, tibble(
    cell_type = cell_type,
    pi1_obrien = pi1_result$pi1,
    prop_replicating_obrien = pi1_result$prop_replicating,
    prop_same_direction_obrien = pi1_result$prop_same_direction,
    overlap_count_obrien = pi1_result$overlap_count
  ))
  
  # Collect overlapping eQTLs
  overlapping <- merge(my_eqtl, obrien_full, by = c("variant_id", "phenotype_id")) %>%
    mutate(cell_type = cell_type, ref_dataset = "O'Brien")
  overlapping_eqtls <- bind_rows(overlapping_eqtls, overlapping)
  
  # Generate plot if computation succeeded
  if (!is.na(pi1_result$pi1)) {
    plot_data <- data.frame(pvalue = pi1_result$pvals)
    max_density <- max(ggplot_build(ggplot(plot_data, aes(x = pvalue)) + geom_histogram(aes(y = ..density..), bins = 30))$data[[1]]$density)
    p <- ggplot(plot_data, aes(x = pvalue)) +
      geom_histogram(aes(y = ..density..), bins = 30, fill = custom_palette[cell_type], color = "black", alpha = 0.8) +
      annotate("text", x = 0.15, y = 0.90 * max_density, label = paste(
        "pi1 =", round(pi1_result$pi1, 2),
        "\nProp. replicating at 5% FDR:", round(pi1_result$prop_replicating, 2),
        "\nProp. same direction:", round(pi1_result$prop_same_direction, 2)
      ), hjust = 0) +
      labs(
        title = cell_type,
        x = "P-value",
        y = "Density"
      ) +
      theme_minimal(base_size = 12) +
      theme(
        plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
        axis.title = element_text(face = "bold"),
        axis.text = element_text(color = "black"),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 0.5),
        plot.caption = element_text(size = 10, hjust = 0)
      )
    plot_list[[cell_type]] <- p
  }
  
  message("O'Brien enrichment completed for ", cell_type)
}

#### Plotting
# Generate combined plot
if (length(plot_list) > 0) {
  combined_plot <- plot_grid(plotlist = plot_list, ncol = 3)
  print(combined_plot)
}


#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------



