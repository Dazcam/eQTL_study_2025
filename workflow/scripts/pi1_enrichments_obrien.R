#--------------------------------------------------------------------------------------
#
#    pi1 enrichments between sig. sn-eQTL & O'brien 2018 bulk brain eQTL
#   
#--------------------------------------------------------------------------------------

# Questions: What is the lower threshold for overlaps to calc pi1?
# MG and OPC currently failing as too few SNPs to judge pi1

# pi1 process
# For each cell type:
#   1. Extract sig. eQTL @ FDR < 0.05
#   2. Extract pvals in ref eQTL data for SNPs ID'ed in step 1
#   3. Compute pi1 stat using the pvals from step 2 (i.e. ref data)

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

message("\n\nCalculating pi1 enrichments: sig. sn-eQTL & O'brien 2018 bulk brain eQTL ...")

# -------------------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(tidyverse)
  library(qvalue)
  library(cowplot)
})

# Input and output paths
obrien_all_qtl <- snakemake@input[['obrien_all']]
obrien_top_qtl <- snakemake@input[['obrien_top']]
qtl_all <- snakemake@input[['qtl_all']]
qtl_top <- snakemake@input[['qtl_top']]
cell_type <- snakemake@wildcards[["cell_type"]]
exp_pc <- snakemake@wildcards[["exp_pc"]]
gen_pc <- snakemake@wildcards[["geno_pc"]]
norm_method <- snakemake@wildcards[["norm_method"]]

# Check variable assignment
message("\nVariables")
cat("============================")
tibble(
  variable = c("obrien_all_qtl", "obrien_top_qtl", "qtl_all", "qtl_top",
               "cell_type", "exp_pc", "gen_pc", "norm_method"),
  value    = c(obrien_all_qtl, obrien_top_qtl, qtl_all, qtl_top, 
               cell_type, exp_pc, gen_pc, norm_method)) |> 
  knitr::kable(format = "simple", align = "l") |>
  print()
message("\n============================\n")


# Functions  --------------------------------------------------------------------------
# Helper
compute_pi1 <- function(query_eqtl, full_eqtl, ref_name, min_overlap = 50) {
  message("Computing pi1 statistic for ", ref_name)
  
  merged <- left_join(query_eqtl, full_eqtl, by = c("variant_id", "phenotype_id"))
  valid_merged <- merged %>%
    filter(!is.na(pval), is.finite(pval), pval >= 0, pval <= 1)
  valid_pvals <- valid_merged$pval
  overlap_count <- length(valid_pvals)
  
  message("Number of overlapping eQTLs: ", overlap_count)
  
  if (overlap_count < min_overlap) {
    message("Too few overlapping eQTLs (<", min_overlap, ") for pi1 calculation.")
    return(list(pi1 = NA, prop_replicating = NA,
                prop_same_direction = NA,
                overlap_count = overlap_count,
                pvals = valid_pvals))
  }
  
  qobj <- tryCatch(
    qvalue(p = valid_pvals, pi0.method = "smoother"),
    error = function(e) {
      message("Smoother failed, trying bootstrap ...")
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
    message("pi1 calculated: ", round(pi1, 3))
  } else {
    pi1 <- NA
    prop_replicating <- NA
    prop_same_direction <- NA
  }
  
  return(list(pi1 = pi1,
              prop_replicating = prop_replicating,
              prop_same_direction = prop_same_direction,
              overlap_count = overlap_count,
              pvals = valid_pvals))
}

#------------------------------------------
# Main function
run_pi1_enrichment <- function(cell_type, obrien_all_qtl, obrien_top_qtl,
                               qtl_all, qtl_top, output) {
  
  # Load datasets
  message("Loading O'Brien all eQTL ...")
  obrien_full <- read_tsv(obrien_all_qtl, show_col_types = FALSE) %>%
    select(variant_id, phenotype_id = gene_id,
           pval = pval_nominal, slope_ref = slope)
  
  message("Loading O'Brien top eQTL (FDR < 0.05) ...")
  obrien_top <- read_tsv(obrien_top_qtl, show_col_types = FALSE) %>%
    filter(qval < 0.05) %>%
    select(variant_id, phenotype_id = gene_id, slope_my = slope)
  
  message("Loading sig. eQTL: ", cell_type)
  query_eqtl <- read_tsv(qtl_top, show_col_types = FALSE) %>%
    filter(qval < 0.05) %>%
    select(variant_id, phenotype_id, slope_my = slope)
  
  message("Loading all eQTL: ", cell_type)
  full_cell <- read_tsv(qtl_all, show_col_types = FALSE) %>%
    select(variant_id, phenotype_id,
           pval = pval_nominal, slope_ref = slope)
  
  # Initialize results tibble
  enrichment_results <- tibble(
    cell_type = cell_type,
    pi1_obrien = NA_real_,
    prop_replicating_obrien = NA_real_,
    prop_same_direction_obrien = NA_real_,
    overlap_count_obrien = NA_integer_,
    pi1_cell_ref = NA_real_,
    prop_replicating_cell_ref = NA_real_,
    prop_same_direction_cell_ref = NA_real_,
    overlap_count_cell_ref = NA_integer_
  )
  
  # Forward: cell sig > O’Brien full
  if (nrow(query_eqtl) > 0) {
    message("\n--- Forward enrichment: cell sig → O'Brien full ---")
    pi1_result <- compute_pi1(query_eqtl, obrien_full, "O'Brien")
    enrichment_results <- enrichment_results %>%
      mutate(
        pi1_obrien = pi1_result$pi1,
        prop_replicating_obrien = pi1_result$prop_replicating,
        prop_same_direction_obrien = pi1_result$prop_same_direction,
        overlap_count_obrien = pi1_result$overlap_count
      )
  } else {
    message("No significant eQTLs found for ", cell_type, " (forward).")
  }
  
  # Reverse: O’Brien sig > cell full
  if (nrow(full_cell) > 0) {
    message("\n--- Reverse enrichment: O'Brien sig → cell full ---")
    pi1_result <- compute_pi1(obrien_top, full_cell, cell_type)
    enrichment_results <- enrichment_results %>%
      mutate(
        pi1_cell_ref = pi1_result$pi1,
        prop_replicating_cell_ref = pi1_result$prop_replicating,
        prop_same_direction_cell_ref = pi1_result$prop_same_direction,
        overlap_count_cell_ref = pi1_result$overlap_count
      )
  } else {
    message("No full eQTL data found for ", cell_type, " (reverse).")
  }
  
  # Always write result (even if NA)
  write_rds(enrichment_results, output)
  message("\nResults written to: ", output)
}
#------------------------------------------

#------------------------------------------
# Run with Snakemake params
if (exists("snakemake")) {
  run_pi1_enrichment(
    cell_type       = snakemake@wildcards[["cell_type"]],
    obrien_all_qtl  = snakemake@input[["obrien_all"]],
    obrien_top_qtl  = snakemake@input[["obrien_top"]],
    qtl_all         = snakemake@input[["qtl_all"]],
    qtl_top         = snakemake@input[["qtl_top"]],
    output          = snakemake@output[[1]]
  )
}
#------------------------------------------



# # Load O'Brien full eQTL results
# message("Loading O'Brien all eQTL ...")
# obrien_full <- read_tsv(obrien_all_qtl) %>%
#   dplyr::select(variant_id, phenotype_id = gene_id, pval = pval_nominal, slope_ref = slope)
# 
# # Load O'Brien top eQTL (FDR < 0.05)
# message("Loading O'Brien top eQTL ...")
# obrien_top <- read_tsv(obrien_top_qtl) %>%
#   filter(qval < 0.05) %>%
#   select(variant_id, phenotype_id = gene_id, slope_my = slope)
# 
# compute_pi1 <- function(query_eqtl, full_eqtl, ref_name, min_overlap = 50) {
#   message("Computing pi1 statistic for ", ref_name)
#   
#   # Merge and extract p-values
#   merged <- left_join(query_eqtl, full_eqtl, by = c("variant_id", "phenotype_id"))
#   valid_merged <- merged[!is.na(merged$pval) & is.finite(merged$pval) & merged$pval >= 0 & merged$pval <= 1, ]
#   valid_pvals <- valid_merged$pval
#   overlap_count <- length(valid_pvals)
#   
#   message("Number of overlapping eQTLs: ", overlap_count)
#   
#   # Check if too few p-values
#   if (overlap_count < min_overlap) {
#     message("Too few overlapping eQTLs (<50) for pi1 calculation.")
#     return(list(pi1 = NA, prop_replicating = NA, prop_same_direction = NA, overlap_count = overlap_count, pvals = valid_pvals))
#   }
#   
#   # Attempt qvalue calculation
#   qobj <- tryCatch(
#     qvalue(p = valid_pvals, pi0.method = "smoother"),
#     error = function(e) {
#       message("Smoother method failed, switching to bootstrap")
#       tryCatch(
#         qvalue(p = valid_pvals, pi0.method = "bootstrap"),
#         error = function(e2) {
#           message("Bootstrap also failed: ", e2$message)
#           NULL
#         }
#       )
#     }
#   )
#   
#   if (!is.null(qobj)) {
#     pi1 <- 1 - qobj$pi0
#     prop_replicating <- mean(qobj$qvalues < 0.05)
#     rep_idx <- which(qobj$qvalues < 0.05)
#     replicating <- valid_merged[rep_idx, ]
#     prop_same_direction <- if (nrow(replicating) > 0) {
#       mean(sign(replicating$slope_my) == sign(replicating$slope_ref), na.rm = TRUE)
#     } else {
#       NA
#     }
#     message("pi1 calculated: ", pi1)
#   } else {
#     pi1 <- NA
#     prop_replicating <- NA
#     prop_same_direction <- NA
#   }
#   
#   return(list(pi1 = pi1, prop_replicating = prop_replicating, prop_same_direction = prop_same_direction, overlap_count = overlap_count, pvals = valid_pvals))
# }
# 
# # Results storage
# enrichment_results <- tibble()
# overlapping_eqtls <- tibble()
# plot_list <- list()
# 
# # Loop for O'Brien enrichment (pi1) - forward (sig. sn-cell-types as query, O'Brien full as ref)
# message("\n--- Starting O'Brien pi1 Enrichment Analysis (forward: cell sig query, O'Brien ref) ---")
# 
# message("\nProcessing cell type: ", cell_type)
# 
# query_eqtl <- read_tsv(qtl_top) %>%
#   filter(qval < 0.05) %>%
#   dplyr::select(variant_id, phenotype_id, slope_my = slope)
# 
# if (nrow(query_eqtl) == 0) {
#   message("No significant eQTLs found for ", cell_type)
#   next
# }
# 
# pi1_result <- compute_pi1(query_eqtl, obrien_full, "O'Brien")
# 
# enrichment_results <- bind_rows(enrichment_results, tibble(
#   cell_type = cell_type,
#   pi1_obrien = pi1_result$pi1,
#   prop_replicating_obrien = pi1_result$prop_replicating,
#   prop_same_direction_obrien = pi1_result$prop_same_direction,
#   overlap_count_obrien = pi1_result$overlap_count
# ))
# 
# # Collect overlapping eQTLs
# overlapping <- merge(query_eqtl, obrien_full, by = c("variant_id", "phenotype_id")) %>%
#   mutate(cell_type = cell_type, ref_dataset = "O'Brien")
# overlapping_eqtls <- bind_rows(overlapping_eqtls, overlapping)
# 
# # Generate plot if computation succeeded - Move to report file
# # if (!is.na(pi1_result$pi1)) {
# #   plot_data <- data.frame(pvalue = pi1_result$pvals)
# #   max_density <- max(ggplot_build(ggplot(plot_data, aes(x = pvalue)) + geom_histogram(aes(y = ..density..), bins = 30))$data[[1]]$density)
# #   p <- ggplot(plot_data, aes(x = pvalue)) +
# #     geom_histogram(aes(y = ..density..), bins = 30, fill = custom_palette[cell_type], color = "black", alpha = 0.8) +
# #     annotate("text", x = 0.15, y = 0.90 * max_density, label = paste(
# #       "pi1 =", round(pi1_result$pi1, 2),
# #       "\nProp. replicating at 5% FDR:", round(pi1_result$prop_replicating, 2),
# #       "\nProp. same direction:", round(pi1_result$prop_same_direction, 2)
# #     ), hjust = 0) +
# #     labs(
# #       title = cell_type,
# #       x = "P-value",
# #       y = "Density"
# #     ) +
# #     theme_minimal(base_size = 12) +
# #     theme(
# #       plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
# #       axis.title = element_text(face = "bold"),
# #       axis.text = element_text(color = "black"),
# #       panel.grid.minor = element_blank(),
# #       panel.border = element_rect(color = "black", fill = NA, size = 0.5),
# #       plot.caption = element_text(size = 10, hjust = 0)
# #     )
# #   plot_list[[cell_type]] <- p
# # }
# 
# message("O'Brien enrichment completed for ", cell_type)
# 
# # Reverse analysis: O'Brien sig as query, cell full as ref
# message("\n--- Starting O'Brien pi1 Enrichment Analysis (reverse: O'Brien sig query, cell full ref) ---")
# plot_list_reverse <- list()
# overlapping_eqtls_reverse <- tibble()
#   
# message("\nProcessing cell type as ref: ", cell_type)
# 
# full_cell_path <- paste0(qtl_all)
# full_cell <- read_tsv(full_cell_path) %>%
#   select(variant_id, phenotype_id, pval = pval_nominal, slope_ref = slope)
# 
# if (nrow(full_cell) == 0) {
#   message("No full eQTL data for ", cell_type)
#   next
# }
# 
# pi1_result <- compute_pi1(obrien_top, full_cell, cell_type)
# 
# # Add reverse results to the enrichment table (update existing rows)
# enrichment_results <- enrichment_results %>%
#   mutate(
#     pi1_cell_ref = ifelse(cell_type == !!cell_type, pi1_result$pi1, NA),
#     prop_replicating_cell_ref = ifelse(cell_type == !!cell_type, pi1_result$prop_replicating, NA),
#     prop_same_direction_cell_ref = ifelse(cell_type == !!cell_type, pi1_result$prop_same_direction, NA),
#     overlap_count_cell_ref = ifelse(cell_type == !!cell_type, pi1_result$overlap_count, NA)
#   )
# 
# # Collect overlapping eQTLs for reverse
# overlapping <- merge(obrien_top, full_cell, by = c("variant_id", "phenotype_id")) %>%
#   mutate(cell_type = cell_type, ref_dataset = "Cell")
# overlapping_eqtls_reverse <- bind_rows(overlapping_eqtls_reverse, overlapping)
# 
# # Generate plot if computation succeeded - Move to report file
# # if (!is.na(pi1_result$pi1)) {
# #   plot_data <- data.frame(pvalue = pi1_result$pvals)
# #   max_density <- max(ggplot_build(ggplot(plot_data, aes(x = pvalue)) + geom_histogram(aes(y = ..density..), bins = 30))$data[[1]]$density)
# #   p <- ggplot(plot_data, aes(x = pvalue)) +
# #     geom_histogram(aes(y = ..density..), bins = 30, fill = custom_palette[cell_type], color = "black", alpha = 0.8) +
# #     annotate("text", x = 0.15, y = 0.90 * max_density, label = paste(
# #       "pi1 =", round(pi1_result$pi1, 2),
# #       "\nProp. replicating at 5% FDR:", round(pi1_result$prop_replicating, 2),
# #       "\nProp. same direction:", round(pi1_result$prop_same_direction, 2)
# #     ), hjust = 0) +
# #     labs(
# #       title = paste(cell_type, "(cell ref)"),
# #       x = "P-value",
# #       y = "Density"
# #     ) +
# #     theme_minimal(base_size = 12) +
# #     theme(
# #       plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
# #       axis.title = element_text(face = "bold"),
# #       axis.text = element_text(color = "black"),
# #       panel.grid.minor = element_blank(),
# #       panel.border = element_rect(color = "black", fill = NA, size = 0.5),
# #       plot.caption = element_text(size = 10, hjust = 0)
# #     )
# #   plot_list_reverse[[cell_type]] <- p
# #}
# 
# message("O'Brien enrichment (reverse) completed for ", cell_type)
# 
# 
# write_rds(enrichment_results, paste0(output))
# 
# 
# 
# #### Plotting - move to report file
# # Generate combined plot
# # combined_plot <- plot_grid(plotlist = plot_list, ncol = 3)
# # print(combined_plot)
# # 
# # 
# # # Generate combined plot for Bryois reverse
# # combined_plot_reverse <- plot_grid(plotlist = plot_list_reverse, ncol = 3)
# # print(combined_plot_reverse)


#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------



