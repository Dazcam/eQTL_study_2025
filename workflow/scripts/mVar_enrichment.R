#--------------------------------------------------------------------------------------
#
#    mVAR enrichment tests
#
#--------------------------------------------------------------------------------------

## Info  ------------------------------------------------------------------------------

# 1. Extract FDR sig. genes (qval < 0.05) and assoc. gene-specific pval_nominal_threshold from perm file
# 2. Filt nominal file to inc. all SNPs in FDR sig. genes' in cis-region with pval_nominal < pval_nominal_threshold

##  Load Packages, functions and variables  -------------------------------------------
## Set up logging for smk  ------------------------------------------------------------
if (exists("snakemake")) {
  log_smk <- function() {
    if (exists("snakemake") & length(snakemake@log) != 0) {
      log <- file(snakemake@log[1][[1]], open = "wt")
      sink(log, append = TRUE)
      sink(log, append = TRUE, type = "message")
    }
  }
}
log_smk()
message('\n\nRunning mVar enrichment tests ...')

library(tidyverse)
library(readxl)
library(readr)

# Define paths and cell types
# Input and output paths
mvar_path <- snakemake@input[['mvars']] 

cell_types <- c(
  "ExN-DL", "ExN-DL-0", "ExN-DL-1", "ExN-DL-2", "ExN-DL-3",
  "ExN-UL", "ExN-UL-0", "ExN-UL-1", "ExN-UL-2"
)

# Function to read perm file for a cell type
read_perm <- function(ct) {
  file_path <- str_glue("../results/05TENSORQTL/tensorqtl_perm/{ct}_quantile_genPC_4_expPC_40/{ct}_quantile_perm.cis_qtl.txt.gz")
  df <- read_tsv(file_path, show_col_types = FALSE) |>
    filter(qval < 0.05) |> # FDR-significant only for overlaps/Fisher
    pull(variant_id) |> 
    unique() |> 
    as.character()
  
}

# Get all FDR sig SNP / gene pairs from tensorQTL nominal file 
read_all_sig_snps <- function(ct) {
  perm_path <- str_glue("../results/05TENSORQTL/tensorqtl_perm/{ct}_quantile_genPC_4_expPC_40/{ct}_quantile_perm.cis_qtl.txt.gz")
  nom_path <- str_glue("../results/10SMR/smr_input/{ct}/{ct}_nom.cis_qtl_pairs.tsv")  # Your original nominal path
  
  # Step 1: Read perm file, filter sig genes, extract thresholds
  message('\nReading permutation file for ', ct, ' ...')
  thresholds <- read_tsv(perm_path, show_col_types = FALSE) %>%
    filter(qval < 0.05) %>%
    select(phenotype_id, pval_nominal_threshold) %>%
    mutate(phenotype_id = as.character(phenotype_id))
  
  n_genes <- nrow(thresholds)
  if (n_genes == 0) {
    message("  → No significant eGenes")
    return(list(snps = character(0), gene_counts = integer(0), n_genes = 0))
  }
  
  # Step 2: Read nominal file (all SNP-gene pairs) and filter
  message('Found ', nrow(thresholds), ' significant eGenes. Reading nominal file ...')
  sig_pairs <- read_tsv(nom_path, show_col_types = FALSE) %>%
    select(phenotype_id, variant_id, pval_nominal) %>%
    mutate(phenotype_id = as.character(phenotype_id),
           variant_id   = as.character(variant_id),
           pval_nominal = as.numeric(pval_nominal)) %>%
    inner_join(thresholds, by = "phenotype_id") %>%
    filter(pval_nominal <= pval_nominal_threshold)
  
  # 3. Per-gene SNP counts
  gene_counts <- sig_pairs %>%
    count(phenotype_id, name = "n_snps_per_gene") %>%
    pull(n_snps_per_gene)
  
  # 4. Unique SNPs
  unique_snps <- unique(sig_pairs$variant_id)
  
  message(ct, ": ", length(unique_snps), " FDR-sig SNPs across ", n_genes, 
          " eGenes (median ", median(gene_counts), " SNPs/gene)")
  
  return(list(
    snps        = unique_snps,
    gene_counts = gene_counts,
    n_genes     = n_genes
  ))
}

results_list <- map(cell_types, read_all_sig_snps)
names(results_list) <- cell_types

# Summary table
message("\n=== Summary of FDR-sig. eQTLs across 9 developing ExN pops ===============")
eqtl_summary <- tibble(
  cell_type            = names(results_list),
  n_sig_genes          = map_int(results_list, "n_genes"),
  n_sig_snps           = map_int(results_list, ~ length(.x$snps)),
  mean_snps_per_gene   = round(n_sig_snps / n_sig_genes, 1),
  median_snps_per_gene = map_int(results_list, ~ as.integer(median(.x$gene_counts)))
) %>%
  arrange(desc(n_sig_genes)) |>
  print()
message("==========================================================================\n")

# Read mVar data


# Test 1: Raw overlaps
cat("\n======= PART 1: Raw overlaps (mVar ↔ FDR-significant eQTLs) =======\n\n")
mvar_df <- read_excel(mvar_path) %>%
  select(rsid = RSID, mDisorder) 
mvar_rsids <- pull(mvar_df, rsid) %>% unique()

cat("Number of unique mVar SNPs:", length(mvar_rsids), "\n\n")

# Extract SNP vectors from the optimized results_list
sig_snps_list <- map(results_list, "snps")
names(sig_snps_list) <- cell_types

# Per-cell-type overlaps
overlaps_per_ct <- tibble(cell_type = cell_types) %>%
  mutate(
    n_sig_snps_in_pop = map_int(sig_snps_list, length),
    overlap_count     = map_int(sig_snps_list, ~ length(intersect(mvar_rsids, .x))),
    overlap_percent   = round(overlap_count / length(mvar_rsids) * 100, 2)
  ) %>%
  arrange(desc(overlap_count)) %>%
  select(cell_type, overlap_count, overlap_percent, n_sig_snps_in_pop)

print(overlaps_per_ct, n = Inf)

# Combined (union across all 9 populations)
all_sig_snps_union <- reduce(sig_snps_list, union)
combined_overlap   <- length(intersect(mvar_rsids, all_sig_snps_union))
combined_percent   <- round(combined_overlap / length(mvar_rsids) * 100, 2)

cat("\nCombined (union across all 9 ExN populations):\n")
cat(sprintf("   → %d of %d mVar SNPs (%.2f%%) are FDR-significant eQTLs in ≥1 population\n",
            combined_overlap, length(mvar_rsids), combined_percent))
message("==========================================================================\n")


# Part 2: Fisher's exact test - Are mVars enriched amoung FDR sig. eQTLs per cell type?
get_total_tested_snps <- function(ct) {
  perm_path <- str_glue("../results/05TENSORQTL/tensorqtl_perm/{ct}_quantile_genPC_4_expPC_40/{ct}_quantile_perm.cis_qtl.txt.gz")
  
  read_tsv(perm_path, show_col_types = FALSE) %>%
    summarise(total_tested = sum(num_var, na.rm = TRUE)) %>%
    pull(total_tested)
}

# Apply to all cell types
total_tested_per_ct <- tibble(
  cell_type = cell_types,
  total_tested_snps = map_dbl(cell_types, get_total_tested_snps)
)



cat("\n=================== PART 2: Fisher's Exact Tests ====================\n")

message('\nUniverse (total SNPs) per cell type:\n')
print(total_tested_per_ct)

message('\nFishers Exact test results per cell type:\n')
# Function for Fisher's per cell type
fisher_results <- overlaps_per_ct %>%
  left_join(total_tested_per_ct, by = "cell_type") %>%
  mutate(
    a = overlap_count,                                      # mVar & sig eQTL
    b = length(mvar_rsids) - overlap_count,                 # mVar & not sig
    c = n_sig_snps_in_pop - overlap_count,                  # not mVar & sig
    d = total_tested_snps - (a + b + c),                    # not mVar & not sig
    or = (a/b) / (c/d),
    p_value = pmap_dbl(list(a, b, c, d), ~ {
      matrix <- matrix(c(..1, ..2, ..3, ..4), nrow = 2)
      fisher.test(matrix, alternative = "greater")$p.value
    })
  ) %>%
  select(cell_type, a, b, c, d, or, p_value, everything()) |>
  mutate(
    p_adjusted_bonf = pmin(p_value * 9, 1),
    rank = row_number(),
    p_adjusted_fdr  = p.adjust(p_value, method = "BH")
  ) %>%
  arrange(desc(a))

print(fisher_results %>%
        select(cell_type, a, or, p_value, p_adjusted_fdr, p_adjusted_bonf), n = Inf)

message("==========================================================================\n")


# Part 3: Are Pleiotropic mVars more active ExN pops than single-trait mVars?
cat("\n=== PART 3: Wicoon test  ===\n")

# Count in how many of the 9 cell types each mVar is an FDR-significant eQTL
mvar_activity <- tibble(rsid = mvar_rsids) %>%
  mutate(n_active_pops = rowSums(sapply(sig_snps_list, function(x) mvar_rsids %in% x))) %>%
  left_join(mvar_df %>% select(rsid, mDisorder), by = "rsid") %>%
  filter(!is.na(mDisorder)) %>%
  mutate(
    group = case_when(
      mDisorder >= 3                    ~ "Pleiotropic (≥3 GWAS traits)",
      mDisorder == 1                    ~ "Single-trait",
      TRUE                              ~ "Excluded (mDisorder = 0 or 2)"
    )
  ) %>%
  filter(group != "Excluded (mDisorder = 0 or 2)") |>
  arrange(desc(n_active_pops)) 

# Summary statistics
summary_pleio <- mvar_activity %>%
  group_by(group) %>%
  summarise(
    n_snps          = n(),
    median_pops     = median(n_active_pops),
    mean_pops       = round(mean(n_active_pops), 2),
    Q25             = quantile(n_active_pops, 0.25),
    Q75             = quantile(n_active_pops, 0.75),
    .groups = "drop"
  ) %>%
  mutate(across(c(median_pops, Q25, Q75), as.integer))

message('\nSumstats for Wilcoxon input:\n')
print(summary_pleio)

# One-sided Wilcoxon rank-sum test (pleiotropic > single-trait)
wilcox_result <- wilcox.test(n_active_pops ~ group,
                             data = mvar_activity,
                             alternative = "greater")

cat("\nWilcoxon rank-sum test (one-sided): pleiotropic > single-trait\n")
cat("W = ", wilcox_result$statistic[[1]], 
    "  p-value = ", format.pval(wilcox_result$p.value, digits = 3), "\n\n")

message("==========================================================================\n")

# Save everything
#write_csv(mvar_activity, "results/mVar_activity_across_9_ExN_pops.csv")
#write_csv(summary_pleio,   "results/mVar_pleiotropy_summary.csv")

# Quick violin/boxplot (optional – looks great in papers)
message("Plotting and saving to file ... ")
ggplot(mvar_activity, aes(x = group, y = n_active_pops, fill = group)) +
  geom_violin(trim = FALSE, alpha = 0.6) +
  geom_boxplot(width = 0.25, outlier.alpha = 0, alpha = 0.3) +
  scale_y_continuous(breaks = 0:10, limits = c(-0.5, 9.5)) +
  scale_fill_brewer(palette = "Set2") +
  labs(title = "Pleiotropic and single-trait operating mVars in ExN pops",
       x = "", y = "No. of ExN pops with mVar",
       fill = "") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none")
ggsave("../results/20MVARS/mVar_pleiotropy_vs_celltype_breadth.pdf", width = 8, height = 6)

enrichment_for_plot <- overlaps_per_ct %>%
  left_join(total_tested_per_ct, by = "cell_type") %>%
  mutate(prop_mvar = overlap_count / length(mvar_rsids) * 100) %>%
  select(cell_type, overlap_count, prop_mvar) %>%
  mutate(cell_type = fct_reorder(cell_type, prop_mvar))

ggplot(enrichment_for_plot, aes(x = cell_type, y = prop_mvar)) +
  geom_col(fill = "#1f78b4", alpha = 0.9, width = 0.7) +
  geom_text(aes(label = paste0(overlap_count, "/", length(mvar_rsids))),
            vjust = -0.5, size = 3.5) +
  labs(x = "", y = "% of mVars that are FDR-significant eQTLs",
       title = "Enrichment of mVars as eQTLs across 9 ExN populations") +
  theme_minimal(base_size = 13) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("../results/20MVARS/mVars_enrichment_bars.pdf", width = 7.5, height = 4.5)

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------