# --------------------------------------------------------------------------------------
#
#    Proportion overlap Wen
#
#--------------------------------------------------------------------------------------

#  1. Extract eQTL (top SNP) from all eGenes in L1 cell types (will have dup genes)
#  2. Calc proportion of 1 that are FDR 5% sig in Wen for the same gene 
#  3. Report single percentage overlap for all genes

#--------------------------------------------------------------------------------------

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

message("\n\nCalculating pi1 enrichments: sig. sn-eQTL & public reference eQTL ...")

suppressPackageStartupMessages({
  library(tidyverse)
})

# Input and output paths)
wen_all_qtl <- snakemake@input[['public_all']]
cell_types <- c("ExN-UL", "ExN-DL", "InN", "RG", "MG", "OPC", "Endo-Peri")
public_top_qtl <- snakemake@input[['public_top']]
qtl_all <- snakemake@input[['qtl_all']]
qtl_top <- snakemake@input[['qtl_top']]
output_enrich <- snakemake@output[["enrich"]]
output_pi1 <- snakemake@output[["pi1"]]


library(tidyverse)

cell_types <- c("ExN-UL", "ExN-DL", "InN", "RG", "MG", "OPC", "Endo-Peri")
tensor_dir_pattern <- "../results/05TENSORQTL/tensorqtl_perm/%s_quantile_genPC_4_expPC_40/%s_quantile_perm.cis_qtl.txt.gz"
wen_file <- "../resources/public_datasets/wen_2024/devBrain_eQTL_ALL_nominal_90HCP_FDR_0.05.txt.gz"
lookup <- read_tsv("../resources/sheets/gene_lookup_hg38.tsv") |>
  select(ensembl_gene_id, gene_biotype) |>
  distinct()

# Load Wen et al. FDR < 5% significant SNP-gene pairs -----
message("Loading Wen et al. reference (FDR < 5%) ...")
# For devBrain_eQTL_EUR_nominal_50HCP_FDR_0.05.txt.gz
wen_sig <- read_tsv(wen_file,
  col_types = cols(
    pid = col_character(),
    sid = col_character(),
    dist = col_double(),
    npval = col_double(),
    slope = col_double(),
    fdr = col_double()
  ),
  progress = FALSE
) %>%
  filter(fdr < 0.05) %>%
  transmute(
    phenotype_id = pid,
    variant_id = sid
  ) %>%
  distinct()

# # For eur_perm_eqtl_50HCP_all_assoc.txt.gz
# wen_sig <- read_tsv(wen_file,
#                     col_types   = cols(
#                       pid       = col_character(),
#                       sid       = col_character(),
#                       qval      = col_double(),
#                       .default  = col_skip()
#                     ),
#                     progress = FALSE
# ) %>%
#   filter(qval < 0.05) %>%
#   transmute(
#     phenotype_id = pid,
#     variant_id = sid
#   ) %>%
#   distinct()


# Fast lookup set: "ENSG..._rs..."
wen_sig_set <- paste(wen_sig$phenotype_id, wen_sig$variant_id, sep = "_")
message("Wen FDR < 5% eQTLs loaded: ", nrow(wen_sig))
rm(wen_sig)

# Extract top FDR < 5% eQTL from all 7 cell types  -----
message("Processing TensorQTL results for 7 cell types...")
all_top_eqtls <- map_dfr(cell_types, function(ct) {
  file <- sprintf(tensor_dir_pattern, ct, ct)
  
  message("   → ", ct)
  
  if (!file.exists(file)) {
    warning("File not found: ", file)
    return(tibble())
  }
  
  read_tsv(file,
    col_types = cols(
      phenotype_id = col_character(),
      variant_id    = col_character(),
      qval          = col_double(),
      pval_beta     = col_double(),
      .default      = col_skip()
    ),
    progress = FALSE
  ) %>%
    filter(qval < 0.05) %>%         
    group_by(phenotype_id) %>%
    slice_min(pval_beta, n = 1, with_ties = FALSE) %>%  # Make sure top SNP per gene
    ungroup() %>%
    mutate(cell_type = ct)
})

message("Total top eQTLs (FDR < 5%) across all cell types (genes may repeat): ", nrow(all_top_eqtls))

# Check which of these top eQTLs are also FDR < 5% in Wen  -----
all_top_eqtls <- all_top_eqtls %>%
  mutate(
    gene_snp_key = paste(phenotype_id, variant_id, sep = "_"),
    in_wen_fdr5  = gene_snp_key %in% wen_sig_set
  ) |>
  left_join(lookup, by = join_by(phenotype_id == ensembl_gene_id))

# Final summary  -----
summary_tbl <- all_top_eqtls %>%
  summarise(
    total_eqtls         = n(),
    overlap_wen         = sum(in_wen_fdr5),
    no_overlap_wen      = total_eqtls - overlap_wen,
    proportion_overlap  = overlap_wen / total_eqtls,
    percent_overlap     = scales::percent(proportion_overlap, accuracy = 0.01)
  )

cat("\n=== OVERLAP WITH WEN ET AL. (FDR < 5%) ===\n")
print(summary_tbl)

# Get eQTLs by gene biotype -----
replicated_by_biotype <- all_top_eqtls %>%
  filter(in_wen_fdr5) %>%                     
  count(gene_biotype, name = "n_replicated") %>%
  arrange(desc(n_replicated)) %>%          
  mutate(
    percent_of_replicated = n_replicated / sum(n_replicated),
    percent_of_replicated = scales::percent(percent_of_replicated, accuracy = 0.1)
  )

cat("\n=== Replicated top eQTLs (also in Wen FDR<5%) by gene biotype ===\n")
print(replicated_by_biotype, n = Inf)

# Get eQTLs by gene biotype not in Wen  -----
non_replicated_by_biotype <- all_top_eqtls %>%
  filter(!in_wen_fdr5) %>%                     
  count(gene_biotype, name = "n_non_replicated") %>%
  arrange(desc(n_non_replicated)) %>%          
  mutate(
    percent_of_non_replicated = n_non_replicated / sum(n_non_replicated),
    percent_of_non_replicated = scales::percent(percent_of_non_replicated, accuracy = 0.1)
  )

cat("\n=== Non-replicated top eQTLs (not in Wen FDR<5%) by gene biotype ===\n")
print(non_replicated_by_biotype, n = Inf)

# Save results if you want
# write_tsv(all_top_eqtls, "L1_top_eQTLs_with_Wen_overlap_status.tsv")
# write_tsv(summary_tbl,   "L1_Wen_overlap_summary.tsv")
# 
# # Return everything
# list(
#   detailed = all_top_eqtls,
#   summary  = summary_tbl
# )

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------