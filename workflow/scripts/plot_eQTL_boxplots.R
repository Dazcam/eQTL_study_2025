#--------------------------------------------------------------------------------------
#
#    Plot eQTL boxplots
#
#--------------------------------------------------------------------------------------

## Load libraries and variables -------------------------------------------------------
library(tidyverse)
library(cowplot)

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

message('\n\nCreating SLDSR annotation files from susie output ...')


exp_dir <- "../results/05TENSORQTL/prep_input/"
pval_dir <- "../results/10SMR/smr_input/"
gen_prefix <- "../results/05TENSORQTL/prep_input/chrALL_final.filt"
gene_id <- "ENSG00000214435"
snp_id <- "rs11191424"  # Can be rsID; code handles rsID lookup via pvar

gene_id <- snakemake@input[['exp_dir']]
pval_dir <- snakemake@input[['pval_dir']]
gen_prefix <- snakemake@input[['gen_prefix']]
gene_id <- snakemake@input[['gene_id']]
snp_id <- snakemake@input[['snp_id']]
output <- snakemake@output[['output']]


cell_types <- c("ExN-UL", "ExN-DL", "InN", "RG", "MG", "OPC", "Endo-Peri")

cat("============================")
tibble(
  variable = c("gene_id", "pval_dir", "gen_prefix", "gene_id", "snp_id", "output"),
  value    = c(gene_id, pval_dir, gen_prefix, gene_id, snp_id, output)) |> 
  knitr::kable(format = "simple", align = "l") |>
  print()
message("\n============================\n")

# User inputs
# cell_types <- c("ExN-UL", "ExN-DL", "InN", "RG", "MG", "OPC", "Endo-Peri")
# exp_dir <- "../results/05TENSORQTL/prep_input/"
# pval_dir <- "../results/10SMR/smr_input/"
# gen_prefix <- "../results/05TENSORQTL/prep_input/chrALL_final.filt"
# gene_id <- "ENSG00000214435"
# snp_id <- "rs11191424"  # Can be rsID; code handles rsID lookup via pvar

# Load pvar to find SNP details (assumes standard PLINK2 pvar format)
pvar_file <- paste0(gen_prefix, ".pvar")
pvar <- read_tsv(pvar_file, comment = "", col_names = TRUE, col_types = cols(
  `#CHROM` = col_character(),
  POS = col_integer(),
  ID = col_character(),
  REF = col_character(),
  ALT = col_character()
))
colnames(pvar)[1] <- "CHROM"  # Clean up the #CHROM column name

# Find SNP row
snp_row <- pvar %>% filter(ID == snp_id)
if (nrow(snp_row) == 0) {
  stop(paste("SNP", snp_id, "not found in pvar file."))
}
chrom <- snp_row$CHROM[1]
pos <- snp_row$POS[1]
ref <- snp_row$REF[1]
alt <- snp_row$ALT[1]
variant_id <- paste0(chrom, "_", pos, "_", ref, "_", alt)

# Extract genotype dosages using PLINK2 (assumes plink2 is in PATH)
temp_rs <- tempfile(fileext = ".txt")
write_lines(snp_id, temp_rs)
temp_out <- tempfile()
system2("plink2", args = c(
  "--pfile", gen_prefix,
  "--extract", temp_rs,
  "--export", "A-transpose",
  "--out", temp_out
))
dosage_file <- paste0(temp_out, ".A-transpose")

# Read the exported dosages (no header: FID IID dosage)
geno_df <- read_tsv(dosage_file, col_names = c("FID", "IID", "genotype_dosage")) %>%
  select(IID, genotype_dosage) %>%
  rename(sample_id = IID)

# Load p-values for this gene-SNP pair across cell types
pval_map <- tibble(cell_type = character(), pval_nominal = double())
for (ct in cell_types) {
  pval_file <- paste0(pval_dir, ct, "/", ct, "_nom.cis_qtl_pairs.tsv")
  if (file.exists(pval_file)) {
    dt <- read_tsv(pval_file, col_types = cols(
      gene_id = col_character(),
      variant_id = col_character(),
      pval_nominal = col_double()
    ))
    hit <- dt %>% filter(gene_id == !!gene_id & variant_id == !!variant_id)
    if (nrow(hit) > 0) {
      pval_map <- bind_rows(pval_map, tibble(
        cell_type = ct,
        pval_nominal = hit$pval_nominal[1]  # Take first if multiple
      ))
    }
  }
}
pval_map <- pval_map %>% mutate(pval_str = sprintf("P=%.2g", pval_nominal))

# Generate plots per cell type
plots <- list()
for (ct in cell_types) {
  expr_file <- paste0(exp_dir, ct, "_quantile.bed")
  if (!file.exists(expr_file)) next
  
  # Read expression BED file (header starts with #Chr)
  expr_bed <- read_tsv(expr_file, comment = "", col_names = TRUE)
  colnames(expr_bed)[1] <- "Chr"  # Clean up #Chr column name
  
  # Find gene row
  gene_idx <- which(expr_bed$TargetID == gene_id)
  if (length(gene_idx) == 0) next
  gene_row <- expr_bed[gene_idx, ]
  
  # Extract sample names and expression values (columns 5+ are samples)
  sample_names <- colnames(gene_row)[5:ncol(gene_row)]
  expression <- as.numeric(gene_row[1, 5:ncol(gene_row)])
  
  # Build plot data for this cell type
  plot_data_ct <- tibble(
    sample_id = sample_names,
    cell_type = ct,
    expression = expression
  ) %>%
    filter(!is.na(expression)) %>%
    left_join(geno_df, by = "sample_id") %>%
    filter(!is.na(genotype_dosage)) %>%
    mutate(
      genotype = case_when(
        genotype_dosage == 0 ~ paste0(ref, "/", ref),
        genotype_dosage == 1 ~ paste0(ref, "/", alt),
        genotype_dosage == 2 ~ paste0(alt, "/", alt),
        TRUE ~ NA_character_
      ),
      genotype = factor(genotype, levels = c(paste0(ref, "/", ref), paste0(ref, "/", alt), paste0(alt, "/", alt)))
    )
  
  if (nrow(plot_data_ct) == 0) next
  
  # Get p-value for this cell type
  tmp_pval <- pval_map %>% filter(cell_type == ct) %>% pull(pval_str)
  cell_pval <- ifelse(length(tmp_pval) == 0, "P=NA", tmp_pval)
  
  # Compute x-axis labels with sample counts
  present_counts <- plot_data_ct %>% count(genotype)
  geno_levels <- levels(plot_data_ct$genotype)
  x_labels <- sapply(geno_levels, function(g) {
    n <- present_counts$n[present_counts$genotype == g]
    if (length(n) == 0) n <- 0
    paste0(as.character(g), " (", n, ")")
  })
  
  # Create plot
  p_ct <- ggplot(plot_data_ct, aes(x = genotype, y = expression)) +
    geom_boxplot(fill = "lightblue", alpha = 0.7, outlier.shape = NA) +
    geom_jitter(width = 0.2, alpha = 0.6, size = 0.8) +
    scale_x_discrete(labels = x_labels) +
    labs(
      title = ct,
      x = "",
      y = "Normalised Expression",
      subtitle = cell_pval
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 10, face = "bold"),
      plot.subtitle = element_text(size = 9, color = "red"),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8)
    )
  
  plots[[ct]] <- p_ct
}

# Combine plots
if (length(plots) > 0) {
  combined_plot <- plot_grid(plotlist = plots, ncol = 4)
  saveRDS(combined_plot, output)
  
  # Optionally save
  # ggsave("eqtl_plots.png", combined_plot, width = 12, height = 8, dpi = 300)
} else {
  cat("No data available for plotting.\n")
}