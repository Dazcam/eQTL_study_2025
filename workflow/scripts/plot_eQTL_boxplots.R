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

cell_types <- c("ExN-UL", "ExN-DL", "InN", "RG", "MG", "OPC", "Endo-Peri")

cat("============================")
tibble(
  variable = c("exp_dir", "pval_dir", "geno_prefix", "gene_id", "snp_id", "output"),
  value    = c(exp_dir, pval_dir, geno_prefix, gene_id, snp_id, output)) |> 
  knitr::kable(format = "simple", align = "l") |>
  print()
message("\n============================\n")

# PLINK2 full path (from your which output; adjust if needed)
plink2_path <- "/apps/genomics/plink/2.0/el7/AVX512/intel-2018/serial/plink-2.0/build_dynamic/plink2"

# Debug: Check if plink2 path exists inside container
if (!file.exists(plink2_path)) {
  stop(paste("plink2 not found at:", plink2_path, "\nThis suggests the /apps path is not mounted in the singularity container. Check snakemake singularity bind paths or extract genotypes outside the script."))
}

# Set output directory for non-temp files
output_dir <- dirname(output)

# Load pvar to find SNP details (handle VCF-like header with multiple ## lines)
message('Loading pvar file to get SNP info ...')
pvar_file <- paste0(geno_prefix, ".pvar")
pvar_lines <- read_lines(pvar_file)
header_idx <- which(str_starts(pvar_lines, "#CHROM"))
if (length(header_idx) == 0) {
  stop("Could not find #CHROM header line in pvar file.")
}
skip_lines <- header_idx - 1
pvar <- read_tsv(pvar_file, skip = skip_lines, col_names = TRUE, col_types = cols(
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

# Load psam for sample IDs
psam_file <- paste0(geno_prefix, ".psam")
psam <- read_tsv(psam_file, col_names = TRUE, show_col_types = FALSE)

# Extract genotype dosages using PLINK2 (write extract file to output_dir)
extract_file <- file.path(output_dir, paste0(snp_id, "_extract.txt"))
write_lines(snp_id, extract_file)
out_prefix <- file.path(output_dir, "plink_out")
cmd_args <- c(
  "--pfile", geno_prefix,
  "--extract", extract_file,
  "--export", "A-transpose",
  "--out", out_prefix
)
cat("Running PLINK2 with args:", paste(cmd_args, collapse = " "), "\n")
plink_out <- system2(plink2_path, args = cmd_args, stdout = TRUE, stderr = TRUE)
if (length(plink_out) > 0) {
  cat("PLINK2 output:\n")
  cat(plink_out, sep = "\n")
}
dosage_file <- paste0(out_prefix, ".traw")
if (!file.exists(dosage_file)) {
  stop(paste("PLINK2 export failed. Dosage file not created:", dosage_file, "\nCheck PLINK2 output above for errors (e.g., SNP not found, file paths)."))
}

# Read the .traw file (header present, one row for the single variant)
traw <- read_tsv(dosage_file, col_names = TRUE, show_col_types = FALSE)
if (nrow(traw) != 1) {
  stop("Unexpected number of rows in .traw file (expected 1 for single SNP).")
}
dosages <- as.numeric(traw[1, 7:ncol(traw)])
geno_df <- tibble(
  sample_id = psam$IID,
  genotype_dosage = dosages
) %>%
  filter(!is.na(genotype_dosage))

message('Read exported dosages ...')

# Load p-values for this gene-SNP pair across cell types
message('Loading p-values for SNP / Gene pair for TensorQTL nominal file ...') 
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
  message('Generating plot for ', ct, '...') 
  expr_file <- paste0(exp_dir, ct, "_quantile.bed")
  if (!file.exists(expr_file)) next
  
  # Read expression BED file (handle #Chr header)
  expr_bed <- read_tsv(expr_file, comment = "", col_names = TRUE)
  colnames(expr_bed)[1] <- "Chr"  # Clean up #Chr column name
  
  # Find gene row
  gene_idx <- which(expr_bed$TargetID == gene_id)
  if (length(gene_idx) == 0) next
  
  # Extract sample names and expression values (columns 5+ are samples)
  sample_names <- colnames(expr_bed)[5:ncol(expr_bed)]
  expression <- as.numeric(expr_bed[gene_idx, 5:ncol(expr_bed)])
  
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
  message('Saving combined plot ...')
  combined_plot <- plot_grid(plotlist = plots, ncol = 4)
  saveRDS(combined_plot, output)
  
  # Optionally save
  # ggsave("eqtl_plots.png", combined_plot, width = 12, height = 8, dpi = 300)
} else {
  cat("No data available for plotting.\n")
}
message('All done.')