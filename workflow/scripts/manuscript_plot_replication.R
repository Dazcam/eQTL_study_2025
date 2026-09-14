#--------------------------------------------------------------------------------------
#
#    Generate eQTL p1 replication plt for manuscript
#
#--------------------------------------------------------------------------------------

# A: Pie chart of shared eGenes accross cell types
# B: Upset Plot
# C: Internal Pi1 heatmap
# D: Fetal vs. adult Pi1 heatmap (Jang et al. 2026; Glu-A, Glu-B, GABA, NPC only)
# E: Fetal vs. adult beta correlation - Glu-A vs Ext (Jang)
# F: Fetal vs. adult beta correlation - Glu-B vs Ext (Jang)
# G: Fetal vs. adult beta correlation - GABA vs IN (Jang)

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

message("\n\nGenerating replication plot for the manuscript ...")

# -------------------------------------------------------------------------------------
library(tidyverse)
library(cowplot)
library(ggrepel)
library(UpSetR)
library(grid)

# --- Set variables
in_dir <- snakemake@params[['in_dir']]
internal_dir <- snakemake@params[['internal_dir']]
jang_pi1_dir <- snakemake@params[['jang_pi1_dir']]
jang_beta_file <- snakemake@input[['beta_file']]
out_file <- snakemake@output[[1]]

# in_dir <- '../results/05TENSORQTL/tensorqtl_perm/'
# internal_dir <- "../results/06QTL-REPLICATION/internal/"
# jang_pi1_dir <- "../results/19DEV-SPECIFICITY/pi1_jang/"
# jang_beta_file <- "../results/19DEV-SPECIFICITY/<step5_output>.rds"
# out_dir <- "../results/13MANUSCRIPT_PLOTS_TABLES/"

# NOTE: these old-format names (Glu-UL / Glu-DL) drive file paths for the eQTL
# input files (Panel A) and the internal Pi1 comparison (Panel C), which are
# still stored under the old naming on disk. Do not rename these - the
# relabel_cell_type() helper below handles renaming for display only.
cell_types <- c("Glu-UL", "Glu-DL", "NPC", "GABA", "Endo-Peri", "OPC", "MG")

# Exp PC map (old names - keys file paths for Panels A and C)
expPC_map <- c(
  "Glu-UL"     = 50,
  "Glu-DL"     = 40,
  "GABA"       = 30,
  "NPC"        = 30,
  "MG"         = 30,
  "OPC"        = 30,
  "Endo-Peri"  = 30
)

# --- Relabel map for cell type name changes, for display/plotting only
relabel_cell_type <- function(x) {
  x <- str_replace(x, "Glu-UL", "Glu-A")
  x <- str_replace(x, "Glu-DL", "Glu-B")
  x
}

cell_types_new <- relabel_cell_type(cell_types)

# Define custom color palette (new Glu-A/Glu-B keys, since this is applied
# only to already-relabelled data downstream)
custom_palette <- c(
  'Glu-A' = '#4363d8',
  'Glu-B' = '#00B6EB',
  'NPC' = '#FF5959',
  'GABA' = '#3CBB75FF',
  'Endo-Peri' = '#B200ED',
  'MG' = '#F58231',
  'OPC' = '#FDE725FF'
)

eqtl_list <- list()
genPC <- 4
norm_method <- 'quantile'

gene_lookup <- read_tsv('../resources/sheets/gene_lookup_hg38.tsv') |>
  select(gene = ensembl_gene_id, symbol = external_gene_name)


# ----- Pie chart shared eQTL amoung L1 cell types (Panel A - unchanged) -----
for (cell_type in cell_types) {
  
  expPC <- expPC_map[[cell_type]]
  
  f <- file.path(
    in_dir,
    paste0(cell_type, "_", norm_method, "_genPC_", genPC, "_expPC_", expPC),
    paste0(cell_type, "_", norm_method, "_perm.cis_qtl.txt.gz")  
  )
  
  if (file.exists(f)) {
    dat <- read_tsv(f, show_col_types = FALSE) %>%
      filter(!is.na(start_distance)) %>%
      mutate(distance_kb = start_distance / 1000) %>%
      filter(qval < 0.05)   # significance cutoff
    
    if (nrow(dat) > 0) {
      eqtl_list[[cell_type]] <- dat %>% mutate(cell_type = cell_type)
    }
  }
}

eqtls_all <- bind_rows(eqtl_list)

# Relabel here (feeds both Panel A's counts, which are cell-type-agnostic,
# and Panel B's set names, which do need the new labels)
gene_cell <- eqtls_all %>%
  distinct(phenotype_id, cell_type) %>%
  mutate(cell_type = relabel_cell_type(cell_type))

gene_counts <- gene_cell %>%
  count(phenotype_id) %>%  
  count(n)  

label_threshold <- 0.08

pie_dat <- gene_counts %>%
  mutate(label = paste0(n, " cell type", ifelse(n > 1, "s", "")),
         count = nn) %>%
  arrange(n) %>%
  mutate(fraction = count / sum(count),
         ymax = cumsum(fraction),
         ymin = lag(ymax, default = 0),
         midpoint = (ymax + ymin) / 2,
         label = factor(label, levels = label)) |>
  mutate(large_slice = fraction >= label_threshold,
         label_x = ifelse(large_slice, 0.5, 1.15))   

mono_cols <- colorRampPalette(c("#08306B", "#DEEBF7"))(nrow(pie_dat))

connectors <- subset(pie_dat, !large_slice)

connectors <- transform(
  connectors,
  x_start = 1,
  x_end   = 1.12,
  y_start = midpoint,
  y_end   = midpoint
)

pie_chart <- ggplot(pie_dat,aes(ymax = ymax, ymin = ymin, xmax = 1, 
                   xmin = 0, fill = label)) +
  geom_rect(color = "black", linewidth = 0.2) +
  coord_polar(theta = "y") +
  geom_text(data = subset(pie_dat, large_slice),
            aes(x = 0.5,y = midpoint, label = count),
            color = "white", size = 5, fontface = "bold") +
  geom_text_repel(data = subset(pie_dat, !large_slice),
                  aes(x = 1.15, y = midpoint, label = count),
                  size = 5,
                  segment.color = "grey60",
                  show.legend = FALSE,
                  fontface = "bold") +
  scale_fill_manual(values = mono_cols) +
  theme_void(base_size = 15) +
  theme(legend.position = "bottom", 
        legend.title = element_blank(),
        legend.box.margin = margin(t = -50, r = 0, b = 10, l = 0),
        plot.margin = margin(t = 30, r = 30, b = 30, l = 30)) +
  geom_segment(
    data = connectors,
    aes(x = x_start,
        xend = x_end,
        y = y_start,
        yend = y_end
    ),
    inherit.aes = FALSE,
    linewidth = 0.4,
    color = "black"
  )

# --- Upset Plot (Panel B - relabelled set names) ----
gene_by_cell <- gene_cell %>%
  pivot_wider(names_from = cell_type, values_from = cell_type,
              values_fill = 0, values_fn = function(x) 1) %>%
  column_to_rownames("phenotype_id")

# Render to PNG 
tmp_upset <- tempfile(fileext = ".png")

png(tmp_upset, width = 2800, height = 1800, res = 300)

upset(
  gene_by_cell,
  nsets          = length(cell_types_new),
  order.by       = "freq",
  nintersects    = 20,
  sets.bar.color = custom_palette[cell_types_new],
  point.size     = 3.8,
  line.size      = 2,
  # text.scale: c(intersection size title, intersection size tick labels, 
  #               set size title, set size tick labels, 
  #               set names, numbers above bars)
  text.scale     = c(2.2, 2, 1.8, 1.8, 1.8, 1.7)
)

dev.off()

# Re-import as grob for cowplot
img <- png::readPNG(tmp_upset)

upset_grob <- grid::rasterGrob(
  img,
  interpolate = TRUE
)

upset_plt <- ggplotify::as.ggplot(upset_grob) +
  theme(
    plot.margin = margin(t = 40, r = 0, b = 40, l = 0, unit = "pt")
  )

# --- Internal pi1 heatmap (Panel C - relabelled + new axis titles) -----
read_pi1_results <- function(ct, ref_ct) {
  
  expPC <- expPC_map[[ct]]
  ref_expPC <- expPC_map[[ref_ct]]
  
  file_path <- paste0(internal_dir, ct, "_vs_", ref_ct, "_quantile_genPC_4_expPC_", 
                      expPC, "_expPCref_", ref_expPC,"_pi1_results_tbl.rds")
  if (file.exists(file_path)) {
    pi1_list <- read_rds(file_path)
    
    tibble(
      cell_type = ct,
      ref_cell_type = ref_ct,
      pi1_forward = pi1_list$forward$pi1,
      pi1_reverse = pi1_list$reverse$pi1,
      prop_replicating_forward = pi1_list$forward$prop_replicating,
      prop_replicating_reverse = pi1_list$reverse$prop_replicating,
      prop_same_direction_forward = pi1_list$forward$prop_same_direction,
      prop_same_direction_reverse = pi1_list$reverse$prop_same_direction
    )
  } else {
    message("File missing: ", file_path)
    NULL
  }
}

# Generate all combinations (old names - matches file paths on disk)
combinations <- expand.grid(cell_type = cell_types, ref_cell_type = cell_types, stringsAsFactors = FALSE)

# Read all pi1 results
pi1_result_tbl <- map2_dfr(combinations$cell_type, combinations$ref_cell_type, read_pi1_results)

# Prepare the full square data without manual mirroring
pi1_square_tbl <- pi1_result_tbl %>%
  mutate(
    # Use the specific Forward result for every unique combination provided by expand.grid
    pi1_final = ifelse(cell_type == ref_cell_type, 1.0, pi1_forward)
  ) %>%
  select(ref = cell_type, query = ref_cell_type, pi1 = pi1_final) %>%
  mutate(
    ref = relabel_cell_type(ref),
    query = relabel_cell_type(query),
    # rev() ensures the first cell type is at the top (standard matrix view)
    query = factor(query, levels = rev(cell_types_new)), 
    ref = factor(ref, levels = cell_types_new)
  )

# Updated heatmap function
plot_int_heatmap <- function(df) {
  ggplot(df, aes(x = query, y = ref, fill = pi1)) +
    geom_tile(color = "black", lwd = 1.1, linetype = 1) +
    geom_text(aes(label = ifelse(is.na(pi1), "NA", sprintf("%.2f", pi1))),
              color = "black", size = 3.5) +
    scale_fill_gradientn(
      colours = c("white", "yellow", "red"),
      limits = c(0.5, 1.0),
      na.value = "grey80",
      name = expression(pi[1])
    ) +
    coord_equal() +
    labs(x = "Prenatal cell type (replication)", y = "Prenatal cell type (discovery)") +
    theme_minimal(base_size = 13) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
      axis.text.y = element_text(face = "bold"),
      panel.grid = element_blank()
    )
}

# Generate only the forward heatmap (no title)
pi1_int_heatmap <- plot_int_heatmap(pi1_square_tbl)


# --- Jang pi1 heatmap (Panel D - replaces Fugita; restricted to Glu-A, Glu-B,
# GABA, NPC rows in that order, all 7 Jang adult reference types as columns) -----
jang_L1_cell_types <- c("Glu-UL", "Glu-DL", "GABA", "NPC")  # old names - match file names on disk
jang_adult_cell_types <- c("Ast", "End", "Ext", "IN", "MG", "OD", "OPC")

read_jang_pi1_results <- function(ct, ref_ct, jang_dir) {
  file_path <- paste0(jang_dir, ct, "_vs_", ref_ct, "_pi1_results_tbl.rds")
  if (file.exists(file_path)) {
    pi1_list <- read_rds(file_path)
    
    tibble(
      cell_type = ct,
      ref_cell_type = ref_ct,
      pi1_forward = pi1_list$forward$pi1,
      pi1_reverse = pi1_list$reverse$pi1
    )
  } else {
    message("File missing: ", file_path)
    NULL
  }
}

jang_combinations <- expand.grid(cell_type = jang_L1_cell_types, ref_cell_type = jang_adult_cell_types,
                                 stringsAsFactors = FALSE)

pi1_jang_tbl <- map2_dfr(jang_combinations$cell_type, jang_combinations$ref_cell_type,
                         read_jang_pi1_results, jang_dir = jang_pi1_dir)

pi1_jang_square_tbl <- pi1_jang_tbl %>%
  select(cell_type, ref_cell_type, pi1 = pi1_forward) %>%
  mutate(
    cell_type = relabel_cell_type(cell_type),
    # rev() so Glu-A appears at the top (standard matrix view)
    cell_type = factor(cell_type, levels = rev(relabel_cell_type(jang_L1_cell_types))),
    ref_cell_type = factor(ref_cell_type, levels = jang_adult_cell_types)
  )

# Function to generate the Jang heatmap
plot_jang_heatmap <- function(df) {
  ggplot(df, aes(x = ref_cell_type, y = cell_type, fill = pi1)) +
    geom_tile(color = "black", lwd = 1.1, linetype = 1) +
    geom_text(aes(label = ifelse(is.na(pi1), "NA", sprintf("%.2f", pi1))),
              color = "black", size = 3.5) +
    scale_fill_gradientn(
      colours = c("white", "yellow", "red"),
      limits = c(0, 1.0),
      na.value = "grey80",
      name = expression(pi[1])
    ) +
    coord_equal() +
    labs(x = "Adult cell type (replication)", y = "Prenatal cell type (discovery)") +
    theme_minimal(base_size = 13) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
      axis.text.y = element_text(face = "bold"),
      panel.grid = element_blank(),
      plot.margin = margin(50, 50, 50, 50, unit = "pt")
    ) 
}

pi1_jang_heatmap <- plot_jang_heatmap(pi1_jang_square_tbl)


### --- beta correlation plt (Panels E-G - Jang data) -----
jang_beta_list <- read_rds(jang_beta_file)
paired_betas_all <- jang_beta_list$paired_betas

# genes dropped from the table entirely (not just left unlabelled) for E/F/G
exclude_genes <- c("ABCC8", "CLHC1")

make_jang_beta_plot <- function(paired_betas_all, my_ct, jang_ct, gene_lookup,
                                exclude_genes = character(0), label_genes = NULL,
                                ct_labels) {
  
  paired_betas <- paired_betas_all %>%
    filter(my_cell_type == my_ct, jang_cell_type == jang_ct) %>%
    left_join(gene_lookup, by = "gene") %>%
    mutate(label_text = ifelse(is.na(symbol) | symbol == "NA", gene, symbol)) %>%
    filter(!(label_text %in% exclude_genes))
  
  ## Correlation
  cor_val <- cor(
    paired_betas$beta_my,
    paired_betas$beta_jang,
    use = "complete.obs",
    method = "pearson"
  )
  
  cor_label <- sprintf("r = %.2f", cor_val)
  
  # --- Point classification: two-check "significant_jang" convention (see
  # devspec_beta_correlation.R / devspec_report.Rmd), replacing the old
  # sign-mismatch-only "discordant" rule. A pair is only highlighted as
  # concordant/discordant if it is ALSO significant in Jang (i.e. it IS that
  # gene's lead SNP in Jang's top-association file AND qval < 0.05) -
  # everything else, including sign-mismatched pairs that aren't
  # Jang-significant, falls into a single grey "other" bucket.
  paired_betas <- paired_betas |>
    mutate(
      concordant = sign(beta_my) == sign(beta_jang),
      point_category = case_when(
        significant_jang & concordant  ~ "concordant_both_sig",
        significant_jang & !concordant ~ "discordant_both_sig",
        TRUE                            ~ "other"
      ),
      should_label = if (!is.null(label_genes)) label_text %in% label_genes else FALSE
    )
  
  ggplot(paired_betas, aes(x = beta_jang, y = beta_my)) +
    geom_point(
      data = filter(paired_betas, point_category == "other"),
      alpha = 0.2, color = "grey70", size = 1.8
    ) +
    geom_point(
      data = filter(paired_betas, point_category == "discordant_both_sig"),
      color = "#d32f2f", size = 3, alpha = 0.6
    ) +
    geom_point(
      data = filter(paired_betas, point_category == "concordant_both_sig"),
      color = "#1976d2", size = 3, alpha = 0.6
    ) +
    geom_smooth(method = "lm", color = "red", se = TRUE, linewidth = 0.9) +
    geom_abline(
      slope = 1, intercept = 0, linetype = "dashed",
      color = "grey50", linewidth = 0.5
    ) +
    geom_text_repel(
      data = filter(paired_betas, should_label),
      aes(label = label_text),
      size = 3.1,
      fontface = "bold",
      box.padding = 0.45,
      point.padding = 0.5,
      segment.color = "grey50",
      segment.size = 0.25,
      min.segment.length = 0,
      max.overlaps = 25,
      force = 1.5,
      force_pull = 0.8,
      direction = "both",
      seed = 2025
    ) +
    annotate(
      "text", x = -1.25, y = Inf, label = cor_label,
      hjust = 0, vjust = 1.8, size = 5
    ) +
    labs(
      x = substitute("Adult" ~ x ~ beta, list(x = ct_labels[2])),
      y = substitute("Prenatal" ~ x ~ beta, list(x = ct_labels[1]))
    ) +
    coord_cartesian(clip = "off") +
    theme_minimal(base_size = 13) +
    theme(
      plot.margin = margin(25, 50, 22, 50, unit = "pt"),
      axis.title.x = element_text(margin = margin(t = 10)),
      axis.title.y = element_text(margin = margin(r = 10)),
      title = element_blank(),
      subtitle = element_blank(),
    )
}

beta_gluA_plt <- make_jang_beta_plot(paired_betas_all, "Glu-UL", "Ext", gene_lookup,
                                     exclude_genes = exclude_genes,
                                     label_genes = "PLBD2",
                                     ct_labels = c("Glu-A", "Ext"))

beta_gluB_plt <- make_jang_beta_plot(paired_betas_all, "Glu-DL", "Ext", gene_lookup,
                                     exclude_genes = exclude_genes,
                                     label_genes = c("PLBD2", "PTPA"),
                                     ct_labels = c("Glu-B", "Ext"))

beta_gaba_plt <- make_jang_beta_plot(paired_betas_all, "GABA", "IN", gene_lookup,
                                     exclude_genes = exclude_genes,
                                     label_genes = "PLBD2",
                                     ct_labels = c("GABA", "IN"))

### --- plot -----
# Final plot
top_row <- plot_grid(pie_chart, upset_plt, labels = c("A", "B"),
  label_size = 24, ncol = 2,rel_widths = c(1, 1.3))

# Stack heatmaps
heatmaps_stacked <- plot_grid(
  pi1_int_heatmap,
  pi1_jang_heatmap,
  ncol = 1,
  rel_heights = c(1, 0.8),   # C slightly taller
  labels = c("C", "D"),
  label_size = 24,
  align = "v"
)

# Stack betas 
betas_stacked <- plot_grid(
  beta_gluA_plt, beta_gluB_plt, beta_gaba_plt,
  ncol = 1,
  rel_heights = c(1,1,1),
  labels = c("E", "F", "G"),
  label_size = 24
)

# Bottom row: heatmaps + betas 
bottom_row <- plot_grid(
  heatmaps_stacked,
  betas_stacked,
  ncol = 2,
  rel_widths = c(1.1, 0.9), 
  align = "h"
)
final_plt <- plot_grid(
  top_row,
  bottom_row,
  ncol = 1,
  rel_heights = c(1, 1.4)     
)

ggsave(
  filename = out_file,
  plot = final_plt,
  width = 16,
  height = 14,
  units = "in",
  device = "pdf",
  dpi = 300
)

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
