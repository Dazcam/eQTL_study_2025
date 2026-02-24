#--------------------------------------------------------------------------------------
#
#    Generate eQTL p1 replication plt for manuscript
#
#--------------------------------------------------------------------------------------

# A: Pie chart of shared eGenes accross cell types
# B: Upset Plot
# C: Internal Pi1 heatmap
# D: Fetal vs. adult Pi1 heatmap (only Glu and GABA only)
# E: Fetal vs. adult beta correlation 1
# F: Fetal vs. adult beta correlation 2
# G: Fetal vs. adult beta correlation 3

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
fugita_dir <- snakemake@params[['fugita_dir']] 
beta_dir <- snakemake@params[['beta_dir']] 
out_file <- snakemake@output[[1]]

# in_dir <- '../results/05TENSORQTL/tensorqtl_perm/'
# internal_dir <- "../results/06QTL-REPLICATION/internal/"
# fugita_dir <- "../results/06QTL-REPLICATION/fugita/"
# out_dir <- "../results/13MANUSCRIPT_PLOTS_TABLES/"
# beta_dir <- "../results/06QTL-REPLICATION/beta_cor/"

cell_types <- c("Glu-UL", "Glu-DL", "NPC", "GABA", "Endo-Peri", "OPC", "MG")
cell_types_sub <- c("Glu-UL", "Glu-DL", "NPC", "GABA")
fugita_cell_types <- c('Ast', 'Exc', 'Inh', 'Oli', 'OPC', 'Mic', 'End')

# Exp PC map
expPC_map <- c(
  "Glu-UL"     = 50,
  "Glu-DL"     = 40,
  "GABA"       = 30,
  "NPC"        = 30,
  "MG"         = 30,
  "OPC"        = 30,
  "Endo-Peri"  = 30
)

# Define custom color palette
custom_palette <- c(
  'Glu-UL' = '#4363d8',
  'Glu-DL' = '#00B6EB',
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


# ----- Pie chart shared eQTL amoung L1 cell types -----
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

gene_cell <- eqtls_all %>%
  distinct(phenotype_id, cell_type)

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

# --- Upset Plot ----
gene_by_cell <- gene_cell %>%
  pivot_wider(names_from = cell_type, values_from = cell_type,
              values_fill = 0, values_fn = function(x) 1) %>%
  column_to_rownames("phenotype_id")

# ---- Derive intersection ordering EXACTLY as upset() does ----
membership <- UpSetR:::Create_matrix(gene_by_cell)
intersects <- UpSetR:::Get_intersections(membership)

# apply same cutoff + ordering
intersects <- intersects[intersects$freq >= 10, ]
intersects <- intersects[order(-intersects$freq), ]

# ---- Identify unique intersections ----
intersection_cell <- apply(
  intersects[, cell_types],
  1,
  function(x) {
    if (sum(x) == 1) {
      names(which(x == 1))
    } else {
      NA
    }
  }
)

# ---- Colour mapping ----
bar_cols <- ifelse(
  is.na(intersection_cell),
  "grey40",
  custom_palette[intersection_cell]
)

# ---- Render to PNG first (device-safe for Snakemake/pdf()) ----
tmp_upset <- tempfile(fileext = ".png")

png(tmp_upset, width = 2400, height = 2000, res = 300)

upset(
  gene_by_cell,
  nsets          = length(cell_types),
  order.by       = "freq",
  cutoff         = 10,
  main.bar.color = bar_cols,
  matrix.color   = bar_cols,
  sets.bar.color = custom_palette[cell_types],
  point.size     = 3.8,
  line.size      = 1
)

dev.off()

# ---- Re-import as grob for cowplot ----
img <- png::readPNG(tmp_upset)

upset_grob <- grid::rasterGrob(
  img,
  interpolate = TRUE
)

upset_plt <- ggplotify::as.ggplot(upset_grob)

# --- Internal pi1 heatmap -----
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

# Generate all combinations
combinations <- expand.grid(cell_type = cell_types, ref_cell_type = cell_types, stringsAsFactors = FALSE)

# Read all pi1 results
pi1_result_tbl <- map2_dfr(combinations$cell_type, combinations$ref_cell_type, read_pi1_results)

# Prepare data for Forward only (upper triangle including diagonal)
pi1_forward_tbl <- pi1_result_tbl %>%
  select(query = cell_type, ref = ref_cell_type, pi1 = pi1_forward) %>%
  mutate(
    pi1 = ifelse(query == ref, 1.0, pi1),  # diagonal = 1.0
    row_idx = match(query, cell_types),
    col_idx = match(ref, cell_types)
  ) %>%
  filter(row_idx <= col_idx)  # upper triangle including diagonal (opposite direction)

# Force factor levels for correct ordering
pi1_int_tbl <- pi1_forward_tbl %>%
  mutate(
    query = factor(query, levels = cell_types),
    ref = factor(ref, levels = cell_types)
  )

# Updated heatmap function
plot_int_heatmap <- function(df) {
  ggplot(df, aes(x = ref, y = query, fill = pi1)) +
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
    labs(x = "Cell type", y = "Cell type") +
    theme_minimal(base_size = 13) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
      axis.text.y = element_text(face = "bold"),
      panel.grid = element_blank()
    )
}

# Generate only the forward heatmap (no title)
pi1_int_heatmap <- plot_int_heatmap(pi1_int_tbl)


# --- Fugita pi1 heatmap -----
# Function to read pi1 results
read_pi1_results <- function(ct, ref_ct) {
  
  expPC <- expPC_map[[ct]]
  
  file_path <- paste0(fugita_dir, ct, "_vs_", ref_ct, "_quantile_genPC_4_expPC_", expPC ,"_pi1_results_tbl.rds")
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

# Generate all combinations
combinations <- expand.grid(cell_type = cell_types_sub, ref_cell_type = fugita_cell_types, stringsAsFactors = FALSE)

# Read all pi1 results
pi1_result_tbl <- map2_dfr(combinations$cell_type, combinations$ref_cell_type, read_pi1_results)

# Pivot for plotting
pi1_long <- pi1_result_tbl %>%
  select(cell_type, ref_cell_type, pi1_forward, pi1_reverse) %>%
  pivot_longer(cols = c(pi1_forward, pi1_reverse),
               names_to = "direction", values_to = "pi1") %>%
  mutate(direction = recode(direction,
                            pi1_forward = "Forward",
                            pi1_reverse = "Reverse"))

# Split the data by direction and filter to neurons only
pi1_forward_tbl <- pi1_long %>% 
  filter(direction == "Forward",
    (
      (ref_cell_type == "Exc" & (cell_type == "Glu-UL" | cell_type == "Glu-DL" | cell_type == "GABA")) |
        (ref_cell_type == "Inh" & (cell_type == "Glu-UL" | cell_type == "Glu-DL" | cell_type == "GABA"))
    )
  )


# Function to generate a heatmap
plot_fugita_heatmap <- function(df) {
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
    labs(x = "Adult Cell Type", y = "Prenatal Cell type") +
    theme_minimal(base_size = 13) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
      axis.text.y = element_text(face = "bold"),
      panel.grid = element_blank(),
      plot.margin = margin(50, 50, 50, 50, unit = "pt")
    ) 
}

# Generate heatmaps
pi1_fugita_heatmap <- plot_fugita_heatmap(pi1_forward_tbl)

### --- beta correlation plt -----
beta_files <- c(
#  "All"    = paste0(beta_dir, "replication_beta_correlation.tsv"),
  "GABA"   = paste0(beta_dir, "GABA_beta_cor_single_tbl.tsv"),
  "Glu-DL" = paste0(beta_dir, "Glu-DL_beta_cor_single_tbl.tsv"),
  "Glu-UL" = paste0(beta_dir, "Glu-UL_beta_cor_single_tbl.tsv")
)

make_beta_cor_plot <- function(tbl_path, gene_lookup, ct = NULL, label_genes = NULL) {
  
  paired_betas <- read_tsv(tbl_path, show_col_types = FALSE) |>
    separate(key, into = c("snp", "gene")) |>
    left_join(gene_lookup, by = "gene") |>
    distinct()
  
  ## Correlation
  cor_val <- cor(
    paired_betas$beta_my,
    paired_betas$beta_fugita,
    use = "complete.obs",
    method = "pearson"
  )
  
  cor_label <- sprintf("r = %.2f", cor_val)
  
  ## Linear model
  model <- lm(beta_my ~ beta_fugita, data = paired_betas)
  
  paired_betas <- paired_betas |>
    mutate(
      .resid = resid(model),
      .dist  = abs(.resid),
      discordant = sign(beta_my) != sign(beta_fugita) &
        !is.na(beta_my) & !is.na(beta_fugita),
      label_text = ifelse(is.na(symbol) | symbol == "NA", gene, symbol),
      is_strong_outlier = .dist > quantile(.dist, 0.92, na.rm = TRUE),
      should_label = if (!is.null(label_genes)) label_text %in% label_genes else FALSE
    )
  
  ggplot(paired_betas, aes(x = beta_fugita, y = beta_my)) +
    geom_point(alpha = 0.2, color = "grey70", size = 1.8) +
    geom_point(
      data = filter(paired_betas, is_strong_outlier & discordant),
      color = "#d32f2f", size = 3, alpha = 0.5
    ) +
    geom_point(
      data = filter(paired_betas, is_strong_outlier & !discordant),
      color = "#1976d2", size = 3, alpha = 0.5
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
      "text", x = Inf, y = Inf, label = cor_label,
      hjust = -0.1, vjust = 1.4, size = 5
    ) +
    labs(
      title = element_blank(),
      subtitle = element_blank(),
      x = substitute(beta ~ "(Adult" ~ x ~ ")", list(x = ct[2])),
      y = substitute(beta ~ "(Prenatal" ~ x ~ ")", list(x = ct[1]))
    ) +
    theme_minimal(base_size = 13) +
    theme(
      plot.margin = margin(10, 10, 10, 10, unit = "pt")
    )
}

beta_gaba_plt <- make_beta_cor_plot(beta_files[["GABA"]], gene_lookup, 
                                    c("GABA", "Inh"), 'ABCC8')
beta_gluDL_plt <- make_beta_cor_plot(beta_files[["Glu-DL"]], gene_lookup, 
                                     c("Glu-DL", "Exc"), 'ABI3BP')
beta_gluUL_plt <- make_beta_cor_plot(beta_files[["Glu-UL"]], gene_lookup, 
                                     c("Glu-UL", "Exc"), 'ABCC8')

### --- plot -----
# Final plot
top_row <- plot_grid(pie_chart, upset_plt, labels = c("A", "B"),
  label_size = 24, ncol = 2,rel_widths = c(1, 1.3))

# Stack heatmaps
heatmaps_stacked <- plot_grid(
  pi1_int_heatmap, #+ theme(plot.margin = margin(5,10,5,5,"pt")),
  pi1_fugita_heatmap, #+ theme(plot.margin = margin(5,10,5,5,"pt")),
  ncol = 1,
  rel_heights = c(1, 0.8),   # C slightly taller
  labels = c("C", "D"),
  label_size = 24,
  align = "v"
)

# Stack betas 
betas_stacked <- plot_grid(
  beta_gluUL_plt, beta_gluDL_plt, beta_gaba_plt,
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