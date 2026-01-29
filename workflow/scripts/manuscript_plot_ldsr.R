#--------------------------------------------------------------------------------------
#
#    Generate MaxCPP S-LDSR plt for manuscript
#
#--------------------------------------------------------------------------------------

# A: S-LDSR bar chart for L1 populations
# B: S-LDSR bar chart for L2 populations

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

message("\n\nGenerating MaxCPP S-LDSR plot for the manuscript ...")

# -------------------------------------------------------------------------------------
library(tidyverse)
library(cowplot)
library(ggrepel)

# --- Set variables
in_dir <- snakemake@params[['in_dir']] 
out_file <- snakemake@output[[1]] 

# in_dir <- "../results/09SLDSR/strat_bl_v12/"
# out_dir <- "../results/13MANUSCRIPT_PLOTS_TABLES/"

cell_types <- c("Glu-UL", "Glu-DL", "NPC", "GABA",
                "Endo-Peri", "OPC", "MG",
                "Glu-UL-0", "Glu-UL-1", "Glu-UL-2",
                "Glu-DL-0", "Glu-DL-1", "Glu-DL-2",
                "GABA-0", "GABA-1", "GABA-2",
                "NPC-0", "NPC-1", "NPC-2")

L1_cell_types <- c("Glu-UL", "Glu-DL", "NPC", "GABA", "Endo-Peri", "OPC", "MG")
L2_cell_types <- c(
  "Glu-UL-0", "Glu-UL-1", "Glu-UL-2",
  "Glu-DL-0", "Glu-DL-1", "Glu-DL-2",
  "GABA-0", "GABA-1", "GABA-2",
  "NPC-0", "NPC-1", "NPC-2"
)

disorders <- c('scz', 'bpd', 'mdd', 'adhd', 'ocd')


# --- Custom colour palette
custom_palette <- c(
  'Glu-UL' = '#4363d8',
  'Glu-DL' = '#00B6EB',
  'NPC' = '#FF5959',
  'GABA' = '#3CBB75FF',
  'Endo-Peri' = '#B200ED',
  'MG' = '#F58231',
  'OPC' = '#FDE725FF'
)


# --- Read and prepare table
ldsr_tbl <- read_tsv(paste0(in_dir, 'ldsr_strat_hg38_bl_v12.maxCPP.summary.tsv')) |>
  separate(Category, into = c("annot", "file_name"), sep = "/", remove = FALSE) |>
  separate(file_name, into = c("cell_type", "suffix"), sep = "_hg38", remove = TRUE) |>
  mutate(disorder = toupper(str_extract(suffix, "scz|bpd|mdd|adhd|ocd"))) |>
  select(-Category, -suffix) |>
  relocate(cell_type, disorder) |>
  mutate(ldsr = if_else(`Coefficient_z-score` > 0,
                        -log10(pnorm(`Coefficient_z-score`, lower.tail = FALSE)), 0)) |>
  mutate(level = if_else(str_detect(cell_type, "\\d"), 2L, 1L)) |>
  filter(disorder != 'PTSD') |>
  mutate(disorder = recode(disorder,
                           "SCZ" = "Schizophrenia",
                           "BPD" = "Bipolar Disorder")) |>
  mutate(disorder = factor(disorder, 
                           levels = c("Schizophrenia", "Bipolar Disorder", "MDD", "ADHD", "OCD")))

# --- Specify y-axis order
cell_order <- c(
  "Glu-UL",   sort(grep("^Glu-UL-",  ldsr_tbl$cell_type, value = TRUE)),
  "Glu-DL",   sort(grep("^Glu-DL-",  ldsr_tbl$cell_type, value = TRUE)),
  "GABA",     sort(grep("^GABA-",    ldsr_tbl$cell_type, value = TRUE)),
  "NPC",      sort(grep("^NPC-",     ldsr_tbl$cell_type, value = TRUE)),
  "OPC", "MG", "Endo-Peri"
) |> unique()

# --- Assign main cluster to subclusters
ldsr_tbl <- ldsr_tbl |>
  mutate(cell_type = factor(cell_type, levels = rev(cell_order))) |>
  mutate(main_type = case_when(
    str_detect(cell_type, "Glu-UL") ~ "Glu-UL",
    str_detect(cell_type, "Glu-DL") ~ "Glu-DL",
    str_detect(cell_type, "GABA") ~ "GABA",
    str_detect(cell_type, "NPC") ~ "NPC",
    str_detect(cell_type, "OPC") ~ "OPC",
    str_detect(cell_type, "MG") ~ "MG",
    str_detect(cell_type, "Endo-Peri") ~ "Endo-Peri",
    TRUE ~ cell_type
  ))

# --- Significance thresholds
thresh_level1 <- -log10(0.05 / length(L1_cell_types)) 
thresh_level2 <- -log10(0.05 / length(L1_cell_types))  
thresh_nominal <- -log10(0.05)

# --- Base theme
base_theme <- theme_minimal(base_size = 12) +
  theme(
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(size = 10),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold"),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, face = "bold"),
    panel.spacing.x = unit(2, "lines"),
  )

# --- Plot A: Level 1 
plot_A <- ldsr_tbl |>
  filter(level == 1) |>
  ggplot(aes(x = ldsr, y = cell_type, fill = main_type)) +  
  geom_col(width = 0.7, colour = 'black') +
  facet_grid(. ~ disorder, scales = "free_y", space = "free_y") +
  theme_minimal(base_size = 12) +
  base_theme +
  geom_vline(xintercept = 0, color = "black", linewidth = 0.6) +
  geom_vline(xintercept = thresh_level1, linetype = "dashed", color = "black") +
  geom_vline(xintercept = thresh_nominal, linetype = "dotted", color = "black") +
  base_theme +
  scale_fill_manual(values = custom_palette) +
  labs(
    x = expression(-log[10](P)),
    y = "Cell Type"
  ) +
  coord_cartesian(xlim = c(0, 6))

# --- Plot B: Level 2
plot_B <- ldsr_tbl |>
  filter(level == 2) |>
  ggplot(aes(x = ldsr, y = cell_type, fill = main_type)) + 
  geom_col(width = 0.7, colour = 'black') +
  facet_grid(. ~ disorder, scales = "free_y", space = "free_y") +
  scale_fill_manual(values = custom_palette) +
  geom_vline(xintercept = 0, color = "black", linewidth = 0.6) +
  geom_vline(xintercept = thresh_level2, linetype = "dashed", color = "black") +
  geom_vline(xintercept = thresh_nominal, linetype = "dotted", color = "black") +
  theme_minimal(base_size = 12) +
  base_theme +
  scale_fill_manual(values = custom_palette) +
  labs(
    x = expression(-log[10](P)),
    y = "Cell Type"
  ) +
  coord_cartesian(xlim = c(0, 6))

# --- Combine into final figure
final_plot <- plot_grid(
  plot_A, plot_B,
  labels = c("A", "B"),
  label_size = 20,
  ncol = 1,
  align = "v",
  rel_heights = c(0.4, 0.6)  # adjust if Level 2 has many more rows
)

# 1. PDF (best for publications: vector, scalable, crisp text)
ggsave(
  filename = paste0(out_file),
  plot = final_plot,
  width = 10,      # inches, good for full-page width in papers
  height = 6,     # adjust if needed (depends on number of subclusters)
  units = "in",
  device = "pdf",
  dpi = 300
)

# 2. High-resolution PNG (for quick viewing, presentations, or supplements)
# ggsave(
#   filename = paste0(out_dir, "ldsr_stratified_levels_publication.png"),
#   plot = final_plot,
#   width = 10,
#   height = 6,
#   units = "in",
#   dpi = 300,
#   device = "png",
#   bg = "white"
# )

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------