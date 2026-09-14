#--------------------------------------------------------------------------------------
#
#    Generate MaxCPP S-LDSR plt for manuscript
#
#--------------------------------------------------------------------------------------

# A: S-LDSR bar chart for L1 populations
# B: S-LDSR bar chart for L2 populations
# C: S-LDSR bar chart for pseudotime trajectory bins

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

# NOTE: unused elsewhere in this script (kept from the original for reference);
# left in old-name form since it doesn't drive any file path or plot here.
cell_types <- c("Glu-UL", "Glu-DL", "NPC", "GABA",
                "Endo-Peri", "OPC", "MG",
                "Glu-UL-0", "Glu-UL-1", "Glu-UL-2",
                "Glu-DL-0", "Glu-DL-1", "Glu-DL-2",
                "GABA-0", "GABA-1", "GABA-2",
                "NPC-0", "NPC-1", "NPC-2")

# --- L1/L2/trajectory group definitions (new Glu-A/Glu-B names), used only to
# size the per-facet Bonferroni thresholds below
L1_cell_types <- c("Glu-A", "Glu-B", "NPC", "GABA", "Endo-Peri", "OPC", "MG")
L2_cell_types <- c(
  "Glu-A-0", "Glu-A-1", "Glu-A-2",
  "Glu-B-0", "Glu-B-1", "Glu-B-2",
  "GABA-0", "GABA-1", "GABA-2",
  "NPC-0", "NPC-1", "NPC-2"
)
trajectory_cell_types <- c(
  paste0("NPC-to-Glu-A-Bin", 1:4),
  paste0("NPC-to-Glu-B-Bin", 1:4)
)

disorders <- c('scz', 'bpd', 'mdd', 'adhd', 'ocd')


# --- Custom colour palette (Trajectory colours match those used in the SLDSR
# report Rmd, kept distinct from every main-cluster colour)
custom_palette <- c(
  'Glu-A' = '#4363d8',
  'Glu-B' = '#00B6EB',
  'NPC' = '#FF5959',
  'GABA' = '#3CBB75FF',
  'Endo-Peri' = '#B200ED',
  'MG' = '#F58231',
  'OPC' = '#FDE725FF',
  'Trajectory-to-Glu-A' = '#E91E8C',
  'Trajectory-to-Glu-B' = '#2F4F4F'
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
  # --- Trajectory (pseudotime) rows are flagged BEFORE the digit-based L1/L2
  # split, since names like "NPC-to-Glu-UL-Q4-Bin1" contain digits and would
  # otherwise be mis-classified as a level-2 subcluster.
  mutate(is_trajectory = str_detect(cell_type, "-to-")) |>
  mutate(level = case_when(
    is_trajectory ~ 3L,
    str_detect(cell_type, "\\d") ~ 2L,
    TRUE ~ 1L
  )) |>
  filter(disorder != 'PTSD') |>
  mutate(disorder = recode(disorder,
                           "SCZ" = "Schizophrenia",
                           "BPD" = "Bipolar Disorder")) |>
  mutate(disorder = factor(disorder, 
                           levels = c("Schizophrenia", "Bipolar Disorder", "MDD", "ADHD", "OCD")))

# --- Relabel cell type names for plotting: Glu-UL -> Glu-A, Glu-DL -> Glu-B,
# and drop "-Q4-" from pseudotime trajectory bin names, e.g.
# "NPC-to-Glu-UL-Q4-Bin1" -> "NPC-to-Glu-A-Bin1"
relabel_cell_type <- function(x) {
  x <- str_replace(x, "Glu-UL", "Glu-A")
  x <- str_replace(x, "Glu-DL", "Glu-B")
  x <- str_replace(x, "-Q4-", "-")
  x
}

ldsr_tbl <- ldsr_tbl |>
  mutate(cell_type = relabel_cell_type(cell_type))

# --- Specify y-axis order for Plots A/B (main clusters + subclusters).
# Subcluster patterns are anchored ("-[0-9]+$") so they can't accidentally
# re-catch trajectory names (e.g. an unanchored "^NPC-" grep would otherwise
# also match "NPC-to-Glu-A-Bin1").
cell_order_main <- c(
  "Glu-A",   sort(grep("^Glu-A-[0-9]+$",  ldsr_tbl$cell_type, value = TRUE)),
  "Glu-B",   sort(grep("^Glu-B-[0-9]+$",  ldsr_tbl$cell_type, value = TRUE)),
  "GABA",     sort(grep("^GABA-[0-9]+$",    ldsr_tbl$cell_type, value = TRUE)),
  "NPC",      sort(grep("^NPC-[0-9]+$",     ldsr_tbl$cell_type, value = TRUE)),
  "OPC", "MG", "Endo-Peri"
) |> unique()

# --- Specify y-axis order for Plot C (pseudotime trajectory bins):
# Glu-A bins before Glu-B bins, each in Bin1-4 order
cell_order_traj <- c(
  sort(grep("^NPC-to-Glu-A-Bin[0-9]+$", ldsr_tbl$cell_type, value = TRUE)),
  sort(grep("^NPC-to-Glu-B-Bin[0-9]+$", ldsr_tbl$cell_type, value = TRUE))
) |> unique()

cell_order <- c(cell_order_main, cell_order_traj)

# --- Assign main cluster to subclusters / trajectory bins (trajectory check
# happens first, since e.g. "NPC-to-Glu-A-Bin1" would otherwise match the
# "NPC" or "Glu-A" pattern below)
ldsr_tbl <- ldsr_tbl |>
  mutate(cell_type = factor(cell_type, levels = rev(cell_order))) |>
  mutate(main_type = case_when(
    is_trajectory & str_detect(cell_type, "to-Glu-A") ~ "Trajectory-to-Glu-A",
    is_trajectory & str_detect(cell_type, "to-Glu-B") ~ "Trajectory-to-Glu-B",
    str_detect(cell_type, "Glu-A") ~ "Glu-A",
    str_detect(cell_type, "Glu-B") ~ "Glu-B",
    str_detect(cell_type, "GABA") ~ "GABA",
    str_detect(cell_type, "NPC") ~ "NPC",
    str_detect(cell_type, "OPC") ~ "OPC",
    str_detect(cell_type, "MG") ~ "MG",
    str_detect(cell_type, "Endo-Peri") ~ "Endo-Peri",
    TRUE ~ as.character(cell_type)
  ))

# --- Significance thresholds (per-facet Bonferroni, sized from the group
# counts above; trajectory gets its own threshold from its 8 bins)
thresh_level1 <- -log10(0.05 / length(L1_cell_types)) 
thresh_level2 <- -log10(0.05 / length(L2_cell_types))  
thresh_level3 <- -log10(0.05 / length(trajectory_cell_types))
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
  coord_cartesian(xlim = c(0, 4))

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
  base_theme +
  scale_fill_manual(values = custom_palette) +
  labs(
    x = expression(-log[10](P)),
    y = "Cell Type"
  ) +
  coord_cartesian(xlim = c(0, 4))

# --- Plot C: Pseudotime trajectory bins
plot_C <- ldsr_tbl |>
  filter(level == 3) |>
  ggplot(aes(x = ldsr, y = cell_type, fill = main_type)) +
  geom_col(width = 0.7, colour = 'black') +
  facet_grid(. ~ disorder, scales = "free_y", space = "free_y") +
  scale_fill_manual(values = custom_palette) +
  geom_vline(xintercept = 0, color = "black", linewidth = 0.6) +
  geom_vline(xintercept = thresh_level3, linetype = "dashed", color = "black") +
  geom_vline(xintercept = thresh_nominal, linetype = "dotted", color = "black") +
  base_theme +
  scale_fill_manual(values = custom_palette) +
  labs(
    x = expression(-log[10](P)),
    y = "Cell Type"
  ) +
  coord_cartesian(xlim = c(0, 4))

# --- Combine into final figure
# rel_heights are approximate, sized by row count per panel (7 / 12 / 8 rows) -
# adjust if the balance looks off once rendered.
final_plot <- plot_grid(
  plot_A, plot_B, plot_C,
  labels = c("A", "B", "C"),
  label_size = 20,
  ncol = 1,
  align = "v",
  rel_heights = c(0.3, 0.45, 0.32)
)

# 1. PDF 
ggsave(
  filename = out_file,
  plot = final_plot,
  width = 10,      # inches, good for full-page width in papers
  height = 8.5,     # increased to accommodate the new Plot C panel
  units = "in",
  device = "pdf",
  dpi = 300
)

# 2. High-res PNG 
# ggsave(
#   filename = out_file,
#   plot = final_plot,
#   width = 10,
#   height = 8.5,
#   units = "in",
#   dpi = 300,
#   device = "png",
#   bg = "white"
# )

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
