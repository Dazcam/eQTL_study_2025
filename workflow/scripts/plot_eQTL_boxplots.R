# Load libraries
library(ggplot2)
library(cowplot)
library(readr)
library(dplyr)
library(stringr)

# Define palette
custom_palette <- c(
  "RG" = "#FF5959",
  "ExN-UL" = "#00B6EB",
  "InN" = "#3CBB75FF",
  "ExN-DL" = "#CEE5FD",
  "Endo-Peri" = "#B200ED",
  "MG" = "#F58231",
  "OPC" = "#FDE725FF"
)

# Define cell types
cell_types <- c("ExN-UL", "ExN-DL", "InN", "RG", "MG", "OPC", "Endo-Peri")

# Read and combine all available cell-type eQTL CSVs
eqtl_data <- lapply(cell_types, function(ct) {
  file_path <- file.path("~/Desktop/eQTL_study_2025/workflow/reports/11PLOTTING/boxplots_py", 
                         paste0("eqtl_data_", ct, "_rs11191424_AS3MT.csv"))
  if (file.exists(file_path)) {
    read_csv(file_path, show_col_types = FALSE) %>%
      mutate(CellType = ct)
  } else {
    message("File not found for ", ct)
    NULL
  }
}) %>%
  bind_rows()

# Extract REF and ALT alleles from the dataset (they're constant for this SNP)
ref_allele <- unique(eqtl_data$REF)[1]
alt_allele <- unique(eqtl_data$ALT)[1]

# Derive genotype labels (e.g., GG, GA, AA)
eqtl_data <- eqtl_data %>%
  mutate(
    Alleles = str_replace_all(GT, "\\|", "/"),
    Genotype = case_when(
      GT %in% c("0|0", "0/0") ~ paste0(ref_allele, ref_allele),
      GT %in% c("0|1", "1|0", "0/1", "1/0") ~ paste0(ref_allele, alt_allele),
      GT %in% c("1|1", "1/1") ~ paste0(alt_allele, alt_allele),
      TRUE ~ NA_character_
    ),
    Genotype = factor(
      Genotype,
      levels = c(
        paste0(ref_allele, ref_allele),
        paste0(ref_allele, alt_allele),
        paste0(alt_allele, alt_allele)
      )
    )
  )

# Compute counts per genotype per cell type
counts <- eqtl_data %>%
  group_by(CellType, Genotype) %>%
  summarise(n = n(), .groups = "drop")

# Generate labels with subscript counts for the x-axis
# e.g., GG[n==12] → expression(GG[n==12]) for use with parse = TRUE in ggplot
label_df <- counts %>%
  mutate(
    label = paste0(Genotype, "[n==", n, "]")
  )

# Function to get parsed x labels for a given cell type
get_labels_for_celltype <- function(ct) {
  df <- label_df %>% filter(CellType == ct)
  labs <- df$label
  names(labs) <- df$Genotype
  return(labs)
}

# Create one plot per cell type
plots <- lapply(cell_types, function(ct) {
  df <- eqtl_data %>% filter(CellType == ct)
  if (nrow(df) == 0) return(NULL)
  
  ggplot(df, aes(x = Genotype, y = Expression, color = 'black', fill = CellType)) +
    geom_boxplot(width = 0.5, outlier.size = 2) +
   # geom_jitter(width = 0.15, size = 2, alpha = 0.4) +
    scale_color_manual(values = custom_palette) +
    scale_fill_manual(values = custom_palette) +
    scale_x_discrete(
      labels = parse(text = get_labels_for_celltype(ct))
    ) +
    labs(title = ct, x = "Genotype", y = "ASM3T Expression (log2)") +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      axis.title.x = element_blank(),
      axis.title.y = element_text(size = 12),
      axis.text = element_text(size = 10),
      legend.position = "none",
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
      plot.margin = margin(t = 5, r = 5, b = 20, l = 5)
    ) + ylim(0,6) 
})

# Combine all panels
final_plot <- plot_grid(plotlist = plots[!sapply(plots, is.null)], nrow = 1, align = "hv")

final_plot_with_label <- ggdraw(final_plot) +
  draw_label(
    "rs11191424",
    x = 0.5, y = 0.02,                # centered near the bottom
    fontface = "bold", size = 14
  )

# Save output
ggsave("~/Desktop/eQTL_study_2025/workflow/reports/11PLOTTING/AS3MT_eqtl_grouped_plot_genotypes.pdf",
       plot = final_plot, width = 15, height = 5, useDingbats = FALSE)

final_plot


