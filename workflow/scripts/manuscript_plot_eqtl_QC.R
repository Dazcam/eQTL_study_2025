#--------------------------------------------------------------------------------------
#
#    Generate eQTL QC plt for manuscript
#
#--------------------------------------------------------------------------------------

# A: Bar chart L1 and L2 cell counts
# B: Cell counts vs. eGenes
# C: Proportion of eGene sharing across L1 and L2 cell types (main 4 cell types)
# D: eQTL distance to TSS denisty plots
# F: eQTL in OCR Enrichment plot

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

message("\n\nGenerating eQTL QC plot for the manuscript ...")

# -------------------------------------------------------------------------------------
library(tidyverse)
library(cowplot)
library(ggrepel)

# --- Set variables
in_dir <- snakemake@params[['in_dir']] 
ziffra_dir <- snakemake@params[['ziffra_dir']] 
out_file <- snakemake@output[[1]] 

#in_dir <- "../results/05TENSORQTL/tensorqtl_perm/"
#out_dir <- "../results/13MANUSCRIPT_PLOTS_TABLES/"

eqtl_list <- list()
genPC <- 4
norm_method <- 'quantile'

cell_types <- c("Glu-UL", "Glu-DL", "NPC", "GABA",
                "Endo-Peri", "OPC", "MG",
                "Glu-UL-0", "Glu-UL-1", "Glu-UL-2",
                "Glu-DL-0", "Glu-DL-1", "Glu-DL-2",
                "GABA-0", "GABA-1", "GABA-2",
                "NPC-0", "NPC-1", "NPC-2")

cell_counts <-  c(271180, 136070, 128017, 103598, 19842, 13794, 
                  7297, 123028, 87340, 60812, 49120, 44808, 42142, 
                  59350, 47255, 20371, 33789, 32047, 17799)

expPC_map <- c(
  "Glu-UL"     = 50,
  "Glu-DL"     = 40,
  "GABA"       = 30,
  "NPC"        = 30,
  "MG"         = 30,
  "OPC"        = 30,
  "Endo-Peri"  = 30,
  "Glu-UL-0"     = 30,
  "Glu-UL-1"     = 50,
  "Glu-UL-2"     = 40,
  "Glu-DL-0"     = 40,
  "Glu-DL-1"     = 40,
  "Glu-DL-2"     = 30,
  "GABA-0"       = 30,
  "GABA-1"       = 40,
  "GABA-2"       = 40,
  "NPC-0"        = 50,
  "NPC-1"        = 30,
  "NPC-2"        = 20
)

cell_types_L1 <- c("Glu-UL", "Glu-DL", "GABA", "NPC")
cell_types_L2 <- c("Glu-UL-0", "Glu-UL-1", "Glu-UL-2",
                   "Glu-DL-0", "Glu-DL-1", "Glu-DL-2",
                   "NPC-0", "NPC-1", "NPC-2",
                   "GABA-0", "GABA-1", "GABA-2")

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

# --- Specify y-axis order
cell_order <- c(
  "Glu-DL", paste0("Glu-DL-", 0:2), 
  "Glu-UL", paste0("Glu-UL-", 0:2),
  "NPC", paste0("NPC-", 0:2),
  "GABA", paste0("GABA-", 0:2),
  "Endo-Peri", "OPC", "MG"
)

# --- Base theme
base_theme <- theme_minimal(base_size = 14) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold"),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, face = "bold"),
    panel.spacing.x = unit(2, "lines"),
  )


# --- Density plot 
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

# --- Specify y-axis order
eqtls_all <- bind_rows(eqtl_list) |>
  mutate(level = if_else(str_detect(cell_type, "\\d"), 2L, 1L),
         cell_type = fct_relevel(cell_type, cell_order),
         cell_type = fct_rev(cell_type)) 


# --- Set main type for colouring
eqtls_all <- eqtls_all |>
 # mutate(cell_type = factor(cell_type, levels = rev(cell_order))) |>
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

density_plt <- ggplot(eqtls_all, aes(x=distance_kb, color=cell_type)) +
  geom_density() +
  geom_vline(xintercept = 0, linetype="dashed", color = "black") +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.6) +
  #scale_color_manual(values=custom_palette) +
  #scale_fill_manual(values=custom_palette) +
  scale_x_continuous(limits=c(-600,600), breaks=seq(-600,600,200)) +
  labs(x="Distance to TSS (kb)", y="Density") +
  base_theme +
  theme(plot.margin = unit(c(1, 1, 1, 1), "cm"))

# --- eQTL counts
# need to pull these programatically from final object / or check
cell_counts <- tibble(
  cell_type = cell_types,
  n_cells   = cell_counts
)

gene_cell <- eqtls_all %>%
  distinct(phenotype_id, cell_type)

egene_per_cell <- gene_cell %>% 
  count(cell_type) %>%
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
  

eqtl_cnt_plt <- ggplot(egene_per_cell, aes(x=cell_type, y=n, fill=main_type)) +
  geom_col(width = 0.7, colour = 'black') +
  scale_fill_manual(values = custom_palette) +
  labs(x="Cell Type", y="Number of eGenes") +
  coord_flip() +
  base_theme +
  geom_hline(yintercept = 0, color = "black") +
  theme(plot.margin = unit(c(1, 1, 1, 1), "cm"),
        legend.position = "none")

# --- Scatter plot
scatter_dat <- egene_per_cell %>%
  left_join(cell_counts, by="cell_type") |>
  mutate(level = if_else(str_detect(cell_type, "\\d"), 2L, 1L))

# Correlation
cor_test <- cor.test(scatter_dat$n_cells, scatter_dat$n, method = "pearson")
cor_label <- sprintf("r = %.3f\np = %.2e", cor_test$estimate, cor_test$p.value)

scatter_plt <- ggplot(scatter_dat, aes(x=n_cells, y=n)) +
  geom_point(aes(color = main_type, shape = factor(level), size = factor(level))) +
  # geom_text_repel(aes(label=cell_type, color=main_type), 
  #                 size=4, 
  #                 show.legend=FALSE,
  #                 max.overlaps=20, 
  #                 force=2,      # pushes labels apart
  #                 nudge_y=100,  # pushes labels up a bit
  # ) +
  scale_color_manual(values=custom_palette) +
  scale_shape_manual(values = c("1" = 16, "2" = 17), guide = "none") +  # 16: filled circle, 17: filled triangle
  scale_size_manual(values = c("1" = 4, "2" = 2.5), guide = "none") +
  scale_x_continuous(
    breaks = c(0, 1e5, 2e5, 3e5),
    labels = c("0", "100K", "200K", "300K")
  ) +
  labs(x="Nuclei count", y="eGenes (FDR < 0.05)") +
  annotate("text", x = -Inf, y = Inf, label = cor_label,
           hjust = -0.2, vjust = 1.5, size = 3, fontface = "bold") +
  base_theme +
  theme(
    plot.margin = unit(c(1, 1, 1, 1), "cm"),
    panel.grid.minor = element_blank(),
    legend.position = "none")


# --- Combined
read_eGenes <- function(cell_type, in_dir, norm_method, genPC, expPC) {
  log_dir  <- file.path(in_dir, paste0(cell_type, "_", norm_method, "_genPC_", genPC, "_expPC_", expPC))
  log_file <- file.path(log_dir, paste0(cell_type, "_", norm_method, "_perm.cis_qtl.txt.gz"))
  if (!file.exists(log_file)) return(tibble(cell_type = cell_type, eGene = character(0)))
  
  read_tsv(log_file, show_col_types = FALSE) %>%
    filter(!is.na(qval) & qval < 0.05) %>%
    distinct(phenotype_id) %>%
    mutate(cell_type = cell_type)
}

eqtl_all <- map_dfr(c(cell_types_L1, cell_types_L2),
                    read_eGenes,
                    in_dir = in_dir,
                    norm_method = norm_method,
                    genPC = genPC,
                    expPC = expPC)

# Separate into lists
egenes_L1 <- eqtl_all %>% filter(cell_type %in% cell_types_L1)
egenes_L2 <- eqtl_all %>% filter(cell_type %in% cell_types_L2)

# Group L2 subtypes under their parent L1 type
egenes_L2 <- egenes_L2 %>%
  mutate(parent = sub("-[^-]+$", "", cell_type))

# Compare sets per parent
comparison_tbl <- map_dfr(cell_types_L1, function(ct) {
  egenes_ct_L1 <- egenes_L1 %>% filter(cell_type == ct) %>% pull(phenotype_id) %>% unique()
  egenes_ct_L2 <- egenes_L2 %>% filter(parent == ct) %>% pull(phenotype_id) %>% unique()
  
  shared     <- intersect(egenes_ct_L1, egenes_ct_L2)
  unique_L1  <- setdiff(egenes_ct_L1, egenes_ct_L2)
  unique_L2  <- setdiff(egenes_ct_L2, egenes_ct_L1)
  
  tibble(
    cell_type = ct,
    unique_L1 = length(unique_L1),
    shared    = length(shared),
    unique_L2 = length(unique_L2)
  )
})

# Plot
comparison_long <- comparison_tbl %>%
  pivot_longer(cols = c(unique_L1, shared, unique_L2),
               names_to = "category",
               values_to = "count") %>%
  mutate(category = factor(category, 
                           levels = c("unique_L2", "shared", "unique_L1"))) %>%  # <-- set order
  group_by(cell_type) %>%
  mutate(prop = count / sum(count))

combined_plt <- ggplot(comparison_long, aes(x = prop, y = fct_rev(cell_type), fill = category)) +
  geom_col(position = "stack", color = "black", width = 0.7) +
  geom_text(aes(label = count),
            position = position_stack(vjust = 0.5),
            size = 3, color = "black", fontface = "bold") +
  scale_fill_manual(values = c("unique_L1" = "#66C2A5",
                               "shared"    = "#B2DF8A",
                               "unique_L2" = "gray70"),
                    labels = c("unique_L1" = "Cell type",
                               "shared"    = "Shared",
                               "unique_L2" = "Sub-cluster"),
                    breaks = c("unique_L1", "shared", "unique_L2")) +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(x = "Proportion of eGenes",
       y = 'Cell Type',
       fill = NULL) +
  base_theme +
  theme(legend.position = "top",
        plot.margin = unit(c(1, 1, 1, 1), "cm"))

# --- Ziffra
ziffra_tbl <- read_tsv(paste0(ziffra_dir, 'ziffra_overlaps_primary.tsv')) |>
  mutate(clean_peak = str_replace_all(peak_cell_type, "_MACSpeaks", ""),
         test = paste(eqtl_cell_type, clean_peak, sep = ' in '))
  
# --- Specify y-axis order
cell_order <- c(
  sort(grep("^Glu-UL",  ziffra_tbl$test, value = TRUE)),
  sort(grep("^Glu-DL",  ziffra_tbl$test, value = TRUE)),
  sort(grep("^GABA",    ziffra_tbl$test, value = TRUE)),
  sort(grep("^NPC",     ziffra_tbl$test, value = TRUE)),
  sort(grep("^OPC",       ziffra_tbl$test, value = TRUE)),
  sort(grep("^MG",       ziffra_tbl$test, value = TRUE)),
  sort(grep("^Endo-",       ziffra_tbl$test, value = TRUE))
) |> unique()

# --- Assign main cluster to subclusters
ziffra_tbl <- ziffra_tbl |>
  mutate(test = factor(test, levels = cell_order)) |>
  mutate(main_type = case_when(
    str_detect(test, "Glu-UL") ~ "Glu-UL",
    str_detect(test, "Glu-DL") ~ "Glu-DL",
    str_detect(test, "GABA") ~ "GABA",
    str_detect(test, "NPC") ~ "NPC",
    str_detect(test, "OPC") ~ "OPC",
    str_detect(test, "MG") ~ "MG",
    str_detect(test, "Endo-Peri") ~ "Endo-Peri",
    TRUE ~ test
  ))

# --- Base theme
ziffra_plt <- ziffra_tbl |>
  ggplot(aes(x = test, y = fold_enrichment, fill = main_type)) +  
  geom_col(width = 0.7, colour = 'black') +
  theme_minimal(base_size = 14) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.6) +
#  geom_vline(yintercept = 1, linetype = "dashed", color = "black") +
  geom_hline(yintercept = 1, linetype = "dotted", color = "black") +
  scale_fill_manual(values = custom_palette) +
  labs(
    x = "Enrichment Test",
    y = "Fold Enrichment"
  ) +
  coord_cartesian(ylim = c(0, 7)) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
    legend.position = "none",
    plot.margin = margin(t = 10, r = 60, b = 10, l = 50)
  )

# --- Combine A-D into a 2x2 grid
top_rows <- plot_grid(
  eqtl_cnt_plt, scatter_plt, 
  combined_plt, density_plt, 
  labels = c('A', 'B', 'C', 'D'), 
  ncol = 2, 
  label_size = 24
)

# --- Final assembly: Top rows + Ziffra on bottom
final_plt <- plot_grid(
  top_rows, 
  ziffra_plt, 
  labels = c('', 'E'), # Label E for the bottom plot
  ncol = 1, 
  rel_heights = c(2, 1), # Adjust this ratio to give Ziffra more/less vertical space
  label_size = 24
)

ggsave(
  filename = out_file,
  plot = final_plt,
  width = 10,  
  height = 14,     
  units = "in",
  device = "pdf",
  dpi = 300
)

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
