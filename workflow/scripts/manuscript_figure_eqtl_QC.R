#--------------------------------------------------------------------------------------
#
#    eQTL discovery and QC across cell populations (Fig 2)
#
#--------------------------------------------------------------------------------------
#
# Pipeline:  17MANUSCRIPT-FIGURES | Rule: eqtl_qc_fig
#            Upstream:   05 tensorqtl_perm (final norm method, genPC and expPCs);
#                        07 ziffra_summary (OCR enrichment of lead eQTL);
#                        16 cell_counts_tbl (nuclei counts)
#            Downstream: Fig 2
#
# Purpose:   Panels for the 19 L1 and L2 cell populations:
#            A  eGenes (FDR < 0.05) per cell population
#            B  eGenes against nuclei count, with Pearson r
#            C  eGenes unique to L1, unique to L2 or shared, for the 4 major L1
#               types and their L2 subpopulations
#            D  Density of lead eQTL distance to the TSS
#            E  Fold enrichment of lead eQTL in Ziffra 2021 snATAC-seq peaks
#            Manuscript cell type labels (config manuscript_labels)
#
# Inputs:    in_dir         tensorQTL permutation results directory
#            ziffra_tbl     07 Ziffra OCR enrichment summary
#            nuclei_counts  16 cell_counts_tbl (nuclei per L1 and L2 cell type)
#            cell_types, l1_types, exp_pc_map, geno_pc, norm_method, eqtl_fdr  (config)
#
# Outputs:   out_file       Fig 2 (PDF)
#
# Notes:     Every panel reads each cell type's permutation results at its
#            exp_pc_map value.
#
#--------------------------------------------------------------------------------------

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

message("\n\nPlotting eQTL QC figure ...")

##  Load packages, functions and variables  -------------------------------------------
library(tidyverse)
library(readxl)
library(cowplot)

# Load functions (relabel_cell_types)
source(file.path(snakemake@scriptdir, "functions.R"))

# Inputs, outputs and settings
in_dir             <- snakemake@params[['in_dir']]
ziffra_tbl_file    <- snakemake@input[['ziffra_tbl']]
nuclei_counts_file <- snakemake@input[['nuclei_counts']]
manuscript_labels  <- snakemake@params[["labels"]]
expPC_map          <- unlist(snakemake@params[["exp_pc_map"]])
genPC              <- snakemake@params[["geno_pc"]]
norm_method        <- snakemake@params[["norm_method"]]
eqtl_fdr           <- as.numeric(snakemake@params[["eqtl_fdr"]])
out_file           <- snakemake@output[[1]]

# L1 and L2 cell types (pipeline labels, config order); pseudotime bins excluded
cell_types    <- unlist(snakemake@params[["cell_types"]])
cell_types    <- cell_types[!str_detect(cell_types, "-to-")]
l1_types      <- unlist(snakemake@params[["l1_types"]])
cell_types_L2 <- setdiff(cell_types, l1_types)
# L1 types with L2 subclusters (panel C)
cell_types_L1 <- l1_types[l1_types %in% sub("-[^-]+$", "", cell_types_L2)]

eqtl_list <- list()

# Print what each variable is set to
message("\nVariables")
message("============================")
tibble(
  variable = c("in_dir", "ziffra_tbl_file", "nuclei_counts_file", "labels", "norm_method",
               "genPC", "eqtl_fdr", "cell_types", "cell_types_L1 (panel C)",
               "cell_types_L2", "expPC_map", "out_file"),
  value    = c(in_dir, ziffra_tbl_file, nuclei_counts_file, paste(names(manuscript_labels), unlist(manuscript_labels), sep = " -> ", collapse = ", "),
               norm_method, genPC, eqtl_fdr, paste(cell_types, collapse = ", "),
               paste(cell_types_L1, collapse = ", "), paste(cell_types_L2, collapse = ", "),
               paste(names(expPC_map), expPC_map, sep = ":", collapse = ", "), out_file)) |>
  knitr::kable(format = "simple") |>
  cat(sep = "\n")
message("============================\n")

# Custom colour palette
custom_palette <- c(
  'Glu-UL' = '#4363d8',
  'Glu-DL' = '#00B6EB',
  'NPC' = '#FF5959',
  'GABA' = '#3CBB75FF',
  'Endo-Peri' = '#B200ED',
  'MG' = '#F58231',
  'OPC' = '#FDE725FF'
)

# y-axis order (pipeline labels; relabelled for plotting below)
cell_order <- c(
  "Glu-UL", paste0("Glu-UL-", 0:2),
  "Glu-DL", paste0("Glu-DL-", 0:2), 
  "GABA", paste0("GABA-", 0:2),
  "NPC", paste0("NPC-", 0:2),
  "OPC", "MG", "Endo-Peri"
)

# Relabelled versions of the order/palette vectors, for use in plots only
cell_order_new <- relabel_cell_types(cell_order, manuscript_labels)
custom_palette_new <- custom_palette
names(custom_palette_new) <- relabel_cell_types(names(custom_palette_new), manuscript_labels)

# Base theme
base_theme <- theme_minimal(base_size = 14) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold"),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, face = "bold"),
    panel.spacing.x = unit(2, "lines"),
  )

##  Panel D: distance to TSS (and data for A and B)  ----------------------------------
for (cell_type in cell_types) {
  
  expPC <- expPC_map[[cell_type]]
  
  f <- file.path(
    in_dir,
    paste0(cell_type, "_", norm_method, "_genPC_", genPC, "_expPC_", expPC),
    paste0(cell_type, "_", norm_method, "_perm.cis_qtl.txt.gz")  
  )
  
  if (file.exists(f)) {
    dat <- read_tsv(f, show_col_types = FALSE) |>
      filter(!is.na(start_distance)) |>
      mutate(distance_kb = start_distance / 1000) |>
      filter(qval < eqtl_fdr)
    
    if (nrow(dat) > 0) {
      eqtl_list[[cell_type]] <- dat |> mutate(cell_type = cell_type)
    }
  }
}

# y-axis order
eqtls_all <- bind_rows(eqtl_list) |>
  mutate(level = if_else(str_detect(cell_type, "\\d"), 2L, 1L),
         cell_type = fct_relevel(cell_type, cell_order),
         cell_type = fct_rev(cell_type)) 

# Main type for colouring
eqtls_all <- eqtls_all |>
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

# Manuscript labels for plotting (cell_type and main_type)
eqtls_all <- eqtls_all |>
  mutate(cell_type = relabel_cell_types(as.character(cell_type), manuscript_labels),
         cell_type = fct_relevel(cell_type, cell_order_new),
         cell_type = fct_rev(cell_type),
         main_type = relabel_cell_types(main_type, manuscript_labels))

density_plt <- ggplot(eqtls_all, aes(x=distance_kb, color=cell_type)) +
  geom_density() +
  geom_vline(xintercept = 0, linetype="dashed", color = "black") +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.6) +
  scale_x_continuous(limits=c(-600,600), breaks=seq(-600,600,200)) +
  labs(x="Distance to TSS (kb)", y="Density") +
  base_theme +
  theme(plot.margin = unit(c(1, 1, 1, 1), "cm"))

##  Panel A: eGenes per cell population  ----------------------------------------------
# Nuclei per cell type from 16 cell_counts_tbl (manuscript labels; L2 names are
# prefixed with their L1 cluster if the sheet does not already include it)
nuclei_l1 <- read_excel(nuclei_counts_file, sheet = "L1") |>
  transmute(label = `L1 cluster`, n_cells = Nuclei)
nuclei_l2 <- read_excel(nuclei_counts_file, sheet = "L2") |>
  transmute(label = if_else(str_starts(`L2 subcluster`, fixed(paste0(`L1 cluster`, "-"))),
                            `L2 subcluster`, paste0(`L1 cluster`, "-", `L2 subcluster`)),
            n_cells = Nuclei)
nuclei_lookup <- deframe(bind_rows(nuclei_l1, nuclei_l2))

cell_counts <- tibble(
  cell_type = cell_types,
  n_cells   = unname(nuclei_lookup[relabel_cell_types(cell_types, manuscript_labels)])
)
if (anyNA(cell_counts$n_cells)) {
  stop("No nuclei count for: ", paste(cell_counts$cell_type[is.na(cell_counts$n_cells)], collapse = ", "))
}
message("Nuclei counts per cell type:")
cell_counts |> knitr::kable(format = "simple") |> cat(sep = "\n")

gene_cell <- eqtls_all |>
  distinct(phenotype_id, cell_type)

egene_per_cell <- gene_cell |> 
  count(cell_type) |>
  mutate(main_type = case_when(
    str_detect(cell_type, "Glu-A") ~ "Glu-A",
    str_detect(cell_type, "Glu-B") ~ "Glu-B",
    str_detect(cell_type, "GABA") ~ "GABA",
    str_detect(cell_type, "NPC") ~ "NPC",
    str_detect(cell_type, "OPC") ~ "OPC",
    str_detect(cell_type, "MG") ~ "MG",
    str_detect(cell_type, "Endo-Peri") ~ "Endo-Peri",
    TRUE ~ cell_type
  ))
  

eqtl_cnt_plt <- ggplot(egene_per_cell, aes(x=cell_type, y=n, fill=main_type)) +
  geom_col(width = 0.7, colour = 'black') +
  scale_fill_manual(values = custom_palette_new) +
  labs(x="Cell Type", y="Number of eGenes") +
  coord_flip() +
  base_theme +
  geom_hline(yintercept = 0, color = "black") +
  theme(plot.margin = unit(c(1, 1, 1, 1), "cm"),
        legend.position = "none")

##  Panel B: eGenes against nuclei count  ----------------------------------------------
# cell_counts uses pipeline labels, so relabel it to join egene_per_cell
scatter_dat <- egene_per_cell |>
  left_join(cell_counts |> mutate(cell_type = relabel_cell_types(cell_type, manuscript_labels)), by="cell_type") |>
  mutate(level = if_else(str_detect(cell_type, "\\d"), 2L, 1L))

# Correlation
cor_test <- cor.test(scatter_dat$n_cells, scatter_dat$n, method = "pearson")
cor_label <- sprintf("r = %.3f\np = %.2e", cor_test$estimate, cor_test$p.value)

scatter_plt <- ggplot(scatter_dat, aes(x=n_cells, y=n)) +
  geom_point(aes(color = main_type, shape = factor(level), size = factor(level))) +
  scale_color_manual(values=custom_palette_new) +
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


##  Panel C: eGenes unique to L1, unique to L2 or shared  -----------------------------
read_eGenes <- function(cell_type, in_dir, norm_method, genPC, expPC) {
  log_dir  <- file.path(in_dir, paste0(cell_type, "_", norm_method, "_genPC_", genPC, "_expPC_", expPC))
  log_file <- file.path(log_dir, paste0(cell_type, "_", norm_method, "_perm.cis_qtl.txt.gz"))
  if (!file.exists(log_file)) return(tibble(cell_type = cell_type, eGene = character(0)))
  
  read_tsv(log_file, show_col_types = FALSE) |>
    filter(!is.na(qval) & qval < eqtl_fdr) |>
    distinct(phenotype_id) |>
    mutate(cell_type = cell_type)
}

# Each cell type at its own exp_pc_map value
eqtl_all <- map_dfr(c(cell_types_L1, cell_types_L2),
                    ~ read_eGenes(.x, in_dir = in_dir, norm_method = norm_method,
                                  genPC = genPC, expPC = expPC_map[[.x]]))

# Separate into lists
egenes_L1 <- eqtl_all |> filter(cell_type %in% cell_types_L1)
egenes_L2 <- eqtl_all |> filter(cell_type %in% cell_types_L2)

# Group L2 subtypes under their parent L1 type
egenes_L2 <- egenes_L2 |>
  mutate(parent = sub("-[^-]+$", "", cell_type))

# Compare sets per parent
comparison_tbl <- map_dfr(cell_types_L1, function(ct) {
  egenes_ct_L1 <- egenes_L1 |> filter(cell_type == ct) |> pull(phenotype_id) |> unique()
  egenes_ct_L2 <- egenes_L2 |> filter(parent == ct) |> pull(phenotype_id) |> unique()
  
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

# Manuscript labels, in the order Glu-A, Glu-B, GABA, NPC
comparison_tbl <- comparison_tbl |>
  mutate(cell_type = relabel_cell_types(cell_type, manuscript_labels),
         cell_type = factor(cell_type, levels = relabel_cell_types(cell_types_L1, manuscript_labels)))

# Plot
comparison_long <- comparison_tbl |>
  pivot_longer(cols = c(unique_L1, shared, unique_L2),
               names_to = "category",
               values_to = "count") |>
  mutate(category = factor(category, 
                           levels = c("unique_L2", "shared", "unique_L1"))) |>
  group_by(cell_type) |>
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

##  Panel E: enrichment in Ziffra 2021 snATAC-seq peaks  ------------------------------
# Relabel the source columns before building the compound `test` label
ziffra_tbl <- read_tsv(ziffra_tbl_file, show_col_types = FALSE) |>
  mutate(clean_peak = str_replace_all(peak_cell_type, "_MACSpeaks", ""),
         eqtl_cell_type = relabel_cell_types(eqtl_cell_type, manuscript_labels),
         clean_peak = relabel_cell_types(clean_peak, manuscript_labels),
         test = paste(eqtl_cell_type, clean_peak, sep = ' in '))
  
# y-axis order
cell_order <- c(
  sort(grep("^Glu-A",  ziffra_tbl$test, value = TRUE)),
  sort(grep("^Glu-B",  ziffra_tbl$test, value = TRUE)),
  sort(grep("^GABA",    ziffra_tbl$test, value = TRUE)),
  sort(grep("^NPC",     ziffra_tbl$test, value = TRUE)),
  sort(grep("^OPC",       ziffra_tbl$test, value = TRUE)),
  sort(grep("^MG",       ziffra_tbl$test, value = TRUE)),
  sort(grep("^Endo-",       ziffra_tbl$test, value = TRUE))
) |> unique()

# Main type for colouring
ziffra_tbl <- ziffra_tbl |>
  mutate(test = factor(test, levels = cell_order)) |>
  mutate(main_type = case_when(
    str_detect(test, "Glu-A") ~ "Glu-A",
    str_detect(test, "Glu-B") ~ "Glu-B",
    str_detect(test, "GABA") ~ "GABA",
    str_detect(test, "NPC") ~ "NPC",
    str_detect(test, "OPC") ~ "OPC",
    str_detect(test, "MG") ~ "MG",
    str_detect(test, "Endo-Peri") ~ "Endo-Peri",
    TRUE ~ test
  ))

ziffra_plt <- ziffra_tbl |>
  ggplot(aes(x = test, y = fold_enrichment, fill = main_type)) +  
  geom_col(width = 0.7, colour = 'black') +
  theme_minimal(base_size = 14) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.6) +
  geom_hline(yintercept = 1, linetype = "dotted", color = "black") +
  scale_fill_manual(values = custom_palette_new) +
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

##  Assemble and save  ----------------------------------------------------------------
# A-D in a 2x2 grid
top_rows <- plot_grid(
  eqtl_cnt_plt, scatter_plt, 
  combined_plt, density_plt, 
  labels = c('A', 'B', 'C', 'D'), 
  ncol = 2, 
  label_size = 24
)

# A-D above E
final_plt <- plot_grid(
  top_rows, 
  ziffra_plt, 
  labels = c('', 'E'),
  ncol = 1, 
  rel_heights = c(2, 1),
  label_size = 24
)

message("Writing: ", out_file)
ggsave(
  filename = out_file,
  plot = final_plt,
  width = 10,  
  height = 14,     
  units = "in",
  device = "pdf",
  dpi = 300
)

message("Done.")

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------

