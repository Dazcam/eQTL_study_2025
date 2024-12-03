#--------------------------------------------------------------------------------------
#
#    Parse_testing - Basic Seuart processing
#
#--------------------------------------------------------------------------------------

## Info  ------------------------------------------------------------------------------

#  Using Seurat 5 primarily, but also bioconductor packages for QC 

# Load packages  ----------------------------------------------------------------------
library(BPCells)
library(Seurat)
library(SeuratObject)
library(SeuratDisk)
#library(Azimuth)
library(dplyr)
library(Matrix)
library(ggplot2)
library(scater)
library(scDblFinder)
library(tidyverse)
library(harmony)
library(clustifyr)

library(scCustomize)

library(dlookr)
library(flextable)
library(egg)

# Set variables  ----------------------------------------------------------------------
root_dir <- '/Users/darren/Desktop/eQTL_study_2025/'
config_dir <- paste0(root_dir, 'config/')
results_dir <- paste0(root_dir, 'results/02PARSE/')
seurat_dir <- paste0(root_dir, 'results/03SEURAT/')
resources_dir <- paste0(root_dir, 'resources/')
sheets_dir <- paste0(resources_dir, 'sheets/')
bp_cells_dir <- paste0(seurat_dir, 'bp_cells_dir/')
combine_dir <- paste0(results_dir, 'combine_plate1/')
mat_dir <- paste0(combine_dir, 'all-sample/DGE_filtered/')
report_dir <- paste0(combine_dir, 'all-sample/report/')
script_dir <- paste0(root_dir, 'workflow/scripts/')
markdown_doc <- paste0(script_dir, 'eQTL_2025_seurat.Rmd')
markdown_html <- paste0(results_dir, 'eQTL_2025_seurat.html')
markdown_dir <- paste0(root_dir, 'workflow/reports/markdown/')
public_dir <- paste0(resources_dir, 'public_datasets/')
cameron_dir <- paste0(public_dir, 'cameron_2023/')

resolution_set <- c(0.3, 0.5, 0.8)

# Load functions ----------------------------------------------------------------------
source('~/Desktop/parse_testing/workflow/scripts/parse_testing_functions.R')
#source('~/Desktop/eQTL_study_2025/workflow/scripts/eQTL_2025_functions.R')
options(ggrepel.max.overlaps = Inf)
options(future.globals.maxSize = 8000 * 1024^2) # Required for label transfer

## Load data  -------------------------------------------------------------------------
# mat <- ReadParseBio(mat_dir) # Crashes locally with 1M cells
mat <- ReadParseBio(mat_dir)

# Load additional metadata
sample_fcx_meta <- read_csv(paste0(sheets_dir, 'FC_samples_pcw.csv'))
sample_meta <- readxl::read_excel(paste0(sheets_dir, 'Fetal single cell eQTL final samples 29-11-24.xlsx'), 
                                  sheet = 'Final cortex samples') |>
  dplyr::rename(sample = ...1) |>
  dplyr::select(sample, PCW, Sex) |>
  mutate(sample = str_replace_all(sample, " \\((A|a)\\)", "_A")) |>
  print(n=Inf)

# ID Claire's samples and add sex and pcw info to meta data
parse_meta <- read.csv(paste0(mat_dir, "cell_metadata.csv"), row.names = 1) |>
  mutate(region = case_when(
    str_detect(sample, pattern = "_WGE") ~ "GE",
    str_detect(sample, pattern = "_Hip") ~ "Hip",
    str_detect(sample, pattern = "_Thal") ~ "Tha",
    str_detect(sample, pattern = "18184") ~ "FC_11pcw", # Rm and not 2nd Trim
    .default = 'CTX'),
    sample_id = str_replace(sample, "sample_", ""),
    claire_sample = ifelse(str_detect(region, pattern = "GE|Hip|Tha|FC_11pcw"), T, F)) |> 
  left_join(sample_fcx_meta, by = join_by('sample_id' == 'sample')) |>
  left_join(sample_meta, by = join_by('sample_id' == 'sample')) %>%
  mutate(PCW_combined = coalesce(as.character(pcw), as.character(PCW))) %>%
  dplyr::select(-pcw, -PCW) |>
  dplyr::rename(pcw = PCW_combined)

seurat_obj <- CreateSeuratObject(mat, meta.data = parse_meta)
rm(mat, parse_meta, sample_meta, sample_fcx_meta)



# brain_data <- open_matrix_anndata_hdf5(
#   path = paste0(mat_dir, 'anndata.h5ad')
# )
# 
# # Write the matrix to a directory
# write_matrix_dir(
#   mat = brain_data,
#   dir = paste(seurat_dir, 'bp_cells_dir'))
# 
# brain_mat <- open_matrix_dir(dir = paste(seurat_dir, 'bp_cells_dir'))
# 
# # Gene IDs have an hg38 suffix which affects Azimuth
# rownames(brain_mat) <- str_replace(rownames(brain_mat), '_hg38', '')
# 
# # Azimuth not needed? Removed tons of genes anyway
# # brain.mat <- Azimuth:::ConvertEnsembleToSymbol(mat = brain.mat, species = "human") 
# 
# 
# meta <- read.csv(paste0(mat_dir, "cell_metadata.csv"), row.names = 1)
# seurat_obj <- CreateSeuratObject(brain_mat, meta.data = meta)

## Parse QC ---------------------------------------------------------------------------
# Agg table
# agg_samp_tbl <- read_csv(paste0(combine_dir, 'agg_samp_ana_summary.csv')) 
# agg_samp_tbl <- as_tibble(cbind(nms = names(agg_samp_tbl), t(agg_samp_tbl))) 
# agg_samp_tbl <- agg_samp_tbl |> 
#   janitor::row_to_names(row_number = 1) |>
#   dplyr::rename(sample_id = statistic) |>
#   mutate_at(vars(-sample_id), as.numeric) |>
#   dplyr::select(!ends_with("Q30") & !ends_with("valid") & !ends_with("fraction") & !ends_with("read1")) |>
#   dplyr::select(!sample_well_count) |>
#   mutate(sample_id = str_replace(sample_id, "sample_", "")) |>
#   relocate(sample_id, number_of_cells, hg38_median_tscp_per_cell, hg38_median_genes_per_cell) |>
#   arrange(hg38_median_genes_per_cell) 
# 
# 
# # Step 1: Calculate Z-scores for number_of_tscp
# agg_samp_tbl <- agg_samp_tbl %>%
#   mutate(z_score = (number_of_tscp - mean(number_of_tscp)) / sd(number_of_tscp),
#          outlier = abs(z_score) > 2)  # Define outliers as Z-scores greater than 2 (can adjust threshold)
# 
# # Wells round 1
# sample_key_tbl <- read_delim(paste0(config_dir, 'sample-list.txt'), col_names = FALSE) |>
#   dplyr::rename(sample = X1, well = X2)
# sample_loading_tbl <- read_tsv(paste0(config_dir, 'sample_loading.tsv'), col_names = FALSE) |>
#   dplyr::rename(well = X1, loadings = X2)
# counts_rnd1_tbl <- read_csv(paste0(report_dir, 'cell_counts_by_rnd1_well.csv')) |>
#   dplyr::rename(cell_cnt = count)
# tscp_median_rnd1_tbl <- read_csv(paste0(report_dir, 'tscp_median_by_rnd1_well.csv')) |>
#   dplyr::rename(tscp_med_cnt = count)
# 
# rnd1_tbl <- sample_key_tbl |>
#   inner_join(counts_rnd1_tbl) |>
#   inner_join(tscp_median_rnd1_tbl) |>
#   inner_join(sample_loading_tbl) |>
#   mutate(Row = gsub("[^a-zA-Z]", "", well), .before=2) |>
#   mutate(Column = parse_number(well), .before=2) 
# 
# # Plot heatmaps
# rnd1_hmap <- plot_parse_rnd1_heatmaps(rnd1_tbl)
# 
# # Pull out outliers if necessary
# rnd1_tbl |>
#   filter(sample %in% c(14402, 14635, 11654, 17109))

## QC ----------------------------------------------------------------------------------
counts_preFilt_tbl <- tibble('Stage' = 'Pre-filter', 
                             'Cells' = c(ncol(seurat_obj)),
                             'Genes' = c(nrow(seurat_obj)))

# # ID Claire's samples
# seurat_obj@meta.data <- seurat_obj@meta.data |>
#   mutate(region = case_when(
#     str_detect(sample, pattern = "_WGE") ~ "GE",
#     str_detect(sample, pattern = "_Hip") ~ "Hip",
#     str_detect(sample, pattern = "_Thal") ~ "Tha",
#     str_detect(sample, pattern = "18184") ~ "FC_11pcw", # Rm and not 2nd Trim
#     .default = 'CTX'),
#     sample_id = str_replace(sample, "sample_", ""),
#     claire_sample = ifelse(str_detect(region, pattern = "GE|Hip|Tha|FC_11pcw"), T, F))

# Remove Claire's samples
seurat_obj <- subset(seurat_obj, subset = claire_sample == FALSE)
counts_FiltClaire_tbl <- tibble('Stage' = 'CT-filt', 
                                'Cells' = c(ncol(seurat_obj)),
                                'Genes' = c(nrow(seurat_obj)))

# Median counts per sample
tscp_boxplot <- seurat_obj@meta.data |>
  as_tibble() |>
  arrange(desc(sample_id)) |>
  mutate(sample_id = as_factor(sample_id)) |>
  ggplot(aes(x = log10(nCount_RNA), y = sample_id)) +
  geom_boxplot(outlier.size = 0.1, outlier.alpha = 0.1) +
  theme_linedraw() +
  theme(legend.position = "none") +
  ylab(NULL) 

gene_boxplot <- seurat_obj@meta.data |>
  as_tibble() |>
  arrange(desc(sample_id)) |>
  mutate(sample_id = as_factor(sample_id)) |>
  ggplot(aes(x = log10(nFeature_RNA), y = sample_id)) +
  geom_boxplot(outlier.size = 0.1, outlier.alpha = 0.1) +
  theme_linedraw() +
  theme(axis.text.y = element_blank()) +
  geom_vline(xintercept = log10(500), linetype="dotted", color = "red") +
  geom_vline(xintercept = log10(1000), linetype="dotted", color = "red") +
  ylab(NULL) 

tscp_gene_boxplt <- egg::ggarrange(tscp_boxplot, gene_boxplot, ncol = 2)

median_cnt_tbl <- seurat_obj@meta.data |> 
  tibble::as_tibble() |> 
  dplyr::group_by(sample_id) |> 
  summarise(genes_per_cell_median = median(nFeature_RNA),
            tscp_per_cell_median = median(nCount_RNA)) |>
  arrange(genes_per_cell_median) |>
  print(n = Inf) |>
  flextable() |> 
  color(~ genes_per_cell_median < 1500, ~ genes_per_cell_median, color = "red") |>
  set_caption(caption = "Median counts table") |>
  theme_vanilla()

## Apply Bryois filters - Several samples with 
# Cell level filtering - 36,718
samples_low_cnts_tbl <- seurat_obj@meta.data |> 
  as_tibble() |>
  arrange(nCount_RNA) |>
  filter(nCount_RNA < 500 | nFeature_RNA < 300) |>
  group_by(sample_id) |>
  summarize(count = n()) |>
  arrange(desc(count)) |>
  print(n = Inf) |>
  flextable() |> 
  set_caption(caption = "Cells to rm (per sample) with reads < 500 OR genes < 300 ") |>
  theme_vanilla()

# Cell counts per sample after filter
cell_cnts_post_filt_tbl <- seurat_obj@meta.data |> 
  as_tibble() |>
  arrange(nCount_RNA) |>
  filter(nCount_RNA >= 500 & nFeature_RNA >= 300) |>
  group_by(sample_id) |>
  summarize(count = n()) |>
  arrange(count) |>
  print(n = Inf) |>
  flextable() |> 
  set_caption(caption = "Counts per sample after filter") |>
  theme_vanilla()

# Median counts after filter
median_cnt_post_filt_tbl <- seurat_obj@meta.data |> 
  as_tibble() |>
  arrange(nCount_RNA) |>
  filter(nCount_RNA >= 500 & nFeature_RNA >= 300) |>
  dplyr::group_by(sample_id) |> 
  summarise(genes_per_cell_median = median(nFeature_RNA),
            tscp_per_cell_median = median(nCount_RNA)) |>
  arrange(genes_per_cell_median) |>
  print(n = Inf) |>
  flextable() |> 
  set_caption(caption = "Median counts after filter") |>
  theme_vanilla()


# Median counts per sample boxplots post filter
tscp_post_filt_boxplot <- seurat_obj@meta.data |>
  as_tibble() |>
  filter(nCount_RNA >= 500 & nFeature_RNA >= 300) |>
  arrange(desc(sample_id)) |>
  mutate(sample_id = as_factor(sample_id)) |>
  ggplot(aes(x = log10(nCount_RNA), y = sample_id)) +
  geom_boxplot(outlier.size = 0.1, outlier.alpha = 0.1) +
  theme_linedraw() +
  theme(legend.position = "none") +
  geom_vline(xintercept = log10(500), linetype="dotted", color = "red") +
  geom_vline(xintercept = log10(1000), linetype="dotted", color = "red") +
  ylab(NULL) 

gene_post_filt_boxplot <- seurat_obj@meta.data |>
  as_tibble() |>
  filter(nCount_RNA >= 500 & nFeature_RNA >= 300) |>
  arrange(desc(sample_id)) |>
  mutate(sample_id = as_factor(sample_id)) |>
  ggplot(aes(x = log10(nFeature_RNA), y = sample_id)) +
  geom_boxplot(outlier.size = 0.1, outlier.alpha = 0.1) +
  theme_linedraw() +
  theme(axis.text.y = element_blank()) +
  geom_vline(xintercept = log10(500), linetype="dotted", color = "red") +
  geom_vline(xintercept = log10(1000), linetype="dotted", color = "red") +
  ylab(NULL) 

tscp_gene_post_filt_boxplt <- egg::ggarrange(tscp_post_filt_boxplot, gene_post_filt_boxplot, ncol = 2)

# Subset 1
seurat_obj <- subset(seurat_obj, subset = nFeature_RNA >= 300 & nCount_RNA >= 500)
counts_subset1_tbl <- tibble('Stage' = 'Sub 1', 
                             'Cells' = c(ncol(seurat_obj)),
                             'Genes' = c(nrow(seurat_obj)))

# Convert to single cell experiment
sce_obj <- SingleCellExperiment(list(counts = as(seurat_obj[["RNA"]]$counts, "dgCMatrix")),
                                colData = seurat_obj@meta.data)

# Identify cell outliers and apply filters to sce object
sce_obj <- get_cell_outliers(sce_obj, 3, 'higher', 5, 5, 
                             'nCount_RNA', 'nFeature_RNA', 'sample')
outlier_plt <- create_outlier_plots(sce_obj, sce_obj$umi_outlier, sce_obj$genes_outlier, sce_obj$mito_outlier,
                                    sce_obj$ribo_outlier, 'nCount_RNA', 'nFeature_RNA', 'sample')

# Identify doublets
sce_obj <- scDblFinder(sce_obj, dbr = 0.0003, samples = sce_obj$sample_id)
table(sce_obj$scDblFinder.class)
sce_obj$doublet_logic <- ifelse(sce_obj$scDblFinder.class == "doublet", TRUE, FALSE)
sce_obj$all_outliers <- sce_obj@colData |>
  as_tibble() |>
  mutate(all_outlier = as.logical(pmax(cell_outlier, doublet_logic))) |>
  dplyr::select(all_outlier) %>%
  pull()

outlier_cnts_tbl <- tibble(
  measure = c('umi', 'genes', 'mito', 'ribo', 'doublets', 'total'), 
  count = c(sum(sce_obj$umi_outlier), sum(sce_obj$genes_outlier), 
            sum(sce_obj$mito_outlier),sum(sce_obj$ribo_outlier), 
            sum(sce_obj$doublet_logic), sum(sce_obj$all_outliers))  
) 

# Subset Seurat object
seurat_obj <- subset_seurat_object(seurat_obj, sce_obj$all_outliers)
counts_subset2_tbl <- tibble('Stage' = 'Sub 2', 
                             'Cells' = c(ncol(seurat_obj)),
                             'Genes' = c(nrow(seurat_obj)))
filt_all_counts_tbl <- rbind(counts_preFilt_tbl, 
                             counts_FiltClaire_tbl, 
                             counts_subset1_tbl, 
                             counts_subset2_tbl) |>
  flextable() |> 
  set_caption(caption = "Counts after filters") |>
  theme_vanilla()
rm(sce_obj)

# Final counts after QC
# Median counts after filter
median_cnt_final_tbl <- seurat_obj@meta.data |> 
  as_tibble() |>
  arrange(nCount_RNA) |>
  dplyr::group_by(sample_id) |> 
  summarise(genes_per_cell_median = median(nFeature_RNA),
            tscp_per_cell_median = median(nCount_RNA)) |>
  arrange(genes_per_cell_median) |>
  print(n = Inf) |>
  flextable() |> 
  set_caption(caption = "Median counts final") |>
  theme_vanilla()


# Median counts per sample boxplots post filter
tscp_final_boxplot <- seurat_obj@meta.data |>
  as_tibble() |>
  arrange(desc(sample_id)) |>
  mutate(sample_id = as_factor(sample_id)) |>
  ggplot(aes(x = log10(nCount_RNA), y = sample_id)) +
  geom_boxplot(outlier.size = 0.1, outlier.alpha = 0.1) +
  theme_linedraw() +
  theme(legend.position = "none") +
  geom_vline(xintercept = log10(500), linetype="dotted", color = "red") +
  geom_vline(xintercept = log10(1000), linetype="dotted", color = "red") +
  ylab(NULL) 

gene_final_boxplot <- seurat_obj@meta.data |>
  as_tibble() |>
  arrange(desc(sample_id)) |>
  mutate(sample_id = as_factor(sample_id)) |>
  ggplot(aes(x = log10(nFeature_RNA), y = sample_id)) +
  geom_boxplot(outlier.size = 0.1, outlier.alpha = 0.1) +
  theme_linedraw() +
  theme(axis.text.y = element_blank()) +
  geom_vline(xintercept = log10(500), linetype="dotted", color = "red") +
  geom_vline(xintercept = log10(1000), linetype="dotted", color = "red") +
  ylab(NULL) 

tscp_gene_final_boxplt <- egg::ggarrange(tscp_post_filt_boxplot, gene_post_filt_boxplot, ncol = 2)

## Seurat Processing ------------------------------------------------------------------
# Run basic pipeline
seurat_obj <- run_seurat_process(seurat_obj)

# Post pipeline QCs
# Check PCs
#pca_plot <- DimHeatmap(seurat_obj, dims = 20:40, cells = 500, balanced = TRUE)

# Create resolution plots for comparison
res_plotlist <- create_resolution_plotlist(seurat_obj, resolution = resolution_set,
                                           meta_id = 'sample_id', gene_set = general_genes)

# Plot QCs
cluster_qc_plot <- create_cluster_qc_plot(seurat_obj, 30, 'sample_id', 'RNA_snn_res.0.3')
meta_col_cluster_cnts <- get_meta_col_counts(seurat_obj, 'RNA_snn_res.0.3') |>
  arrange(desc(n)) |>
  print(n = Inf) |>
  flextable() |> 
  set_caption(caption = "Seurat Clusters: cell counts") |>
  theme_vanilla()

## Integration - using Harmony
## Is this necessary? 
## Integration - using Harmony on sploit object / maybe not necessary
## Crashing near the end of process (maybe need to split before noraml Seurat process)
# seurat_obj[["RNA"]] <- split(seurat_obj[["RNA"]], f = seurat_obj$sample)
# seurat_object <- run_integration(seurat_obj, 'harmony', 30, resolution_set)
# seurat_obj <- JoinLayers(seurat_obj)
# ####
seurat_obj <- run_harmony(seurat_obj, 'harmony', 'sample_id', 30, resolution_set)

#save.image(paste0(results_dir, 'parse_tetsing.RData'))
#load(paste0(results_dir, 'parse_tetsing.RData'))

integration_plotlist <- create_integration_plotlist(seurat_obj, 
                                                    meta_id = 'sample_id', 
                                                    dims = 30,
                                                    reduction = resolution_set) 

meta_col_harmony_cnts <- get_meta_col_counts(seurat_obj, 'harmony_clusters_0.3') |>
  arrange(desc(n)) |>
  print(n = Inf) |>
  flextable() |> 
  set_caption(caption = "Harmony Clusters: cell counts") |>
  theme_vanilla()

## Compare FC sample cluster distribution ---------------------------------------------
seurat_obj@meta.data <- seurat_obj@meta.data |>
  mutate(fcx_sample = case_when(
    str_detect(sample_id, pattern = "_FC") ~ TRUE,
    .default = FALSE))

fc_sample_ids <- seurat_obj@meta.data |> filter(fcx_sample == TRUE) |> distinct(sample_id) |> pull()
non_fc_sample_ids <- seurat_obj@meta.data |> filter(fcx_sample == FALSE) |> distinct(sample_id) |> pull()
seurat_obj$sample_id_fcx_fact <- seurat_obj@meta.data |>
  mutate(sample_id = factor(sample_id, levels = c(fc_sample_ids, non_fc_sample_ids))) |>
  dplyr::select(sample_id)

seurat_obj$pcw <- seurat_obj@meta.data |>
  mutate(pcw = factor(pcw)) 

fcx_umap_plt <- DimPlot(seurat_obj, reduction = "umap.harmony",
                        label = TRUE, group.by = "harmony_clusters_0.3", cols = discreet_cols_n23) 
fcx_feat_plt <- FeaturePlot(seurat_obj, 'fcx_sample', reduction = 'umap.harmony') 
cowplot::plot_grid(fcx_umap_plt, fcx_feat_plt)
  
fcx_fact_barplt <- create_proportion_barplot(seurat_obj, 'sample_id_fcx_fact', meta_id = 'harmony_clusters_0.3') +
  scale_fill_manual(values = discreet_cols_n23) +
  theme(axis.text.x = element_text(angle = 90)) 

fcx_pcw_barplt <- create_proportion_barplot(seurat_obj, 'pcw', meta_id = 'harmony_clusters_0.3') +
  scale_fill_manual(values = discreet_cols_n23) +
  theme(axis.text.x = element_text(angle = 90)) 

seurat_obj@meta.data |>
  select(harmony_clusters_0.3, sample_id) |>
  group_by()

meta_col_cluster_cnts <- get_meta_col_counts(seurat_obj, 'RNA_snn_res.0.3') |>
  arrange(desc(n)) |>
  print(n = Inf) |>
  flextable() |> 
  set_caption(caption = "Harmony Clusters: cell counts") |>
  theme_vanilla()

seurat_obj@meta.data %>%
  tibble::as_tibble(rownames = 'cell_id') %>%
  dplyr::select(any_of(c(meta_col))) %>%
  dplyr::group_by(.data[[meta_col]]) %>%
  tally() %>%
  mutate(freq = n / sum(n))


## Cell type assignment  --------------------------------------------------------------
# Gene lists  -----
biol_psych_vln_plot <- egg::ggarrange(create_stacked_vln_plot(seurat_obj, 'harmony_clusters_0.3', pfc_features, 'PFC harmony 0.3'),
                                  create_stacked_vln_plot(seurat_obj, 'harmony_clusters_0.3', hip_features, 'Hip harmony 0.3'),
                                  create_stacked_vln_plot(seurat_obj, 'harmony_clusters_0.3', wge_features, 'GE harmony 0.3'),
                                  create_stacked_vln_plot(seurat_obj, 'harmony_clusters_0.3', tha_features, 'Tha harmony 0.3'),
                                  create_stacked_vln_plot(seurat_obj, 'harmony_clusters_0.3', cer_features, 'Cer harmony 0.3'))

public_vln_plot <- cowplot::plot_grid(create_stacked_vln_plot(seurat_obj, 'harmony_clusters_0.3', Nowakowski_Fig4A_genes, 'Nowak 2017 harmony 0.3'),
                                      create_stacked_vln_plot(seurat_obj, 'harmony_clusters_0.3', Pouliodakis_fig1c_genes, 'Poulio 2019 harmony 0.3'))


ct_genes_vln <- create_stacked_vln_plot(seurat_obj, set_ident = 'harmony_clusters_0.3', genes = claire_genes)
poulio_genes_vln <- create_stacked_vln_plot(seurat_obj, set_ident = 'harmony_clusters_0.3', genes = Pouliodakis_fig1c_genes)
nowak_genes_vln <- create_stacked_vln_plot(seurat_obj, set_ident = 'harmony_clusters_0.3', genes = Nowakowski_Fig4A_genes)
camer_genes_vln <- create_stacked_vln_plot(seurat_obj, set_ident = 'harmony_clusters_0.3', genes = pfc_features)
fcx_adult_genes_vln <- create_stacked_vln_plot(seurat_obj, set_ident = 'harmony_clusters_0.3', genes = fcx_genes)
general_genes_vln <- create_stacked_vln_plot(seurat_obj, set_ident = 'harmony_clusters_0.3', genes = general_genes)
exN_genes_vln <- create_stacked_vln_plot(seurat_obj, set_ident = 'harmony_clusters_0.3', genes = c(exN_genes, inN_genes))
InN_genes_vln <- create_stacked_vln_plot(seurat_obj, set_ident = 'harmony_clusters_0.3', genes = cholinergic_genes)
glia_genes_vln <- create_stacked_vln_plot(seurat_obj, set_ident = 'harmony_clusters_0.3', genes = c(olig_genes, opc_genes, r_glia_genes, mg_genes))
proj_genes_vln <- create_stacked_vln_plot(seurat_obj, set_ident = 'harmony_clusters_0.3', genes = c(cyc_pro_genes, ipc_genes, 'FURIN'))

violin_plt_list <- list(ct_genes_vln, poulio_genes_vln, nowak_genes_vln, camer_genes_vln, fcx_adult_genes_vln,
                        general_genes_vln, exN_genes_vln, InN_genes_vln, glia_genes_vln, proj_genes_vln)

# Cell mapping  -----------------------------------------------------------------------
# Biol_psych  --------------
 for (cam_region in c('fcx', 'ge', 'cer', 'tha', 'hip')) {
  
  message('\nRunning label transfer using Biol Psych region: ', cam_region, '\n') 
  
  seurat_cam <- readRDS(paste0(cameron_dir, 'seurat_', cam_region, '_fetal.rds'))

  # Set anchors based on main clusters
  study.anchors <- FindTransferAnchors(reference = seurat_cam, query = seurat_obj, dims = 1:30)
  predictions <- TransferData(anchorset = study.anchors, refdata = seurat_cam$cellIDs, dims = 1:30)

  cell_type <- paste0('cam_', cam_region, '_pred_id')
  max_score <- paste0('cam_', cam_region, '_max_score')
  predictions <- predictions |>
    dplyr::rename(!!cell_type := predicted.id) |>
    dplyr::rename(!!max_score := prediction.score.max) 
  seurat_obj <- AddMetaData(seurat_obj, metadata = predictions)
  
}

for (cam_region in c('fcx', 'ge', 'cer', 'tha', 'hip')) {
  

  cam_umap <- DimPlot(seurat_obj, reduction = "umap.harmony",
                      label = TRUE, group.by = paste0("cam_", cam_region, "_pred_id")) +
    ggtitle(paste0("Biol Psych ", cam_region, " 2023"))
  cam_score <- FeaturePlot(seurat_obj, 
                               paste0('cam_', cam_region, '_max_score'), 
                               reduction = 'umap.harmony') +
    ggtitle(paste0("Max label transfer score")) 
  
  assign(paste0('cam_', cam_region, '_plt'), cowplot::plot_grid(cam_umap, cam_score))

}

cowplot::plot_grid(cam_fcx_plt, cam_hip_plt, cam_ge_plt, cam_tha_plt, cam_cer_plt)

# cam_umap <- DimPlot(seurat_obj, reduction = "umap.harmony",
#                     label = TRUE, group.by = "cam_fcx_pred_id") +
#   ggtitle(paste0("Biol Psych FCX 2023"))

## Poulioudakis  ------------
load(paste0(public_dir, "sc_dev_cortex_geschwind/raw_counts_mat.rdata"))
gersch_meta <-  as.data.frame(read_csv(paste0(public_dir, "sc_dev_cortex_geschwind/cell_metadata.csv"))) 
gersch_meta <- column_to_rownames(gersch_meta, var = "Cell")
rna_gersch <- CreateSeuratObject(counts = as.matrix(raw_counts_mat), meta.data = gersch_meta)
rna_gersch <- NormalizeData(rna_gersch) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
  ScaleData()

# Set anchors based on main clusters
study.anchors <- FindTransferAnchors(reference = rna_gersch, query = seurat_obj, 
                                     dims = 1:30)
predictions <- TransferData(anchorset = study.anchors, refdata = rna_gersch$Cluster, 
                            dims = 1:30)
predictions <- predictions |>
  dplyr::rename(poulio_predict_id = predicted.id) |>
  select(poulio_predict_id)
seurat_obj<- AddMetaData(seurat_obj, metadata = predictions)

poulio_umap <- DimPlot(seurat_obj, reduction = "umap.harmony",
                       label = TRUE, group.by = "poulio_predict_id") +
  ggtitle("Polioudakis 2018")

harm_umap <- DimPlot(seurat_obj, reduction = "umap.harmony",
                     label = TRUE, group.by = "harmony_clusters_0.3") +
  ggtitle("Harmony - res 0.3")

cowplot::plot_grid(harm_umap, cam_umap, poulio_umap)

## Nowakowski via clustifyr  ----
# Note slight changes for Seurat 5 compatibility: https://github.com/rnabioco/clustifyr/issues/403
load(file = paste0(public_dir, "Nowakowski/Nowakowski2018_cortex_dev_clustifyR.rda"))

# Set varibles genes to consider - 500-1000 recommended 
var_genes <- VariableFeatures(seurat_obj)[1:1000]

res <- clustifyr::clustify(
  input = LayerData(seurat_obj, layer = "data", assay = "RNA"),       # a Seurat object
  ref_mat = ref_cortex_dev,  
  metadata = seurat_obj@meta.data, # matrix of RNA-seq expression data for each cell type
  cluster_col = "harmony_clusters_0.3", # name of column in meta.data containing cell clusters
  obj_out = TRUE,# output Seurat object with cell type inserted as "type" column
  query_genes = var_genes
)

# Plots
nowak_mappings <- cor_to_call(res, cluster_col = "harmony_clusters_0.3")
seurat_obj$nowak_predict_id <- call_to_metadata(
  res = nowak_mappings,                     # data.frame of called cell type for each cluster
  metadata = seurat_obj@meta.data,          # original meta.data table containing cell clusters
  cluster_col = "harmony_clusters_0.3"      # name of column in meta.data containing cell clusters
) |>
  select(type) |>
  dplyr::rename(nowak_predict_id = type) |>
  pull(nowak_predict_id)

nowakowski_umap <- DimPlot(seurat_obj, reduction = "umap.harmony", 
                           label = TRUE, group.by = 'nowak_predict_id') +
  ggtitle("Nowakowski 2017")

## Try subsetting out poor clusters / samples
gene_cnt_plt <- seurat_obj@meta.data |>
  tibble::as_tibble() |>
  dplyr::group_by(sample) |>
  dplyr::tally(nFeature_RNA) |>
  ggplot(aes(x = sample, y = log10(n))) +
  geom_point() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

tscp_cnt_plt <- seurat_obj@meta.data |>
  tibble::as_tibble() |>
  dplyr::group_by(sample) |>
  dplyr::tally(nCount_RNA) |>
  arrange(n) |>
  ggplot(aes(x = sample, y = log10(n))) +
  geom_point() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  
seurat_obj@meta.data |>
  tibble::as_tibble() |>
  dplyr::group_by(sample) |>
  dplyr::distinct(poulio_predict_id, nowak_predict_id) 

seurat_obj$fcx_sample <- seurat_obj@meta.data |>
  mutate(fcx_sample = ifelse(str_detect(sample, "FC"), TRUE, FALSE)) %>%
  pull(fcx_sample)

## Pull out only FCX cells
fcx_only_umap <- FeaturePlot(seurat_obj, 'fcx_sample', reduction = 'umap.harmony') +
  NoLegend() 

#saveRDS(seurat_obj, paste0(results_dir, 'seurat_parse.rds'))
rmarkdown::render(markdown_doc, output_file = markdown_html, output_dir = markdown_dir)

##  ------------------------------------------------------------------------------------
##  ------------------------------------------------------------------------------------




