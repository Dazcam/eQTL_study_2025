#--------------------------------------------------------------------------------------
#
#    Parse_testing - Basic Seuart processing
#
#--------------------------------------------------------------------------------------

## Info  ------------------------------------------------------------------------------

#  Using Seurat 5 primarily, but also bioconductor packages for QC 

# Load packages  ----------------------------------------------------------------------
library(Seurat)
library(dplyr)
library(Matrix)
library(ggplot2)
library(scater)
library(scDblFinder)
library(tidyverse)
library(harmony)
library(clustifyr)

# Set variables  ----------------------------------------------------------------------
root_dir <- '/Users/darren/Desktop/parse_testing/'
config_dir <- paste0(root_dir, 'config/')
results_dir <- paste0(root_dir, 'results/01PARSE/')
combine_dir <- paste0(results_dir, 'combine/')
mat_dir <- paste0(combine_dir, 'all-sample/DGE_filtered/')
report_dir <- paste0(combine_dir, 'all-sample/report/')
script_dir <- paste0(root_dir, 'workflow/scripts/')
markdown_doc <- paste0(script_dir, 'parse_testing_qc_and_seurat.Rmd')
markdown_html <- paste0(results_dir, 'parse_testing.html')
public_dir <- paste0('~/Desktop/fetal_brain_snRNAseq_110122/resources/public_datasets/')
cameron_dir <- paste0(root_dir, 'resources/public_data/cameron_2023/')

resolution_set <- c(0.3, 0.5, 0.8)

# Set features
pfc_features <- c("GAD1", "SLC17A7", "EOMES", "GLI3", "OLIG1", 
                  "MKI67", "C3", "ITM2A", "SST", "CALB2", 
                  "SCGN", "TLE3", "FEZF2", "CRYM", "LHX2")
wge_features <- c("GAD1", "SLC17A7", "EOMES", "GLI3", "OLIG1", 
                  "MKI67", "C3", "ITM2A", "LHX6", "SIX3", 
                  "PROX1", "TSHZ1", "DLX1", "SCGN")
hip_features <- c("NEUROD1", "GRIK4", "EOMES", "GLI3", "OLIG1", 
                  "MKI67", "C3", "ITM2A", "SLC17A6", "ADARB2",
                  "GAD2", "TNC", "PROX1", "RELN", "LHX6")
tha_features <- c("EOMES", "GLI3", "OLIG1", "MKI67", "C3", 
                  "ITM2A", "SLC1A6", "LHX9", "TNC", "GAD2", 
                  "PAX6", "SLC17A6")
cer_features <- c("GAD1", "EOMES", "GLI3", "OLIG1", "MKI67", 
                  "C3", "ITM2A", "CA8", "ITPR1", "RBFOX3", 
                  "RELN")
claire_genes <- c(
  "PDGFRA",  #OPC
  "C3", #Microglia
  "GAD1", "GAD2", #IN
  "RELN", #Cajal-Retzius cell
  "SLC1A3", #Astrocytes
  "TNC", "GLI3", #RG -> putative glioblasts
  "FN1", "FLT1", #Endothelia
  "PAX6", #NPC
  "MKI67", #IP
  "RGS5", #Pericytes
  "TJP1", #Neuroepithelial cells
  "TLE3", "LHX2", #upper cortical layers
  "PROX1", #CGE
  'BCL11B', # Sub-cortical projecting maturing Ns (Nowak 2017)
  'SATB2', # Intercortically projecting maturing Ns (Nowak 2017)
  'CPNE8', # Upper Layers ExNs (Nowak 2017)
  'HCRTR2', # Subplate ExNs (Nowak 2017)
  'KCNJ6', # Maturing ExNs (Nowak 2017)
  "RBFOX3")

Nowakowski_Fig4A_genes <- c(
  'C1QC', # MG 
  'TTR', # Choroid Plexus
  'PECAM1', # Endothelial Cells
  'TBX18', # Mural Cells inc. pericytes
  'OMG', # OPCs
  'CSPG5', # Astrocytes 
  'GHR', # MGE-RG1 
  'FZD8', # Pallial radial glia (MGE-RG1)
  'HJURP', # Dividing cells (RG-div1, IPC-div1, IPC-div2, MGE-div, MGE-IPC4, MGE-IPC45),
  'CRYAB', # Truncated radial glia in G1 (tRG)
  'HEPACAM', # Outer radial glia and astrocytes
  'LMO1', # MGE progenitors (MGE-RG1-2, MGE-div, MGE-IPC1-3)
  'LEF1', # Early progenitor
  'VEPH1', # RG-like clusters among MGE progenitors (MGE-RG1-2) and CTX RG (oRG)
  'NEUROG1', # IPC
  'NRN1', # IPC
  'BEST3',
  'PIF1',
  'ASCL1',
  'EOMES',
  'CALB2',
  'SP8',
  'SST',
  'LINCO1305',
  'TAC3',
  'LHX6',
  'NDST4',
  'SLN',
  'RSPO3',
  'KCNJ6',
  'NEFM',
  'NEFL',
  'HCRTR2',
  'CRYM',
  'TRPM3',
  'PKD1',
  'UNC5D',
  'NHLH1'
)

Pouliodakis_fig1c_genes <- c(
  'RGS5', # Pericyte
  'CCL3', 'AIF1', # MG
  'ITM2A', 'CLDN5', 'ESAM', # Endo
  'SOX5', # ExN deep layer
  'DLX1', 'DLX2', 'LHX6', 'DLX5', # InNs
  'STMN2', 'NEUROD6', 'SATB2', # Maturning Ns upper layer
  'PPP1R17', 'SSTR2', 'EOMES', 'PENK', # IP
  'SATB2', 'STMN2', 'NEUROD6', # ExN maturing
  'NEUROD6', 'POU2F2', # Migrating ExN
  'HMGB2', 'SOX2', 'MKI67', 'PCNA', # Cyclic Progenitors
  'PTPRZ1', 'OLIG1', 'PDGFRA', # OPCs
  'VIM', 'PTPRZ1', 'SOX2', 'SLC1A3', 'HES1', 'HOPX'
)
# Load functions ----------------------------------------------------------------------
source(paste0(root_dir, 'workflow/scripts/parse_testing_functions.R'))
options(ggrepel.max.overlaps = Inf)

## Parse QC - pre-Seurat  -------------------------------------------------------------
# Agg table
agg_samp_tbl <- read_csv(paste0(combine_dir, 'agg_samp_ana_summary.csv')) 
agg_samp_tbl <- as_tibble(cbind(nms = names(agg_samp_tbl), t(agg_samp_tbl))) 
agg_samp_tbl <- agg_samp_tbl |> 
  janitor::row_to_names(row_number = 1) |>
  dplyr::rename(sample_id = statistic) |>
  mutate_at(vars(-sample_id), as.numeric) |>
  
  select(!ends_with("Q30") & !ends_with("valid") & !ends_with("fraction")) |>
  select(!ends_with("read1"))

boxplt_lst <- create_qc_boxplots(agg_samp_tbl, 'sample_id')
cowplot::plot_grid(plotlist = boxplt_lst)

tscp_boxplot <- seurat_obj@meta.data |>
  as_tibble() |>
  mutate(sample_id = str_replace(sample, "sample_", "")) |>
  ggplot(aes(x = sample_id, y = log10(nCount_RNA))) +
  geom_boxplot(outlier.size = 0.1, outlier.alpha = 0.1) +
  theme_linedraw() +
  theme(axis.text.x = element_blank()) +
  xlab(NULL) 

gene_boxplot <- seurat_obj@meta.data |>
  as_tibble() |>
  mutate(sample_id = str_replace(sample, "sample_", "")) |>
  ggplot(aes(x = sample_id, y = nFeature_RNA)) +
  geom_boxplot(outlier.size = 0.1, outlier.alpha = 0.1) +
  theme_linedraw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  xlab(NULL) 

tscp_gene_boxplt <- egg::ggarrange(tscp_boxplot, gene_boxplot)

median_cnt_tbl <- seurat_obj@meta.data |> 
  tibble::as_tibble() |> 
  dplyr::group_by(sample) |> 
  summarise(genes_per_cell_median = median(nFeature_RNA),
            tscp_per_cell_median = median(nCount_RNA)) |>
  arrange(genes_per_cell_median)

gene_by_cells_plot <- ggplot(agg_samp_tbl, aes(x = log10(number_of_cells), y = hg38_median_genes_per_cell)) + 
  geom_point() +
  ggrepel::geom_text_repel(aes(label = sample_id), na.rm = TRUE) +
  theme_minimal()   

# Wells round 1
sample_key_tbl <- read_delim(paste0(config_dir, 'sample-list.txt'), col_names = FALSE) |>
  dplyr::rename(sample = X1, well = X2)
sample_loading_tbl <- read_tsv(paste0(config_dir, 'sample_loading.tsv'), col_names = FALSE) |>
  dplyr::rename(well = X1, loadings = X2)
counts_rnd1_tbl <- read_csv(paste0(report_dir, 'cell_counts_by_rnd1_well.csv')) |>
  dplyr::rename(cell_cnt = count)
tscp_median_rnd1_tbl <- read_csv(paste0(report_dir, 'tscp_median_by_rnd1_well.csv')) |>
  dplyr::rename(tscp_med_cnt = count)

rnd1_tbl <- sample_key_tbl |>
  inner_join(counts_rnd1_tbl) |>
  inner_join(tscp_median_rnd1_tbl) |>
  inner_join(sample_loading_tbl) |>
  mutate(Row = gsub("[^a-zA-Z]", "", well), .before=2) |>
  mutate(Column = parse_number(well), .before=2) 

# Plot heatmaps
rnd1_hmap <- plot_parse_rnd1_heatmaps(rnd1_tbl)

# Pull out outliers if necessary
rnd1_tbl |>
  filter(sample %in% c(14402, 14635, 11654, 17109))

## Seurat Processing ------------------------------------------------------------------
mat <- ReadParseBio(mat_dir)
meta <- read.csv(paste0(mat_dir, "cell_metadata.csv"), row.names = 1)
seurat_obj <- CreateSeuratObject(mat, meta.data = meta)
counts_preFilt_tbl <- tibble('Cells' = c(ncol(seurat_obj)),
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
sce_obj <- scDblFinder(sce_obj, dbr = 0.0003, samples = sce_obj$sample)
table(sce_obj$scDblFinder.class)
sce_obj$doublet_logic <- ifelse(sce_obj$scDblFinder.class == "doublet", TRUE, FALSE)
all_outliers <- sce_obj@colData |>
  as_tibble() |>
  mutate(all_outlier = as.logical(pmax(cell_outlier, doublet_logic))) |>
  select(all_outlier) %>%
  pull

outlier_cnts_tbl <- tibble(
  measure = c('umi', 'genes', 'mito', 'ribo', 'doublets', 'total'), 
  count = c(sum(seurat_obj$umi_outlier), sum(seurat_obj$genes_outlier), 
            sum(seurat_obj$mito_outlier),sum(seurat_obj$ribo_outlier), 
            sum(seurat_obj$doublet_logic), sum(seurat_obj$all_outliers))  
) 

# Subset Seurat object
seurat_obj <- subset_seurat_object(seurat_obj, all_outliers)
counts_postFilt_tbl <- tibble('Cells' = c(ncol(seurat_obj)),
                             'Genes' = c(nrow(seurat_obj)))

# Run basic pipeline
seurat_obj <- run_seurat_process(seurat_obj)

# Post pipeline QCs
# Check PCs
pca_plot <- DimHeatmap(seurat_obj, dims = 20:40, cells = 500, balanced = TRUE)

# Create resolution plots for comparison
res_plotlist <- create_resolution_plotlist(seurat_obj, resolution = resolution_set,
                                           meta_id = 'sample', gene_set = claire_genes)



# Plot QCs
cluster_qc_plot <- create_cluster_qc_plot(seurat_obj, 30, 'sample', 'RNA_snn_res.0.3')
meta_col_cluster_cnts <- get_meta_col_counts(seurat_obj, 'RNA_snn_res.0.3')

## Integration - using Harmony
## Is this necessary? 
## Integration - using Harmony on sploit object / maybe not necessary
## Crashing near the end of process (maybe need to split before noraml Seurat process)
# seurat_obj[["RNA"]] <- split(seurat_obj[["RNA"]], f = seurat_obj$sample)
# seurat_object <- run_integration(seurat_obj, 'harmony', 30, resolution_set)
# seurat_obj <- JoinLayers(seurat_obj)
# ####

seurat_obj <- run_harmony(seurat_obj, 'harmony', 'sample', 30, resolution_set)

#save.image(paste0(results_dir, 'parse_tetsing.RData'))
#load(paste0(results_dir, 'parse_tetsing.RData'))

integration_plotlist <- create_integration_plotlist(seurat_obj, 
                                                    meta_id = 'sample', 
                                                    dims = 30,
                                                    reduction = resolution_set) 
  
## Cell type assignment  --------------------------------------------------------------
# Gene lists  -----
biol_psych_vln_plot <- egg::ggarrange(create_stacked_vln_plot(seurat_obj, 'harmony_clusters_0.3', pfc_features, 'PFC harmony 0.3'),
                                  create_stacked_vln_plot(seurat_obj, 'harmony_clusters_0.3', hip_features, 'Hip harmony 0.3'),
                                  create_stacked_vln_plot(seurat_obj, 'harmony_clusters_0.3', wge_features, 'GE harmony 0.3'),
                                  create_stacked_vln_plot(seurat_obj, 'harmony_clusters_0.3', tha_features, 'Tha harmony 0.3'),
                                  create_stacked_vln_plot(seurat_obj, 'harmony_clusters_0.3', cer_features, 'Cer harmony 0.3'))

public_vln_plot <- cowplot::plot_grid(create_stacked_vln_plot(seurat_obj, 'harmony_clusters_0.3', Nowakowski_Fig4A_genes, 'Nowak 2017 harmony 0.3'),
                                      create_stacked_vln_plot(seurat_obj, 'harmony_clusters_0.3', Pouliodakis_fig1c_genes, 'Poulio 2019 harmony 0.3'))


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

cam_umap <- DimPlot(seurat_obj, reduction = "umap.harmony",
                    label = TRUE, group.by = "cam_fcx_pred_id") +
  ggtitle(paste0("Biol Psych FCX 2023"))

## Poulioudakis  ------------
load(paste0(public_dir, "sc_dev_cortex_geschwind/raw_counts_mat.rdata"))
gersch_meta <-  as.data.frame(read_csv(public_dir, "sc_dev_cortex_geschwind/cell_metadata.csv")) 
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
  metadata = seurat_obj, # matrix of RNA-seq expression data for each cell type
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
  arrange(n)
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


 
##  ------------------------------------------------------------------------------------
##  ------------------------------------------------------------------------------------

rmarkdown::render(markdown_doc, output_file = markdown_html, output_dir = results_dir)



