#--------------------------------------------------------------------------------------
#
#    Parse_testing - Basic Seuart processing
#
#--------------------------------------------------------------------------------------

## Info  ------------------------------------------------------------------------------

#  Using Seurat 5 primarily, but also bioconductor packages for QC 

##  Load Packages, functions and variables  -------------------------------------------
message('Setting environment variables ...')
if (Sys.info()[["nodename"]] == "Darrens-iMac-2.local") {
  
  root_dir <- '~/Desktop/eQTL_study_2025/'
  
  source(paste0(root_dir, 'workflow/scripts/eQTL_2025_functions.R'))
  source(paste0(root_dir, 'workflow/scripts/eQTL_2025_Renv.R'))
  source(paste0(root_dir, 'workflow/scripts/eQTL_2025_gene_lists.R'))
  
} else {
  
  source('scripts/eQTL_2025_functions.R')
  source('scripts/eQTL_2025_Renvs.R')
  source('scripts/eQTL_2025_gene_lists.R')
  
}

# Load functions ----------------------------------------------------------------------
source('~/Desktop/parse_testing/workflow/scripts/parse_testing_functions.R')

## Load data  -------------------------------------------------------------------------
# mat <- ReadParseBio(mat_dir) # Crashes locally with 1M cells
mat <- ReadParseBio(plate1_dir)
metadata_lst <- prep_parse_meta(meta_dir = plate1_dir)
#metadata_lst[['parse_meta']] |> column_to_rownames(sample)

seurat_obj <- CreateSeuratObject(mat, meta.data = metadata_lst[['parse_meta']])
metadata_lst[['seurat_meta']] <- seurat_obj@meta.data |> head(12)
rm(mat, parse_meta, sample_meta, sample_fcx_meta)

## QC ----------------------------------------------------------------------------------
counts_preFilt_tbl <- tibble('Stage' = 'Pre-filter', 
                             'Cells' = c(ncol(seurat_obj)),
                             'Genes' = c(nrow(seurat_obj)))

# Remove Claire's samples
seurat_obj <- subset(seurat_obj, subset = claire_sample == FALSE)
counts_FiltClaire_tbl <- tibble('Stage' = 'CT-filt', 
                                'Cells' = c(ncol(seurat_obj)),
                                'Genes' = c(nrow(seurat_obj)))

# Create QC boxplot
tscp_gene_boxplt <- create_seurat_qc_boxplt(seurat_obj)

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
# Cells that will be removed by filter
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

# Cell counts retained per sample after filter 1
cell_cnts_post_filt_tbl <- seurat_obj@meta.data |> 
  as_tibble() |>
  arrange(nCount_RNA) |>
  filter(nCount_RNA >= 500 & nFeature_RNA >= 300) |>
  group_by(sample_id) |>
  summarize(count = n()) |>
  arrange(count) |>
  print(n = Inf) |>
  flextable() |> 
  set_caption(caption = "Counts per sample after filter 1") |>
  theme_vanilla()

# Median counts retained after filter
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
  set_caption(caption = "Median counts after filter 1") |>
  theme_vanilla()


tscp_gene_post_filt_boxplt <- create_seurat_qc_boxplt(seurat_obj = seurat_obj, 
                                                      add_filt = TRUE)

# Subset based on filter 1 
seurat_obj <- subset(seurat_obj, subset = nFeature_RNA >= 300 & nCount_RNA >= 500)
counts_subset1_tbl <- tibble('Stage' = 'Filt 1', 
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

# Subset Seurat object: filter 2
seurat_obj <- subset_seurat_object(seurat_obj, sce_obj$all_outliers)
counts_subset2_tbl <- tibble('Stage' = 'Filt 2', 
                             'Cells' = c(ncol(seurat_obj)),
                             'Genes' = c(nrow(seurat_obj)))

# Final counts after QC
filt_all_counts_tbl <- rbind(counts_preFilt_tbl, 
                             counts_FiltClaire_tbl, 
                             counts_subset1_tbl, 
                             counts_subset2_tbl) |>
  flextable() |> 
  set_caption(caption = "Counts after filters") |>
  theme_vanilla()
rm(sce_obj)

# Cell counts after filter 2
cell_cnts_final_tbl <- seurat_obj@meta.data |> 
  as_tibble() |>
  arrange(nCount_RNA) |>
  group_by(sample_id) |>
  summarize(count = n()) |>
  arrange(count) |>
  print(n = Inf) |>
  flextable() |> 
  set_caption(caption = "Counts per sample after filter 2") |>
  theme_vanilla()

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


# Boxplots post filter 2
tscp_gene_final_boxplt <- create_seurat_qc_boxplt(seurat_obj)

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
fcx_umap_split_plt <- DimPlot(seurat_obj, reduction = "umap.harmony", group.by = "harmony_clusters_0.3",
                              split.by = "harmony_clusters_0.3", ncol = 4, cols = discreet_cols_n23)
fcx_feat_plt <- FeaturePlot(seurat_obj, 'fcx_sample', reduction = 'umap.harmony') 
cowplot::plot_grid(fcx_umap_plt, fcx_feat_plt)
  
fcx_fact_barplt <- create_proportion_barplot(seurat_obj, 'sample_id_fcx_fact', meta_id = 'harmony_clusters_0.3') +
  scale_fill_manual(values = discreet_cols_n23) +
  theme(axis.text.x = element_text(angle = 90)) 

fcx_pcw_barplt <- create_proportion_barplot(seurat_obj, 'pcw', meta_id = 'harmony_clusters_0.3') +
  scale_fill_manual(values = discreet_cols_n23) +
  theme(axis.text.x = element_text(angle = 90)) 

fcx_in_clust_tbl <- seurat_obj@meta.data |>
  select(harmony_clusters_0.3, sample_id, fcx_sample) |>
  group_by(harmony_clusters_0.3, fcx_sample) |>
  summarize(count = n()) |>
  mutate(percentage = count / sum(count) * 100) |>
  arrange(as.numeric(as.character(harmony_clusters_0.3))) |>
  print(n = Inf) |>
  flextable() |> 
  set_caption(caption = "Harmony Clusters: FCX sample proportion") |>
  theme_vanilla()

## Cell type assignment  --------------------------------------------------------------
# Violin plots  -----
vln_plt_lst <- create_vln_plot_list(seurat_obj, 'harmony_clusters_0.3')

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
cam_fcx_umap <- DimPlot(seurat_obj, reduction = "umap.harmony",
                    label = TRUE, group.by = paste0("cam_fcx_pred_id")) +
  ggtitle(paste0("Biol Psych fcx 2023"))

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
seurat_obj <- AddMetaData(seurat_obj, metadata = predictions)

poulio_umap <- DimPlot(seurat_obj, reduction = "umap.harmony",
                       label = TRUE, group.by = "poulio_predict_id") +
  ggtitle("Polioudakis 2018")

harm_umap <- DimPlot(seurat_obj, reduction = "umap.harmony",
                     label = TRUE, group.by = "harmony_clusters_0.3") +
  ggtitle("Harmony - res 0.3")



## Nowakowski via clustifyr  ----
# Note slight changes for Seurat 5 compatibility: https://github.com/rnabioco/clustifyr/issues/403
load(file = paste0(public_dir, "Nowakowski/Nowakowski2018_cortex_dev_clustifyR.rda"))

# Set variables genes to consider - 500-1000 recommended 
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

cowplot::plot_grid(harm_umap, poulio_umap, nowakowski_umap)

#saveRDS(seurat_obj, paste0(results_dir, 'seurat_parse.rds'))

rmarkdown::render(markdown_doc, output_file = markdown_html, output_dir = markdown_dir)

##  ------------------------------------------------------------------------------------
##  ------------------------------------------------------------------------------------




