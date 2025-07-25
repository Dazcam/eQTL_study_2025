#--------------------------------------------------------------------------------------
#
#    eQTL study 2025 - Initial QC
#
#--------------------------------------------------------------------------------------

## Info  ------------------------------------------------------------------------------

#  When attempting to combine the data and run QC on both plate at once various errors 
#  were encountered. 

# + Using bp cells:
# + Using manual merging: Error in cbind.Matrix(x, y, deparse.level = 0L) : 
#                         p[length(p)] cannot exceed 2^31-1
# Running QC and then merging

# The latter is a hard limit the inherent limitation on object size in R so decided to 
# run QC on individual plates separately to reduce the object sizes before merging.

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

# QC plates separately  ----------------------------------------------------------------
for (plate_dir in c(plate1_dir, plate2_dir)) {
  
  plate <- str_extract(plate_dir, "plate[1-2]")
  
  ## Load data  -----------------------------------------------------------------------
  message("Loading ", plate, " matrix and metadata ...")
  mat <- ReadParseBio(plate_dir)
  metadata_lst <- prep_parse_meta(meta_dir = plate_dir)

  message("Loading ", plate, " matrix and metadata ...")
  seurat_obj <- CreateSeuratObject(mat, meta.data = metadata_lst[['parse_meta']])
  seurat_obj <- RenameCells(seurat_obj, add.cell.id = paste0(plate, '_'))
  metadata_lst[['seurat_meta']] <- seurat_obj@meta.data |> head(12)
  rm(mat)

  ## QC -------------------------------------------------------------------------------
  counts_preFilt_tbl <- tibble('Stage' = 'Pre-filter', 
                               'Cells' = c(ncol(seurat_obj)),
                               'Genes' = c(nrow(seurat_obj)))
  
  # Remove Claire's samples
  message('Removing non-cortex samples ...')
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
  message('Applying filter 1 ...')
  seurat_obj <- subset(seurat_obj, subset = nFeature_RNA >= 300 & nCount_RNA >= 500)
  counts_subset1_tbl <- tibble('Stage' = 'Filt 1', 
                               'Cells' = c(ncol(seurat_obj)),
                               'Genes' = c(nrow(seurat_obj)))
  
  # Convert to single cell experiment
  message('Creating sce object ...')
  sce_obj <- SingleCellExperiment(list(counts = as(seurat_obj[["RNA"]]$counts, "dgCMatrix")),
                                  colData = seurat_obj@meta.data)
  
  # Identify cell outliers and apply filters to sce object
  message('Identifying outliers ...')
  sce_obj <- get_cell_outliers(sce_obj, 3, 'higher', 5, 5, 
                               'nCount_RNA', 'nFeature_RNA', 'sample')
  outlier_plt <- create_outlier_plots(sce_obj, sce_obj$umi_outlier, sce_obj$genes_outlier, sce_obj$mito_outlier,
                                      sce_obj$ribo_outlier, 'nCount_RNA', 'nFeature_RNA', 'sample')
  
  # Identify doublets
  message('Identifying doublets ...')
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
  message('Applying filter 2 ...')
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
  
  # Set plate specific markdown html name
  markdown_plate_html <- str_replace(markdown_qc_html, '.html', paste0('_', plate, '.html'))
  
  # Run markdown
  message('Generating report ...')
  rmarkdown::render(markdown_qc_doc, output_file = markdown_plate_html, output_dir = markdown_dir)
  
  # Assign plate specific Seurat object
  message('Assigning plate specific Seurat object ...')
  assign(paste0('seurat_obj_', plate), seurat_obj)
  
  # Clean
  rm(seurat_obj, metadata_lst, outlier_plt)
  rm(list=ls(pattern="cnt"))
  rm(list=ls(pattern="tscp"))
  rm(list=ls(pattern="counts"))

}

merged_obj <- merge(x = seurat_obj_plate1, y = seurat_obj_plate2)

data_list <- list()

write_matrix_dir(mat = merged_obj[["RNA"]]$counts.1, dir = paste0(seurat_dir, 'ctx_plate1_BP'))
write_matrix_dir(mat = merged_obj[["RNA"]]$counts.2, dir = paste0(seurat_dir, 'ctx_plate2_BP'))

data_list[['ctx_plate1_BP']] <- open_matrix_dir(dir = paste0(seurat_dir, 'ctx_plate1_BP'))
data_list[['ctx_plate2_BP']] <- open_matrix_dir(dir = paste0(seurat_dir, 'ctx_plate2_BP'))


seurat_obj <- CreateSeuratObject(data_list, merged_obj@meta.data)

mat <- open_matrix_dir(dir = paste0(out_dir, layer_names[i], "_BP"))


rm(seurat_obj_plate1, seurat_obj_plate2)


