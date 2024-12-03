get_cell_outliers <- function(
    
  seurat_object = NULL,
  mad_thresh = 3,
  mad_range = NULL,
  mito_thresh = 5,
  ribo_thresh = 5,
  umi_column = NULL,
  gene_column = NULL,
  sample_column = NULL
  
) {
  
  # Pull out Mito / Ribo gene names
  # mt_genes <- rownames(sce_obj)[grepl("^MT-", rownames(sce_obj))]
  # ribo_genes <- rownames(sce_obj)[grepl("^RP[LS]", rownames(sce_obj))]
  sample_ids <- seurat_object@meta.data[[sample_column]]
  umi_counts <- seurat_object@meta.data[[umi_column]]
  gene_counts <- seurat_object@meta.data[[gene_column]]
  
  message('Add mito ribo and complexity QC metrics ...')
  seurat_object <- scCustomize::Add_Mito_Ribo(object = seurat_object, species = "Human", overwrite = TRUE)
  seurat_object <- scCustomize::Add_Cell_Complexity(object = seurat_object, species = "Human", overwrite = TRUE)
  
  
  # Need log to avoid negative numbers lower threshold
  message('Calculating outlier thresholds for each sample ...')
  umi_outlier <- scuttle::isOutlier(umi_counts, nmads = mad_thresh, 
                                    type = mad_range, batch = sample_ids, log = TRUE)
  genes_outlier <- scuttle::isOutlier(gene_counts,  nmads = mad_thresh, 
                                      type = mad_range, batch = sample_ids, log = TRUE)
  mito_outlier <- seurat_object$percent_mito > mito_thresh
  ribo_outlier <- seurat_object$percent_ribo > ribo_thresh
  
  cell_outliers <- umi_outlier | genes_outlier | mito_outlier | ribo_outlier 
  
  outlier_cnts_tbl <- tibble(
    measure = c('umi', 'genes', 'mito', 'ribo', 'total'), 
    count = c(sum(umi_outlier), sum(genes_outlier), sum(mito_outlier),
              sum(ribo_outlier), sum(cell_outliers))  
  ) 
  
  seurat_object$umi_outlier <- umi_outlier
  seurat_object$genes_outlier <- genes_outlier
  seurat_object$mito_outlier <- mito_outlier
  seurat_object$ribo_outlier <- ribo_outlier
  seurat_object$cell_outlier <- cell_outliers
  
  message('Cell numbers that will be excluded at specified thresholds:')
  message(paste0(capture.output(outlier_cnts_tbl), collapse = "\n"), '\n')
  
  # Plot outliers
  message('Plotting ...')
  # create_outlier_plots(seurat_object, umi_outlier, genes_outlier, 
  #                      mito_outlier, ribo_outlier, umi_column, genes_column,
  #                      sample_column)
  
  return(seurat_object)
  
}
