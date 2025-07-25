#' Wrapper function to create BP Cell h5ad Seurat Object.
#' 
#' Creates a BP cell on disk object for a set of Parse Biosciences
#' h5ad files and returns a single Seurat object.
#' 
#' @param in_dir (vector): of directories where input Parse files are stored.
#'                         anndata.h5ad and cell_metadata.csv required.
#' @param out_dir (string): Output directory path.
#' @param layer_names (vector): strings with unique IDs for each h5ad file. 
#'                              Must be same length as indir.
#' 
#' @returns (Seurat object): Of merged Parse h5ad files.
#' 
#' @examples
#' create_BPCell_h5ad_seurat_object(in_dir = c(plate1_dir, plate2_dir),
#'                                  out_dir = seurat_dir,
#'                                  layer_names = c('Plate1', 'Plate2'))
create_BPCell_h5ad_seurat_object <- function(
    
  in_dir = NULL,
  out_dir = NULL, 
  layer_names = c('plate1', 'plate2')
  
){
  
  # Validate inputs
  if (is.null(in_dir) || is.null(out_dir)) {
    stop("Both 'in_dir' and 'out_dir' must be provided.")
  }
  if (!is.vector(in_dir) || !is.character(in_dir)) {
    stop("'in_dir' must be a character vector of input directories.")
  }
  if (!is.character(out_dir) || length(out_dir) != 1) {
    stop("'out_dir' must be a single string representing the output directory.")
  }
  if (!is.vector(layer_names) || length(layer_names) != length(in_dir)) {
    stop("'layer_names' must be a character vector of the same length as 'in_dir'.")
  }
  
  # Initialise lists
  data_list <- list()
  metadata_list <- list()

  # Loop through h5ad files and output BPCells matrices on-disk
  for (i in seq_along(layer_names)) {
    
    message('Creating BPCell object for ', layer_names[i])
    
    write_matrix_dir(
      mat = open_matrix_anndata_hdf5(paste0(in_dir[i], 'anndata.h5ad')),
      dir = paste0(out_dir, layer_names[i], "_BP")
    )
    
    # Load in BP matrices
    mat <- open_matrix_dir(dir = paste0(out_dir, layer_names[i], "_BP"))
    
    message('Adding ', layer_names[i], ' to colnames ...' )
    colnames(mat) <- paste0(layer_names[i], '_', colnames(mat))
    
    # Get metadata
    #message('Azimuth') # Weird resulst with plate 1 and 2 using Azimuth
    #metadata_list[[i]] <- Azimuth::LoadH5ADobs(paste0(in_dir[i], 'anndata.h5ad'))
    metadata_list[[i]] <- read_csv(paste0(in_dir[i], 'cell_metadata.csv')) |>
      mutate(bc_wells = str_replace(bc_wells, "^", 'plate1_'),
             sample_id = str_replace(sample, 'sample_', ''),
             cell_name = bc_wells) |>
      column_to_rownames('bc_wells')
    
    metadata_list[[i]] <- metadata_list[[i]] |> 
      dplyr::mutate(plate = layer_names[i]) |>
      dplyr::relocate(plate, .before = everything()) 
    data_list[[layer_names[i]]] <- mat
  }
  
  # Name layers
  # message('Adding layer names to matrix list ...')
  # names(data_list) <- layer_names
  
  message('Collapsing metadata list ...')
  metadata <- Reduce(rbind, metadata_list)
  
  message('Creating Seurat Object ...')
  seurat_merged <- CreateSeuratObject(counts = data_list, meta.data = metadata)
  
  message("Removing '-hg38' from rownames ...") 
  rownames(seurat_merged)
  rownames(seurat_merged) <- stringr::str_replace(rownames(seurat_merged), '-hg38', '')
  
  return(seurat_merged)
  
}

#' Load, clean and update Seurat metadata for a sigle plate run.
#' 
#' Removes non cortex samples from seurat object and adds additional 
#' sex and pcw metadata for each sample (if available). 
#' 
#' @param seurat_obj (S4): A seurat object.
#' 
#' @returns (dataframe): A Seurat metadata dataframe.
#' 
#' @examples clean_seurat_meta(seurat_obj = seurat_obj)
#' 
clean_seurat_meta <- function(
    
  seurat_obj = NULL
  
) {
  
  # Load additional metadata files
  sample_fcx_meta <- read_csv(paste0(sheets_dir, 'FC_samples_pcw.csv'))
  sample_meta <- readxl::read_excel(paste0(sheets_dir, 'Fetal single cell eQTL final samples 29-11-24.xlsx'), 
                                    sheet = 'Final cortex samples') |>
    dplyr::rename(sample = ...1) |>
    dplyr::select(sample, PCW, Sex) |>
    mutate(sample = str_replace_all(sample, " \\((A|a)\\)", "_A")) |>
    print(n=Inf)
  
  # ID Claire's samples and add sex and pcw info to meta data
  seurat_meta <- seurat_obj@meta.data |>
    mutate(sample = str_replace_all(sample, "_a", "_A")) |>
    mutate(region = case_when(
      str_detect(sample, pattern = "_WGE") ~ "GE",
      str_detect(sample, pattern = "_Hip") ~ "Hip",
      str_detect(sample, pattern = "_Thal") ~ "Tha",
      str_detect(sample, pattern = "18184") ~ "FC_11pcw", # Rm and not 2nd Trim plate1
      .default = 'CTX'),
      claire_sample = ifelse(str_detect(region, pattern = "GE|Hip|Tha|FC_11pcw"), T, F)) |> 
    left_join(sample_fcx_meta, by = join_by('sample_id' == 'sample')) |>
    left_join(sample_meta, by = join_by('sample_id' == 'sample')) %>%
    mutate(PCW_combined = coalesce(as.character(pcw), as.character(PCW))) %>%
    dplyr::select(-pcw, -PCW) |>
    dplyr::rename(pcw = PCW_combined)
  
  
  return(seurat_meta)
  
}




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


#' Create paired gene and read count QC boxplot for Seurat object. 
#' 
#' Creates a two-panelled set boxplots for the gene and read count 
#' distributions for the samples in a Seurat object metadata.
#' 
#' TODO: Add filter thresholds to params, error handling
#' 
#' @param seurat_obj (S4): A Seurat Object.
#' 
#' @returns (list): A ggplot object.
#' 
#' @examples create_seurat_qc_boxplt(seurat_obj = seurat_obj)
#' 
create_seurat_qc_boxplt <- function(
    
  seurat_obj = NULL,
  add_filt = FALSE
  
) {
  
  if (add_filt == FALSE) {
    
    # tscp per sample
    tscp_boxplot <- seurat_obj@meta.data |>
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
    
    # genes per sample
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
    
  } else {
    
    # Median tscp counts per sample boxplots post filter 1
    tscp_boxplot <- seurat_obj@meta.data |>
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
    
    # Median gene counts per sample boxplots post filter 1
    gene_boxplot <- seurat_obj@meta.data |>
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
    
    
  }
  
  tscp_gene_boxplt <- egg::ggarrange(tscp_boxplot, gene_boxplot, ncol = 2)
  
  return(tscp_gene_boxplt)
  
}

#' Create a list of violin plots for cluster assignment
#' 
#' Create a set of individual viloin plots for a set gene lists. Gene lists
#' need to be pre-loaded.
#' 
#' @param seurat_obj (S4): A Seurat Object.
#' @param meta_col (string): Specify a column is the meta data containing cluster IDs.
#' 
#' @returns (list): A list of violin plots for pre-spscified gene lists.
#' 
#' @examples
#' create_vln_plot_list(seurat_obj = seurat_obj, meta_col = 'seurat_clusters')
create_vln_plot_list <- function(
    
  seurat_obj = NULL,
  meta_col = 'harmony_clusters_0.3'
  
) {
  
  biol_psych_vln_plot <- egg::ggarrange(create_stacked_vln_plot(seurat_obj, meta_col, pfc_features, paste0('PFC ', meta_col)),
                                        create_stacked_vln_plot(seurat_obj, meta_col, hip_features, paste0('Hip ', meta_col)),
                                        create_stacked_vln_plot(seurat_obj, meta_col, wge_features, paste0('GE ', meta_col)),
                                        create_stacked_vln_plot(seurat_obj, meta_col, tha_features, paste0('Tha ', meta_col)),
                                        create_stacked_vln_plot(seurat_obj, meta_col, cer_features, paste0('Cer ', meta_col)))
  
  public_vln_plot <- cowplot::plot_grid(create_stacked_vln_plot(seurat_obj, meta_col, Nowakowski_Fig4A_genes, paste0('Nowak 2017 ', meta_col)),
                                        create_stacked_vln_plot(seurat_obj, meta_col, Pouliodakis_fig1c_genes, paste0('Poulio 2019 ', meta_col)))
  
  fcx_adult_vln <- create_stacked_vln_plot(seurat_obj, set_ident = meta_col, genes = fcx_genes, paste0('FCX adult ', meta_col))
  general_vln <- create_stacked_vln_plot(seurat_obj, set_ident = meta_col, genes = general_genes, paste0('Genral genes ', meta_col))
  neuron_vln <- create_stacked_vln_plot(seurat_obj, set_ident = meta_col, genes = c(exN_genes, inN_genes), paste0('Neurons ', meta_col))
  glia_vln <- create_stacked_vln_plot(seurat_obj, set_ident = meta_col, genes = c(olig_genes, opc_genes, r_glia_genes, mg_genes), paste0('Glia ', meta_col))
  proj_vln <- create_stacked_vln_plot(seurat_obj, set_ident = meta_col, genes = c(cyc_pro_genes, ipc_genes, 'FURIN'), paste0('CycPro-IPC ', meta_col))
  ct_vln <- create_stacked_vln_plot(seurat_obj, set_ident = meta_col, genes = claire_genes, paste0('Claire genes ', meta_col))
  small_pop_vln <- create_stacked_vln_plot(seurat_obj, set_ident = meta_col, genes = small_populations, paste0('Small Populations ', meta_col))
  
  violin_plt_list <- list('biol_psych_vln_plot' = biol_psych_vln_plot, 'public_vln_plot', public_vln_plot, 
                          'fcx_adult_vln' = fcx_adult_vln, 'general_vln' = general_vln, 'neuron_vln' = neuron_vln, 
                          'glia_vln' = glia_vln, 'proj_vln' = proj_vln, 'small_pop_vln' = small_pop_vln, 'ct_vln' = ct_vln)
  
  return(violin_plt_list)
  
  
}


#######
#######
#######

# Single plate functions

#' Load, clean and update Parse metadata for a sigle plate run.
#' 
#' Loads Parse metadata for a single plate and adds pcw and sex 
#' metadata for each sample. 
#' 
#' @param meta_dir (string): Directory where Parse metadata is stored.
#' 
#' @returns (list): 3 metadata objects. 
#' 
#' @examples prep_parse_meta(meta_dir = 'results/plate_1_dir/')
#' 
prep_parse_meta <- function(
    
  meta_dir = NULL
  
) {
  
  # Load additional metadata
  sample_fcx_meta <- read_csv(paste0(sheets_dir, 'FC_samples_pcw.csv'))
  sample_meta <- readxl::read_excel(paste0(sheets_dir, 'Fetal single cell eQTL final samples 29-11-24.xlsx'), 
                                    sheet = 'Final cortex samples') |>
    dplyr::rename(sample = ...1) |>
    dplyr::select(sample, PCW, Sex) |>
    mutate(sample = str_replace_all(sample, " \\((A|a)\\)", "_A")) |>
    print(n=Inf)
  
  # ID Claire's samples and add sex and pcw info to meta data
  parse_meta <- read.csv(paste0(meta_dir, "cell_metadata.csv"), row.names = 1) |>
    mutate(sample = str_replace_all(sample, "_a", "_A"),
           plate = str_extract(meta_dir, 'plate[1-2]')) |>
    mutate(region = case_when(
      str_detect(sample, pattern = "_WGE") ~ "GE",
      str_detect(sample, pattern = "_Hip") ~ "Hip",
      str_detect(sample, pattern = "_Thal") ~ "Tha",
      str_detect(sample, pattern = "18184") ~ "FC_11pcw", # Rm and not 2nd Trim plate1
      .default = 'CTX'),
      sample_id = str_replace(sample, "sample_", ""),
      claire_sample = ifelse(str_detect(region, pattern = "GE|Hip|Tha|FC_11pcw"), T, F)) |> 
    left_join(sample_fcx_meta, by = join_by('sample_id' == 'sample')) |>
    left_join(sample_meta, by = join_by('sample_id' == 'sample')) %>%
    mutate(PCW_combined = coalesce(as.character(pcw), as.character(PCW))) %>%
    dplyr::select(-pcw, -PCW) |>
    dplyr::rename(pcw = PCW_combined)
  
  # Report: make a loop for this
  print(paste0('sample_fcx_meta:'), sep = '\n\n')
  print(sample_fcx_meta |> select(sample) |> arrange(sample) |> distinct() |> pull())
  print(paste0('sample_meta:'), sep = '\n\n')
  print(sample_meta |> select(sample) |> arrange(sample) |> distinct() |> pull())
  print(paste0('parse_meta:'), sep = '\n\n')
  print(parse_meta |> select(sample_id) |> arrange(sample_id) |> distinct() |> pull())
  
  meta_list <- list(
    'parse_meta' = parse_meta,
    'sample_meta' = sample_meta,
    'sample_fcx_meta' = sample_fcx_meta
  )
  
  return(meta_list)
  
}

#' Compare key metadata values across Parse plates 
#' 
#' Compare the raw metadata values from Parse plates and
#' that in the Seurat object. Unsure atm how Seurat deals with
#' reads from the same cells, across two different plates.
#' 
#' @returns (list): Metadata tibbles and common rows between plates.
#' 
compare_cnts_across_plates <- function() {
  
  plate1_meta <- read_csv(paste0(plate1_dir, 'cell_metadata.csv')) |>
    mutate(sample_id = str_replace(sample, "sample_", ""))
  plate2_meta <- read_csv(paste0(plate2_dir, 'cell_metadata.csv')) |>
    mutate(sample_id = str_replace(sample, "sample_", ""))
  all_meta_uniq <- bind_rows(plate1_meta, plate2_meta) 
  seurat_meta <- seurat_obj@meta.data |> as_tibble(rownames = 'rownames')
  
  subset1_wells_tbl <- plate1_meta |> select(bc_wells) |> distinct()
  subset2_wells_tbl <- plate2_meta |> select(bc_wells) |> distinct()
  common_wells_rows <- intersect(subset1_wells_tbl, subset2_wells_tbl) |> pull()
  
  plate1_common_tbl <- plate1_meta |>
    filter(bc_wells %in% (common_wells_rows)) |>
    group_by(sample_id) |>
    summarize(count = n()) |>
    arrange(desc(count))
  
  plate2_common_tbl <- plate2_meta |>
    filter(bc_wells %in% (common_wells_rows)) |>
    group_by(sample_id) |>
    summarize(count = n()) |>
    arrange(desc(count))
  
  
  subset1_samples_tbl <- plate1_meta |> select(sample_id) |> distinct()
  subset2_samples_tbl <- plate2_meta |> select(sample_id) |> distinct()
  common_samples_tbl <- intersect(subset1_samples_tbl, subset2_samples_tbl)
  
  
  cell_cnts_tbl <- tibble('df' = c('plate1_meta', 'plate2_meta', 'plate1&2_meta', 'seurat_meta', 
                                   'common_rows', 'plate1_common', 'plate2_common', 'common_samples'),
                          'cell cnt' = c(nrow(plate1_meta), nrow(plate2_meta), nrow(all_meta_uniq), 
                                         nrow(seurat_meta), length(common_wells_rows), nrow(plate1_common_tbl),
                                         nrow(plate2_common_tbl), nrow(common_samples_tbl))
  )
  
  
  tbl_list <- list('cell_cnts_tbl' = cell_cnts_tbl, 
                   'plate1_meta' = plate1_meta,
                   'plate2_meta' = plate2_meta,
                   'all_meta_uniq' = all_meta_uniq,
                   'seurat_meta' = seurat_meta,
                   'common_rows' = common_rows,
                   'plate1_common_tbl' = plate1_common_tbl,
                   'plate2_common_tbl' = plate2_common_tbl,
                   'common_samples_tbl' = common_samples_tbl)
  
}



