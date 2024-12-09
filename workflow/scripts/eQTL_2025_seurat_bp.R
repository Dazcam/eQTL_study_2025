#--------------------------------------------------------------------------------------
#
#    Seurat processing - all Parse plates
#
#--------------------------------------------------------------------------------------

## Info  ------------------------------------------------------------------------------

#  Plate 1: 861,960 cells
#  Plate 2: 719,652 cells
#  Plate 1&2: 1,537,666

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
#source('~/Desktop/parse_testing/workflow/scripts/eQTL_2025_functions.R')
#source('~/Desktop/eQTL_study_2025/workflow/scripts/eQTL_2025_functions.R')
options(ggrepel.max.overlaps = Inf)
options(future.globals.maxSize = 1e9)

## Load data  -------------------------------------------------------------------------
# mat <- ReadParseBio(mat_dir) # Crashes locally with 1M cells
seurat_obj <- create_BPCell_h5ad_seurat_object(in_dir = c(plate1_dir, plate2_dir),
                                               out_dir = seurat_dir,
                                               layer_names = c('Plate1', 'Plate2'))

# Compare Metadata values across plates
initial_cnts_lst <- compare_cnts_across_plates()

options(knitr.duplicate.label = "allow")
rmarkdown::render(markdown_bp_doc, output_file = markdown_bp_html, output_dir = markdown_dir)



test_df <- tibble(
  harmony_clusters_0.3 = factor(c(0, 0, 1, 1, 2), levels = 0:4),
  fcx_sample = c(TRUE, FALSE, TRUE, FALSE, TRUE),
  sample_id = paste0("sample", 1:5),
  cell_count = c(100, 200, 150, 250, 300),
  percentage = c(33.33, 66.67, 37.50, 62.50, 100.00)
)

test_list[[1]] <- test_df