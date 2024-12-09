#--------------------------------------------------------------------------------------
#
#    R environments 
#
#--------------------------------------------------------------------------------------

## Info  ------------------------------------------------------------------------------

# Load R environment - options are local or cluster

##  Load Packages  --------------------------------------------------------------------
message('Loading packages ...')
# Load packages  ----------------------------------------------------------------------
library(BPCells)
library(Seurat)
library(SeuratObject)
library(SeuratDisk)
library(Azimuth)
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

## Set variables  ---------------------------------------------------------------------
if (exists("snakemake")) { 
  log_smk() 
  root_dir <- snakemake@params[['root_dir']]
  region <- snakemake@params[['region']]
  threads <- snakemake@threads

# Only works atm when snakemake@params[['wgcna']] is present in rule  
#  if (exists(snakemake@params[['wgcna']])) {
#    enableWGCNAThreads(snakemake@threads)
#    cat("WGCNA threads enabled. Set to: ", threads, '\n')
#  } else {
#    cat("Threads are set to: ", threads, '\n')
#    future::plan("multicore", workers = snakemake@threads)
#    plan()
#  }

}

options(ggrepel.max.overlaps = Inf)
options(future.globals.maxSize = 1e9)

config_dir <- paste0(root_dir, 'config/')
results_dir <- paste0(root_dir, 'results/02PARSE/')
seurat_dir <- paste0(root_dir, 'results/03SEURAT/')
resources_dir <- paste0(root_dir, 'resources/')
sheets_dir <- paste0(resources_dir, 'sheets/')
bp_cells_dir <- paste0(seurat_dir, 'bp_cells_dir/')
plate1_dir <- paste0(results_dir, 'combine_plate1/all-sample/DGE_filtered/')
plate2_dir <- paste0(results_dir, 'combine_plate2/all-sample/DGE_filtered/')
#report_dir <- paste0(combine_dir, 'all-sample/report/')
script_dir <- paste0(root_dir, 'workflow/scripts/')
markdown_doc <- paste0(script_dir, 'eQTL_2025_seurat.Rmd')
markdown_html <- paste0(results_dir, 'eQTL_2025_seurat.html')
markdown_bp_doc <- paste0(script_dir, 'eQTL_2025_seurat_bp.Rmd')
markdown_bp_html <- paste0(results_dir, 'eQTL_2025_seurat_bp.html')
markdown_dir <- paste0(root_dir, 'workflow/reports/markdown/')
public_dir <- paste0(resources_dir, 'public_datasets/')
cameron_dir <- paste0(public_dir, 'cameron_2023/')


resolution_set <- c(0.3, 0.5, 0.8)

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
