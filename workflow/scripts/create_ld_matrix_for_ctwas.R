#--------------------------------------------------------------------------------------
#
#    Creating hg38 LD matrices from individual level genotype data for cTWAS
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

message("\n\nCreating hg38 LD matrices from individual level genotype data for cTWAS ...")

# Load packages
library(tidyverse)
library(ctwas)
library(data.table)

# Functions
download_ldetect_file <- function(filepath = "../resources/ctwas/EUR.b38.ldetect.regions.RDS") {
  # Check if file exists
  if (!file.exists(filepath)) {
    message("cTWAS LD detect regions file not found. Creating directory and downloading...")
    
    dir.create(dirname(filepath), recursive = TRUE, showWarnings = FALSE)
    
    # URL to download from (raw file)
    url <- "https://github.com/xinhe-lab/ctwas/raw/refs/heads/multigroup/inst/extdata/ldetect/EUR.b38.ldetect.regions.RDS"
    download.file(url, destfile = filepath, mode = "wb")
    
    # Check if download was successful
    if (file.exists(filepath)) {
      message("Download complete: ", filepath)
    } else {
      stop("Download failed.")
    }
  } else {
    message("File already exists: ", filepath)
  }
  
  region_file <- readRDS("EUR.b38.ldetect.regions.RDS")
  
  return(region_file)
}

# 
ldref_dir <- '../resources/ldsr/ldsr_hg38_refs/plink_files/'
genotype_files <- paste0(ldref_dir, "1000G.EUR.hg38.", 1:22, ".bed")
varinfo_files <- paste0(ldref_dir, "1000G.EUR.hg38.", 1:22, ".bim")
region_info <- download_ldetect_file()

outputdir <- "./cTWAS_LD"
outname <- "1000G_EUR_Phase3_b38"

region_metatable <- convert_geno_to_LD_matrix(region_info,
                                              genotype_files,
                                              varinfo_files,
                                              chrom = 1:22,
                                              outputdir = outputdir,
                                              outname = outname)

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
