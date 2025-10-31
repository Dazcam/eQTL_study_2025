#--------------------------------------------------------------------------------------
#
#    Causal TWAS - Make LD matrices
#
#--------------------------------------------------------------------------------------

## Info  ------------------------------------------------------------------------------

#  Converts PLINK genotype data to LD matrices and SNP info files
#  Saves LD matrices as .RDS files and SNP info as .Rvar files

##  Load Packages, functions and variables  -------------------------------------------
## Set up logging for smk  ------------------------------------------------------------
if (exists("snakemake")) {
  log_smk <- function() {
    if (exists("snakemake") & length(snakemake@log) != 0) {
      log <- file(snakemake@log[1][[1]], open = "wt")
      sink(log, append = TRUE)
      sink(log, append = TRUE, type = "message")
    }
  }
}
log_smk()
message('\n\nCreating LD matrices and variant info files for causal TWAS ...')

library(ctwas)

ld_ref_dir <- snakemake@params[['ld_ref_dir']]
out_dir <- snakemake@params[['out_dir']]
out_name <- "1000G_EUR_Phase3_hg38"

# Reference genotype files in PLINK binary genotype format (.bed)
genotype_files <- file.path(ldref_dir, paste0("1000G.EUR.hg38.", 1:22, ".bed"))

# Reference variant information files in PLINK .bim format. 
varinfo_files <- file.path(ldref_dir, paste0("1000G.EUR.hg38.", 1:22,".bim"))

# Prepare a data frame region_info for LD regions with columns "chr", "start", and "stop"
# The positions should match those in varinfo_files
region_file <- system.file("extdata/ldetect", "EUR.b38.ldetect.regions.RDS", package = "ctwas")
region_info <- readRDS(region_file)

# Save LD matrices
updated_region_info <- convert_geno_to_LD_matrix(region_info, genotype_files, varinfo_files,
                                                 chrom = 1:22,
                                                 outputdir = out_dir, outname = out_name)

# Create file for smk
file.create(paste0(out_dir, "create_ctwas_ld_matrices.done"))
#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------