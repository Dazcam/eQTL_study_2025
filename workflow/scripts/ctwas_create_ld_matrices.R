#--------------------------------------------------------------------------------------
#
#    Causal TWAS - Make LD matrices
#
#--------------------------------------------------------------------------------------
#
# Pipeline:  12CTWAS | Rule: create_ld_matrices
#            Upstream:   1000G EUR hg38 PLINK reference (resources/ldsr)
#            Downstream: run_ctwas (LD matrices and SNP info per LD region)
#
# Purpose:   Convert the per-chromosome 1000G EUR PLINK genotypes into one LD matrix
#            (.RDS) and one SNP info file (.Rvar) per LD region, using the LDetect
#            EUR b38 regions shipped with the ctwas package.
#
# Inputs:    ld_ref_dir  1000G.EUR.hg38.{1..22}.bed / .bim
#
# Outputs:   out_dir     {out_name}_chr{chr}.R_snp.{start}_{stop}.RDS / .Rvar
#                        create_ctwas_ld_matrices.done (Snakemake flag)
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
message('\n\nCreating LD matrices and variant info files for causal TWAS ...')

##  Load packages, functions and variables  -------------------------------------------
library(ctwas)

# Input and output paths
ld_ref_dir <- snakemake@params[['ld_ref_dir']]
out_dir <- snakemake@params[['out_dir']]
out_name <- "1000G_EUR_Phase3_hg38"

# For testing locally
# ld_ref_dir <- "../resources/ldsr/ldsr_hg38_refs/plink_files/"
# out_dir    <- "../results/12CTWAS/ld_mats/"

# Print what each variable is set to
message("\nVariables")
message("============================")
print(data.frame(
  variable = c("ld_ref_dir", "out_dir", "out_name"),
  value    = c(ld_ref_dir, out_dir, out_name)), right = FALSE, row.names = FALSE)
message("============================\n")

## Create LD matrices  ------------------------------------------------------------------
# Reference genotype files in PLINK binary genotype format (.bed)
genotype_files <- file.path(ld_ref_dir, paste0("1000G.EUR.hg38.", 1:22, ".bed"))

# Reference variant information files in PLINK .bim format.
varinfo_files <- file.path(ld_ref_dir, paste0("1000G.EUR.hg38.", 1:22,".bim"))

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

message("Done.")

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
