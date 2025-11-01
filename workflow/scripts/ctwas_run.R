#--------------------------------------------------------------------------------------
#
#    Causal TWAS
#
#--------------------------------------------------------------------------------------

## Info  ------------------------------------------------------------------------------

#  Prep input bed and covarite files for FastQTL 
#  Running locally for now

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
message('\n\nRunning causal TWAS ...')

# Load packages
#library(tidyverse)
#remotes::install_github("xinhe-lab/ctwas",ref = "main")
#BiocManager::install("VariantAnnotation")
#remotes::install_github("mrcieu/gwasvcf")
library(ctwas)
#library(VariantAnnotation)
#library(gwasvcf)
#library(RSQLite)

# Input and output paths
gwas <- snakemake@input[['gwas']] # From 07prep_gwas.smk
weights_dir <- snakemake@params[['weights_dir']] # Folder containing .wgt.RDat files
ld_dir <- snakemake@params[['ld_dir']]
#snp_info <- snakemake@input[['snp_info']]
cell_type <- snakemake@wildcards[['cell_type']]
output <- snakemake@output

# Read in data
message("\nLoading data ...\n")
message("GWAS loaded from: ", gwas)
message("weights_dir loaded from: ", weights_dir)
message("cell_type: ", cell_type)
message("Output will be saved to: ", output)

# Old
#ld_dir <- "../results/12CTWAS/ld_mats/"  # Downloaded LD files
#ref_snp_info_file <- "/path/to/ukb_b38_0.1_allchrs_var_info.Rvar.gz"  # Full ref variant info


# Prep GWAS
message('Loading GWAS ...')
z_snp_raw <- read_tsv(gwas) 
z_snp <- z_snp_raw |>
  select(id = SNP, A1, A2, z = Z)
gwas_n <- as.numeric(names(sort(table(z_snp_raw$N), decreasing = TRUE)[1]))

# Region info: included with package
# local - ctwas_0.1.38
#region_file <- system.file("extdata/ldetect", "EUR.b38.bed", package = "ctwas")
#region_info <- read_tsv(region_file)

# remote - ctwas_0.5.34
message('Loading region info from extdata ...')
region_file <- system.file("extdata/ldetect", "EUR.b38.ldetect.regions.RDS", package = "ctwas")
region_info <- readRDS(region_file)

# Load ref SNP info (for mapping) - think this is pnly used when not using LD matrices
#ref_snp_info <- data.table::fread(ref_snp_info_file, sep = "\t")
#class(ref_snp_info) <- "data.frame"

# Create metadata table for LD files (one row per region)
message('Creating metadata table for LD files ...')
region_metatable <- region_info
LD_filestem <- sprintf("1000G_EUR_Phase3_hg38_chr%s.R_snp.%s_%s", 
                       region_metatable$chrom, 
                       region_metatable$start, 
                       region_metatable$stop)
region_metatable$LD_file <- file.path(ld_dir, paste0(LD_filestem, ".RDS"))
region_metatable$SNP_file <- file.path(ld_dir, paste0(LD_filestem, ".Rvar"))

# Generate snp_map and LD_map
message('Generating snp_map and LD_map for cTWAS ...')
res <- create_snp_LD_map(region_metatable, ncore = 4)
snp_map <- res$snp_map
LD_map <- res$LD_map  # Not directly used in preprocess_weights, but for ctwas() later

# Optional: Harmonize GWAS first (add after reading z_snp)
# message('Preprocessing GWAS for cTWAS ...')
# z_snp <- preprocess_z_snp(z_snp, snp_map, 
#                           drop_multiallelic = TRUE, 
#                           drop_strand_ambig = TRUE)  # Converts to chr:pos_ref_alt

# FUSION weights
message('Preprocessing FUSION weights for cTWAS ...')
weights <- preprocess_weights(weights_dir,
                              region_info,
                              gwas_snp_ids = z_snp$id,
                              snp_map = snp_map, # Need
                              LD_map = LD_map, # Need
                              type = "expression",
                              context = cell_type,                         
                              weight_format = "FUSION",
                              fusion_method = "lasso", # Run Lasso for now
                              fusion_genome_version = "b38",
                              top_n_snps = NULL,
                              drop_strand_ambig = TRUE,
                              filter_protein_coding_genes = FALSE,
                              scale_predictdb_weights = FALSE,
                              load_predictdb_LD = FALSE,
                              ncore = 4)

# Then run ctwas (example for full analysis)
message('Run cTWAS ...')
ctwas_res <- ctwas_sumstats(z_snp, 
                            weights, 
                            region_info, 
                            LD_map, 
                            snp_map, 
                            thin = 1,
                            maxSNP = 20000,
                            min_group_size = 100, 
                            group_prior_var_structure = "shared_all", 
                            min_nonSNP_PIP = 0.5,
                            min_abs_corr = 0.1, 
                            ncore = 6, 
                            ncore_LD = 4,
                            save_cor = TRUE,
                            cor_dir = "./cor_matrix",
                            force_compute_cor = FALSE)

message('Structure of cTWAS res ...')
str(ctwas_res)

# Output: Save ctwas_res$region_pip (gene PIPs per region)
message('Write cTWAS pip output ...')
write_tsv(ctwas_res$region_pip, output)

message('All Done.')
#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------