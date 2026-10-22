#--------------------------------------------------------------------------------------
#
#    Causal TWAS
#
#--------------------------------------------------------------------------------------

## Info  ------------------------------------------------------------------------------

#  Input
#  - GWAS sumstats per trait
#  - TWAS weights per cell type
#  - LD reference file (1000G hg38)

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
library(tidyverse)
library(ctwas)

# Load functions
source(file.path(snakemake@scriptdir, "ctwas_functions.R"))


# Input and output paths
gwas <- snakemake@input[['gwas']] # From 07prep_gwas.smk
weights_dir <- snakemake@params[['weights_dir']] # Folder containing .wgt.RDat files
ld_dir <- snakemake@params[['ld_dir']]
bim_file <- snakemake@input[['bim_file']]
cell_type <- snakemake@wildcards[['cell_type']]
gwas_trait <- snakemake@wildcards[['gwas']]
output <- as.character(snakemake@output[[1]])
out_dir <- dirname(output)
processed_weights_dir <- file.path(out_dir, 'processed_weights/')
cor_dir <- file.path(out_dir, paste0("cor_matrix_", gwas_trait, '_', cell_type))
threads <- snakemake@threads

## Skip known-failing cell_type / gwas combinations  ----------------------------------
# These combinations reliably fail (chiefly OOM during LD computation in
# preprocess_weights(), e.g. "cannot allocate vector of size 2020.4 Gb" on chr7).
# Rather than let snakemake mark the whole pipeline as failed, we write a
# clearly-flagged placeholder .rds and exit cleanly so downstream rules
# (ctwas_report) can detect and report these as skipped rather than crash.
skip_pairs <- tibble::tribble(
  ~cell_type,                 ~gwas,   ~reason,
  "Glu-UL-2",                 "ocd",   "OOM: cannot allocate vector (~2 Tb) computing LD on chr7",
  "GABA-2",                   "ocd",   "OOM: cannot allocate vector (~2 Tb) computing LD on chr7",
  "GABA-2",                   "bpd",   "OOM: cannot allocate vector (~2 Tb) computing LD on chr7",
  "Glu-UL-2",                 "adhd",  "OOM: cannot allocate vector (~2 Tb) computing LD on chr7",
  "GABA-2",                   "adhd",  "OOM: cannot allocate vector (~2 Tb) computing LD on chr7",
  "NPC",		      "adhd",  "Job terminated with no error message after chr17-19 LD step (likely OOM-killed by scheduler before R could report)",	
  "NPC-to-Glu-DL-Q4-Bin1",    "adhd",  "Job terminated with no error message after chr17-19 LD step (likely OOM-killed by scheduler before R could report)",
  "NPC-to-Glu-UL-Q4-Bin2",    "adhd",  "Job terminated with no error message after chr17-19 LD step (likely OOM-killed by scheduler before R could report)"
)

this_pair <- dplyr::filter(skip_pairs, cell_type == !!cell_type, gwas == !!gwas_trait)

if (nrow(this_pair) > 0) {

  message("=== SKIPPING known-failing combination ===")
  message("cell_type: ", cell_type)
  message("gwas: ", gwas_trait)
  message("Reason: ", this_pair$reason[1])

  placeholder <- list(
    status = "skipped",
    cell_type = cell_type,
    gwas_trait = gwas_trait,
    reason = this_pair$reason[1],
    skipped_at = as.character(Sys.time()),
    finemap_res = NULL,
    z_gene = NULL
  )

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  write_rds(placeholder, output)

  message("Placeholder output written to: ", output)
  message("=== Skip complete, exiting cleanly ===")
  quit(save = "no", status = 0)
}

# Read in data
message("\nLoading data ...\n")
message("GWAS loaded from: ", gwas)
message("weights_dir loaded from: ", weights_dir)
message("Bim file loaded from: ", bim_file)
message("cell_type: ", cell_type)
message("Output will be saved to: ", output)
message("Processed weights will be save to: ", processed_weights_dir)
message("Correlation matrices will be saved to: ", cor_dir)

# Prep GWAS
# Note: Need to harmonise alleles in GWAS, Ref and prediction weights
# cTWAS says: A1 is the alternate allele, and A2 is the reference allele
# TensorQTL, input for FUSION prediction weights, has ALT allele as A1
# In most GWAS, REF is A1 so will check and likely flip SNPs
# Note we need to use ../results/07PREP-GWAS/scz_hg38.tsv not LDSR munged sumsats as we need the alleles
message('\nLoading GWAS ...')
gwas_tbl <- read_tsv(gwas) 
gwas_n <- as.numeric(names(sort(table(gwas_tbl$N), decreasing = TRUE)[1]))
message('GWAS N is set to: ', gwas_n)
snp_tbl <- read_tsv(bim_file, col_names = c('Chr', 'SNP', 'CM', 'BP', 'A1', 'A2'))

## GWAS checks  -----------------------------------------------------------------------
gwas_tbl <- clean_gwas_for_ctwas(gwas_tbl)

z_snp <- align_gwas_to_ref(gwas_tbl, snp_tbl)

# Region info: included with package
# local - ctwas_0.1.38
#region_file <- system.file("extdata/ldetect", "EUR.b38.bed", package = "ctwas")
#region_info <- read_tsv(region_file)

# remote - ctwas_0.5.34
message('\nLoading region info from extdata ...')
region_file <- system.file("extdata/ldetect", "EUR.b38.ldetect.regions.RDS", package = "ctwas")
region_info <- readRDS(region_file)

# Load ref SNP info (for mapping) - think this is pnly used when not using LD matrices
#ref_snp_info <- data.table::fread(ref_snp_info_file, sep = "\t")
#class(ref_snp_info) <- "data.frame"

# Create metadata table for LD files (one row per region)
region_metatable <- make_ld_metatable(region_info, ld_dir, prefix = "1000G_EUR_Phase3_hg38")

# Generate snp_map and LD_map
message('Generating snp_map and LD_map for cTWAS ...')
res <- create_snp_LD_map(region_metatable)
snp_map <- res$snp_map
LD_map <- res$LD_map  # Not directly used in preprocess_weights, but for ctwas() later

# Optional: Harmonize GWAS first (add after reading z_snp)
# message('Preprocessing GWAS for cTWAS ...')
# z_snp <- preprocess_z_snp(z_snp, snp_map, 
#                           drop_multiallelic = TRUE, 
#                           drop_strand_ambig = TRUE)  # Converts to chr:pos_ref_alt

# FUSION weights
message('Preprocessing FUSION weights for cTWAS ...')
  
message("\nFirst check if cTWAS processed weight file already exists ... ")
weights_file <- str_glue('../results/12CTWAS/output/processed_weights/processed_weights_', 
                         cell_type, '_', gwas_trait, '.rds')

if (file.exists(weights_file)) {
  message("Yes! Loading processed weight file ...")
  weights <- read_rds(weights_file)
  
} else {
  
  message("No! Generating cTWAS ready weights ...")
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
                                ncore = threads)
  
}


if (cell_type %in%  c("InN-0", "InN-1") &
    gwas_trait %in% c("adhd", "ocd", "ptsd")) {
  
  message("=== Running full gene-by-gene QC before cTWAS ===")
  qc_tbl <- qc_ctwas_weights(weights, z_snp)

  message("\n=== QC Summary ===")
  print(qc_tbl %>% count(status, reason), n = Inf)

  failed_genes <- qc_tbl %>% filter(status == "FAIL") %>% pull(gene)
  passed_genes <- qc_tbl %>% filter(status == "PASS") %>% pull(gene)

  message("\nDropping ", length(failed_genes), " genes failing QC.")
  message("Keeping ", length(passed_genes), " genes.")

  weights_clean <- weights[passed_genes]

  # Save QC table for audit
  qc_file <- file.path(out_dir, paste0("ctwas_weight_qc_", cell_type, "_", gwas_trait, ".tsv"))
  readr::write_tsv(qc_tbl, qc_file)
  message("QC report saved to: ", qc_file)

  # Replace weights with cleaned version
  weights <- weights_clean

  message("=== QC complete. Proceeding to run cTWAS with cleaned weight set ===")
  
}

message('Writing cTWAS weights for plotting later ...')
dir.create(processed_weights_dir)
write_rds(weights, paste0(processed_weights_dir, 'processed_weights_', cell_type, '_', gwas_trait, '.rds'))


### ----------------------------------------------------------
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
                            ncore = threads, 
                            ncore_LD = threads,
                            save_cor = TRUE,
                            cor_dir = cor_dir,
                            force_compute_cor = FALSE)

# Output: Save ctwas_res$region_pip (gene PIPs per region)
message('Write cTWAS pip output ...')
write_rds(ctwas_res, output)

message('All Done.')

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
