#--------------------------------------------------------------------------------------
#
#    Causal TWAS - Run cTWAS
#
#--------------------------------------------------------------------------------------
#
# Pipeline:  12CTWAS | Rule: run_ctwas (one job per cell type x GWAS)
#            Upstream:   07PREP-GWAS sumstats; copy_fusion_weights (11TWAS FUSION
#                        weights); create_ld_matrices
#            Downstream: ctwas_report; 13MANUSCRIPT_TABLES ctwas_sig_tbl,
#                        smr_ctwas_venn_tbl
#
# Purpose:   For one cell type x GWAS:
#            - Skip known-failing combinations, writing a flagged placeholder .rds
#            - Clean the GWAS (non-rsIDs, duplicates, alt contigs) and align its
#              alleles to the reference .bim (see ctwas_functions.R)
#            - Build the SNP / LD maps from the LD matrices
#            - Preprocess the FUSION lasso weights (reused if already processed)
#            - Optionally QC the weights gene by gene (run_weight_qc in config)
#            - Run ctwas_sumstats() and save the results
#
# Inputs:    gwas         {gwas}_hg38.tsv (needs alleles, so not LDSR-munged)
#            weights_dir  FUSION .wgt.RDat files for the cell type
#            ld_dir       LD matrices (.RDS) and SNP info (.Rvar) per LD region
#            bim_file     1000G EUR hg38 .bim used for allele alignment
#            run_weight_qc  TRUE/FALSE: gene-by-gene weight QC before cTWAS
#
# Outputs:   output       ctwas_{cell_type}_{gwas}_ctwas.rds (cTWAS results)
#            Also written to dirname(output):
#            processed_weights/processed_weights_{cell_type}_{gwas}.rds
#            cor_matrix_{gwas}_{cell_type}/  (correlation matrices)
#            ctwas_weight_qc_{cell_type}_{gwas}.tsv  (only when QC runs)
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
message('\n\nRunning causal TWAS ...')

##  Load packages, functions and variables  -------------------------------------------
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
run_weight_qc <- isTRUE(snakemake@params[['run_weight_qc']])

# Print what each variable is set to
message("\nVariables")
message("============================")
print(data.frame(
  variable = c("gwas", "weights_dir", "ld_dir", "bim_file", "cell_type", "gwas_trait",
               "output", "processed_weights_dir", "cor_dir", "threads", "run_weight_qc"),
  value    = c(gwas, weights_dir, ld_dir, bim_file, cell_type, gwas_trait,
               output, processed_weights_dir, cor_dir, threads, run_weight_qc)),
           right = FALSE, row.names = FALSE)
message("============================\n")

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

## Load and prepare GWAS  ---------------------------------------------------------------
# Prep GWAS
# Note: Need to harmonise alleles in GWAS, Ref and prediction weights
# cTWAS says: A1 is the alternate allele, and A2 is the reference allele
# TensorQTL, input for FUSION prediction weights, has ALT allele as A1
# In most GWAS, REF is A1 so will check and likely flip SNPs
# Note we need to use ../results/07PREP-GWAS/scz_hg38.tsv not LDSR munged sumsats as we need the alleles
message('\nLoading GWAS ...')
gwas_tbl <- read_tsv(gwas, show_col_types = FALSE)
snp_tbl <- read_tsv(bim_file, col_names = c('Chr', 'SNP', 'CM', 'BP', 'A1', 'A2'),
                   show_col_types = FALSE)

## GWAS checks  -----------------------------------------------------------------------
gwas_tbl <- clean_gwas_for_ctwas(gwas_tbl)

z_snp <- align_gwas_to_ref(gwas_tbl, snp_tbl)

## LD regions and SNP / LD maps  -------------------------------------------------------
# Region info: included with package as of v0.5.34
message('\nLoading region info from extdata ...')
region_file <- system.file("extdata/ldetect", "EUR.b38.ldetect.regions.RDS", package = "ctwas")
region_info <- readRDS(region_file)

# Create metadata table for LD files (one row per region)
region_metatable <- make_ld_metatable(region_info, ld_dir, prefix = "1000G_EUR_Phase3_hg38")

# Generate snp_map and LD_map
message('Generating snp_map and LD_map for cTWAS ...')
res <- create_snp_LD_map(region_metatable)
snp_map <- res$snp_map
LD_map <- res$LD_map  # Not directly used in preprocess_weights, but for ctwas() later

## FUSION weights  ---------------------------------------------------------------------
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


## Optional weight QC  ------------------------------------------------------------------
# Used to identify the failing cell type / GWAS combinations now in skip_pairs.
# Off by default; switch on with run_ctwas: run_weight_qc: true in the config.
if (run_weight_qc) {

  message("=== Running full gene-by-gene QC before cTWAS ===")
  qc_tbl <- qc_ctwas_weights(weights, z_snp)

  message("\n=== QC Summary ===")
  print(qc_tbl |> count(status, reason), n = Inf)

  failed_genes <- qc_tbl |> filter(status == "FAIL") |> pull(gene)
  passed_genes <- qc_tbl |> filter(status == "PASS") |> pull(gene)

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
dir.create(processed_weights_dir, recursive = TRUE, showWarnings = FALSE)
write_rds(weights, paste0(processed_weights_dir, 'processed_weights_', cell_type, '_', gwas_trait, '.rds'))


## Run cTWAS  --------------------------------------------------------------------------
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
