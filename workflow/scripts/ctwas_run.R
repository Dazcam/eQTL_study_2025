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
# non-rsIDs
non_rsids <- gwas_tbl %>%
  filter(!str_starts(SNP, "rs")) %>%
  distinct(SNP) %>%
  nrow()
if (non_rsids > 0) {
  message(str_glue("Found {non_rsids} non-rsID SNPs. Removing them ..."))
  gwas_tbl <- gwas_tbl %>% filter(str_starts(SNP, "rs"))
} else {
  message("No non-rsID SNPs found.")
}

# duplicate rsIDs
duplicate_rsids <- gwas_tbl %>%
  group_by(SNP) %>%
  filter(n() > 1) %>%
  distinct(SNP) %>%
  nrow()
if (duplicate_rsids > 0) {
  message(str_glue("Found {duplicate_rsids} duplicate rsIDs. Keeping first occurrence ..."))
  gwas_tbl <- gwas_tbl %>% distinct(SNP, .keep_all = TRUE)
} else {
  message("No duplicate rsIDs found.")
}

# Alt Chrs 
message("Checking GWAS chromosomes ...\n")
gwas_tbl |> 
  group_by(CHR) |> 
  count() |> 
  print(n = Inf)

alt_chroms <- gwas_tbl %>%
  filter(str_ends(CHR, "alt")) %>%
  distinct(CHR) %>%
  pull(CHR)

if (length(alt_chroms) > 0) {
  message("\nChromosomes ending with 'alt' found: ", paste(alt_chroms, collapse = ", "))
  message("\nCount of SNPs with 'alt' chromosomes: ", 
          nrow(gwas_tbl %>% filter(str_ends(CHR, "alt"))))
  message("Removing SNPs with 'alt' chromosomes ... ")
  gwas_tbl <- gwas_tbl %>% filter(!str_ends(CHR, "alt"))
} else {
  message("No chromosomes ending with 'alt' found.")
}

message("\nAppending bim ref allele info to GWAS for alignment check...\n")
mrg_tbl <- gwas_tbl %>%
  left_join(snp_tbl, by = c("SNP"), suffix = c("", ".ref"))

mrg_tbl <- mrg_tbl |> 
  mutate(
  # classify allele relationships
  allele_status = case_when(
    A1 == A1.ref & A2 == A2.ref ~ "aligned",
    A1 == A2.ref & A2 == A1.ref ~ "flipped",
    TRUE                       ~ "mismatched"
  ),

  Z = case_when(
    allele_status == "flipped" ~ -Z,  
    TRUE                       ~ Z
  ),
  
  A1 = case_when(
    allele_status == "flipped" ~ A2,
    TRUE                       ~ A1
  ),
  A2 = case_when(
    allele_status == "flipped" ~ A1,
    TRUE                       ~ A2
  )
)

# Report allele flips
message("\n=== GWAS allele alignment summary ===")
mrg_tbl %>%
  count(allele_status) %>%
  mutate(percent = round(100 * n / sum(n), 1)) %>%
  arrange(desc(n)) %>%
  print(n = Inf)
message("========================================\n")

z_snp <- mrg_tbl |>
  select(id = SNP, A1, A2, z = Z) |>
  filter(if_all(everything(), ~ !is.na(.) & . != ""))

# Check for any remaining missing or blank values
if (anyNA(z_snp) || any(z_snp == "", na.rm = TRUE)) {
  stop("Missing values found in GWAS data.")
} else {
  message("No missing values in GWAS data ...")
}

message("Final number of SNPs in GWAS: ", nrow(z_snp) - 1)

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
weights_file <- str_glue('../results/12CTWAS/output/processed_weights/processed_weights_', cell_type, '_', gwas_trait, '.rds')


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
                                ncore = 16)
  
  message('Writing cTWAS weights for plotting later ...')
  dir.create(processed_weights_dir)
  write_rds(weights, paste0(processed_weights_dir, 'processed_weights_', cell_type, '_', gwas_trait, '.rds'))

}

### Adding a diagnostic section to troubleshoot failures. -----------------

message("Checking for genes with no SNPs in weights...")
n_before <- length(weights)

# Filter: Keep only genes where at least one SNP overlaps with this GWAS
weights <- weights[sapply(names(weights), function(gene_key) {
  gene_entry <- weights[[gene_key]]
  
  # Get SNPs for this gene (prefer $snps if present; fallback to rownames(wgt))
  snps_in_gene <- if ("snps" %in% names(gene_entry)) {
    gene_entry$snps
  } else {
    rownames(gene_entry$wgt)
  }
  
  length(intersect(snps_in_gene, z_snp$id)) > 0
})]

n_after <- length(weights)
dropped <- n_before - n_after

if (dropped > 0) {
  message(sprintf("Dropped %d genes with zero SNPs overlapping this GWAS (prevents cTWAS crash)", dropped))
} else {
  message("All genes have at least one SNP—proceeding normally.")
}

message("Checking for problematic R_wgt matrices in weights...")

bad_genes <- character(0)

for (gn in names(weights)) {
  Rmat <- weights[[gn]]$R_wgt
  
  # Case 1: R_wgt is NULL → this is a 1-SNP gene → always safe, keep it
  if (is.null(Rmat)) next
  
  # Case 2: R_wgt exists but is empty or corrupted (should never happen, but has been seen)
  if (!is.matrix(Rmat) || nrow(Rmat) == 0) {
    bad_genes <- c(bad_genes, gn)
    next
  }
  
  # Case 3: matrix ≥2×2 → check if singular
  if (nrow(Rmat) > 1) {
    det_val <- determinant(Rmat, logarithm = TRUE)$modulus
    # det_val can be -Inf, NA, or a very large negative number if singular
    if (is.na(det_val) || is.infinite(det_val) || det_val < -10) {
      bad_genes <- c(bad_genes, gn)
    }
  }
  # 1×1 matrices are always fine → do nothing
}

if (length(bad_genes) > 0) {
  message(sprintf("Dropping %d genes with singular/near-singular R_wgt:", length(bad_genes)))
  print(head(bad_genes, 10))
  weights <- weights[!names(weights) %in% bad_genes]
} else {
  message("No problematic R_wgt matrices found.")
}

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
                            ncore = 16, 
                            ncore_LD = 16,
                            save_cor = TRUE,
                            cor_dir = cor_dir,
                            force_compute_cor = FALSE)

# Output: Save ctwas_res$region_pip (gene PIPs per region)
message('Write cTWAS pip output ...')
write_rds(ctwas_res, output)

message('All Done.')
#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------