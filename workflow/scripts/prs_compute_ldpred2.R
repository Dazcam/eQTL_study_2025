#--------------------------------------------------------------------------------------
#
#    Compute Polygenic Risk Scores using LDpred2
#
#--------------------------------------------------------------------------------------

## Setup logging for Snakemake
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

message("\n\nComputing PRS with LDpred2...")

# -------------------------------------------------------------------------------------

library(tidyverse)
library(bigsnpr)

geno_prefix <- snakemake@params[["geno_prefix"]]
gwas_file <- snakemake@input[["gwas"]]
hm3_file <- snakemake@input[["hm3"]]
out_file <- snakemake@output[[1]]

# ----- 1. Attach genotype -----
message("Loading PLINK genotypes ...")
snp_readBed(paste0(geno_prefix, ".bed"))
obj.bigSNP <- snp_attach(paste0(geno_prefix, ".rds"))

G <- obj.bigSNP$genotypes
CHR <- obj.bigSNP$map$chromosome
POS <- obj.bigSNP$map$physical.pos
SNP_ID <- obj.bigSNP$map$rsid

# ----- 2. Load Data -----
message("Loading HapMap3 SNP list ...")
hm3_snps <- read_tsv(hm3_file, col_names = "rsid", show_col_types = FALSE)

message("Loading GWAS summary stats ...")
gwas_tbl <- read_tsv(gwas_file) |>
  as_tibble() |>
  rename(chr = CHR, pos = BP, beta = BETA, beta_se = SE,
         a1 = A1, a0 = A2, rsid = SNP, n_eff = N) |>
  inner_join(hm3_snps, by = "rsid")

message("GWAS variants: ", nrow(gwas_tbl))

# ----- 3. SNP matching -----
message("Matching GWAS and genotype SNPs ...")
info_snp <- snp_match(
  gwas_tbl,
  obj.bigSNP$map,
  strand_flip = TRUE
)

message("Matched SNPs: ", nrow(info_snp))

# ----- 4. Compute LD matrix -----
message("Computing LD correlations ...")

corr <- snp_cor(
  G,
  infos.pos = POS,
  size = 3 / 1000
)

# ----- 5. LDpred2 auto -----
message("Running LDpred2-auto ...")

beta_auto <- snp_ldpred2_auto(
  corr,
  info_snp$beta,
  info_snp$beta_se,
  n = info_snp$n_eff
)

beta_est <- beta_auto[[1]]$beta_est

# ----- 6. Compute PRS -----
message("Computing PRS scores ...")

prs <- big_prodVec(G, beta_est)
prs_tbl <- tibble(
  sample = obj.bigSNP$fam$sample.ID,
  PRS = prs
)

# ----- 7. Write output -----
message("Writing PRS file ...")
write_tsv(prs_tbl, out_file)

message("PRS computation complete.")

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------