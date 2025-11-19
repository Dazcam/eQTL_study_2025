#--------------------------------------------------------------------------------------
#
#    Post-TensorQTL transQTL analysis
#
#--------------------------------------------------------------------------------------

#  For each cell type 
#   1. Load cis-eQTL with qval < 0.05
#   2. For each SNP in 1, get it’s LD proxies (R^2 == 1) using LDlinkR
#   3. Collate all rsIDs from step 1 and 2 and filter trans eQTL on rsID overlaps 
#   4. Filter BF sig. trans eQTL

# Note: for now just pulling out rsIDs that were lead SNPs in cis run from trans data

#--------------------------------------------------------------------------------------
# Set up logging for Snakemake
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

message("\n\nRunning post-TensorQTL trans-QTL analysis ...")

suppressPackageStartupMessages({
  library(tidyverse)
  library(LDlinkR)
})

# Input and output paths)
perm        <- snakemake@input[['perm']]
trans       <- snakemake@input[['trans']]
proxy_dir   <- snakemake@params[['proxy_dir']]
token_file  <- snakemake@params[["ldlink_token"]]
output      <- snakemake@output[[1]]
cell_type   <- snakemake@wildcards[["cell_type"]]
token       <- read_lines(token_file, n_max = 1) %>% str_trim()

#dir.create(proxy_dir)
#setwd(proxy_dir) # Need to set this as LDproxy_batch spits out files to wd

## Load data -----
message('Loading cis-eQTL for: ', cell_type)
cis_eQTL_tbl <- read_tsv(perm)
lead_variants <- cis_eQTL_tbl %>%
  filter(qval < 0.05) %>%
  pull(variant_id)
n_sig_cis <- length(lead_variants)
message(' ', n_sig_cis, ' FDR sig. eQTL loaded.')


message('Loading trans-eQTL for: ', cell_type)
trans_raw <- read_tsv(trans) 
message(' ', nrow(trans_raw), ' trans eQTL loaded.')

## Get proxies -----
# message('Collecting LD proxies for: ', cell_type)
# LDproxy_batch(snp = lead_variants, 
#               pop = "CEU", 
#               r2d = "r2", 
#               token = token, 
#               genome_build = "grch38")

## Load proxies from file and filter -----
# message('Loading proxies from files ...')
# perfect_proxies <- list.files(proxy_dir, pattern = "\\.txt$", full.names = TRUE) %>%
#   set_names(~ basename(.) %>% str_remove("_grch38\\.txt$")) %>%
#   map_dfr(~ {
#     read_tsv(
#       .x,
#       col_names = c(                     
#         "row_num",
#         "RS_Number", "Coord", "Alleles", "MAF", "Distance",
#         "Dprime", "R2", "Correlated_Alleles", "FORGEdb",
#         "RegulomeDB", "Function"
#       ),
#       col_select = -row_num,
#       col_types = cols(),
#       skip = 1
#     ) %>%
#       mutate(query_snp = basename(.x) %>% str_remove("_grch38\\.txt$"))
#   }, .id = NULL) %>%                       # no extra id column needed now
#   filter(R2 == 1) %>%
#   select(query_snp, RS_Number, everything()) %>%
#   distinct()
# 
# message('Loaded ', nrow(perfect_proxies), ' proxies')
# 
# message('Removing non-biallelic proxies ... ')
# perfect_proxies_rs_ids <- perfect_proxies |>
#   filter(RS_Number != '.') |> # Non-biallelic SNP have rsID == '.'
#   pull(RS_Number)
# 
# message(' ', nrow(perfect_proxies_rs_ids), ' remain.')
# 
# message('Collating cis-eQTL and proxy rsIDs ... ')
# all_rs_ids <- c(lead_variants, perfect_proxies_rs_ids) |>
#   unique()
# message('  ', length(all_rs_ids),'rsIDs in total.')

## Intersect rsIDs with trans data -----
bf_p <- 0.05 / (n_sig_cis * length(lead_variants))
message('BF correction factor set to: ', bf_p)

message('Extracting BF sig. trans eQTL ... ')
trans_sig <- trans_raw %>%
  filter(variant_id %in% lead_variants) |>
  arrange(pval) |>
  filter(pval < bf_p) 

message('Number of BF sig. trans eQTL: ', nrow(trans_sig))

## Annotate table -----
# Add this, or maybe add it in markdown
# message('Add annotations to trans-eQTL tbl ... ')
# gene_lookup <- read_tsv('~/Desktop/eQTL_study_2025/resources/sheets/gene_lookup_hg38.tsv')
# cis_lookup <- cis_eQTL_tbl |>
#   filter(qval < 0.05) |>
#   left_join(gene_lookup, by = join_by(phenotype_id == ensembl_gene_id)) |>
#   select(variant_id, cis_ensembl_id = phenotype_id, cis_symbol = external_gene_name, 
#          cis_chr = chromosome_name)
# 
# trans_sig |>
#   left_join(gene_lookup, by = join_by(phenotype_id == ensembl_gene_id)) |>
#   left_join(cis_lookup, by = join_by(variant_id)) |>
#   select(variant_id, trans_ensembl_id = phenotype_id, trans_symbol = external_gene_name,
#         trans_chr = chromosome_name, cis_ensembl_id, cis_symbol, cis_chr) |>
#   print(n = Inf)
  
## Write output -----
message('Writing file to: ', output)
write_tsv(trans_sig, output)

message('All Done.')




#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------