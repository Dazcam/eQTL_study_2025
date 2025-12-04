#--------------------------------------------------------------------------------------
#
#    Get TensorQTL LD proxies
#
#--------------------------------------------------------------------------------------

#  For each cell type 
#   1. Load cis-eQTL with qval < 0.05
#   2. For each SNP in 1, get it’s LD proxies using LDlinkR
#   3. Collate all rsIDs from step 1 and 2

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

message("\n\nGenerating LD proxies for tensorQTL perm eQTL with LDlinkR ...")

suppressPackageStartupMessages({
  library(tidyverse)
  library(LDlinkR)
})

# Input and output paths)
perm        <- snakemake@input[['perm']]
proxy_dir   <- snakemake@params[['proxy_dir']]
token_file  <- snakemake@params[["ldlink_token"]]
output      <- snakemake@output[[1]]
cell_type   <- snakemake@wildcards[["cell_type"]]
if (!file.exists(token_file)) stop("Token file not found: ", token_file)
token       <- read_lines(token_file, n_max = 1) %>% str_trim()

message(paste("\ncis-eQTL file: ", perm))
message(paste("Token file: ", token_file))
message(paste("Proxy dir: ", proxy_dir))
message(paste("Output will write to: ", output))
message(paste("Cell type: ", cell_type, '\n'))

dir.create(proxy_dir)
setwd(proxy_dir) # Need to set this as LDproxy_batch spits out files to wd

## Load data -----
message('Loading cis-eQTL for: ', cell_type)
cis_eQTL_tbl <- read_tsv(perm)
lead_variants <- cis_eQTL_tbl %>%
  filter(qval < 0.05) %>%
  pull(variant_id)
n_sig_cis <- length(lead_variants)
message(n_sig_cis, ' FDR sig. eQTL loaded.')

# Get proxies -----
message('Collecting LD proxies for: ', cell_type)
LDproxy_batch(snp = lead_variants,
              pop = "CEU",
              r2d = "r2",
              token = token,
              genome_build = "grch38")

# Load proxies from file and filter -----
message('Loading proxies from files ...')
perfect_proxies <- list.files(proxy_dir, pattern = "\\.txt$", full.names = TRUE) %>%
  set_names(~ basename(.) %>% str_remove("_grch38\\.txt$")) %>%
  map_dfr(~ {
    read_tsv(
      .x,
      col_names = c(
        "row_num",
        "RS_Number", "Coord", "Alleles", "MAF", "Distance",
        "Dprime", "R2", "Correlated_Alleles", "FORGEdb",
        "RegulomeDB", "Function"
      ),
      col_select = -row_num,
      col_types = cols(),
      skip = 1
    ) %>%
      mutate(query_snp = basename(.x) %>% str_remove("_grch38\\.txt$"))
  }, .id = NULL) %>%                       # no extra id column needed now
  filter(R2 == 1) %>%
  select(query_snp, RS_Number, everything()) %>%
  distinct()

message('Loaded ', nrow(perfect_proxies), ' proxies')

message('Removing non-biallelic proxies ... ')
perfect_proxies_rs_ids <- perfect_proxies |>
  filter(RS_Number != '.') |> # Non-biallelic SNP have rsID == '.'
  pull(RS_Number)

message(' ', length(perfect_proxies_rs_ids), ' remain.')

message('Collating cis-eQTL and proxy rsIDs ... ')
all_rs_ids <- c(lead_variants, perfect_proxies_rs_ids) |>
  unique()
message('  ', length(all_rs_ids),'rsIDs in total.')

## Write output -----
message('\nWriting file to: ', output)
write_tsv(perfect_proxies, output)

message('All Done.')

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------