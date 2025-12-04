#--------------------------------------------------------------------------------------
#
#    Get TensorQTL LD proxies
#
#--------------------------------------------------------------------------------------

#  For each cell type 
#   1. Load cis-eQTL with qval < 0.05
#   2. For each SNP in 1, get it’s LD proxies using LDlinkR
#   3. Collate all rsIDs from step 1 and 2

# LDlinkR does not like concurrent requests so need to run this sequentially 
# for each cell type in a single script

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

## Load data -----
message('Loading cis-eQTL for: ', cell_type)
cis_eQTL_tbl <- read_tsv(perm)
lead_variants <- cis_eQTL_tbl %>%
  filter(qval < 0.05) %>%
  pull(variant_id)
n_sig_cis <- length(lead_variants)
message(n_sig_cis, ' FDR sig. eQTL loaded.')

# Set wd  -----
# LDproxy_batch spits out files to wd
# Need to do it after cis_eQTL_tbl is loaded
message('WD set to: ', proxy_dir)
dir.create(proxy_dir)
setwd(proxy_dir) 

# Get proxies -----
message('Collecting LD proxies for: ', cell_type)
# SOme SNPs fail when using batch, try another approach
# LDproxy_batch(snp = lead_variants,
#               pop = "CEU",
#               r2d = "r2",
#               token = token,
#               genome_build = "grch38")
failed_snps <- character(0) 
for (snp in lead_variants) {
  message("Submitting request for query variant ", snp, " ..")
  
  # Make sure we don't already have data for SNP
  # Careful here, wd is already proxy_dir
  if (file.exists(paste0(snp, "_grch38.txt"))) {
    message("Skipping ", snp, " (already done).")
    next
  }
  
  tryCatch(
    {
      LDproxy(snp = snp, pop = "CEU", r2d = "r2", token = token, genome_build = "grch38")
    },
    error = function(e) {
      message("Error querying SNP ", snp, ": ", e$message)
      failed_snps <<- c(failed_snps, snp)
      return(NULL)
    }
  )
}

# Save failed SNPs  
if (length(failed_snps) > 0) {
  failed_file <- file.path(proxy_dir, "failed_snps.txt")
  writeLines(failed_snps, failed_file)
  message("Written ", length(failed_snps), " failed SNPs to: ", failed_file)
} else {
  message("No SNPs failed.")
}




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