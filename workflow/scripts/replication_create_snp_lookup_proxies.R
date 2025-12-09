#--------------------------------------------------------------------------------------
#
#    Create SNP lookup table for snATAC-seq enrichment
#
#--------------------------------------------------------------------------------------

# 1. Load sig. eQTL an restrict to FDR < 0.05
# 2. Load LD proxies and restrict to R^2 == 1
# 3. Intersect 1 and 2 to obtain final SNP list
# 4. Use BiomaRt to get hg38 SNP coordinates

# Notes: Needs to run on head node BiomaRt requires net access

## Info  ------------------------------------------------------------------------------

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

message("\n\nGenerate SNP lookup file for cell-specific sig. eQTL and LD proxies ...")

# -------------------------------------------------------------------------------------

library(tidyverse)
library(biomaRt)
library(GenomicRanges)

# Input and output paths
qtl_perm <- snakemake@input[["qtl_perm"]]
proxy_dir <- snakemake@params[[1]]
out_file <- snakemake@output[[1]]
message('qtl: ', qtl_perm)
message('prx: ', proxy_dir)
message('out: ', out_file)

# Check variable assignment
message("\nVariables")
cat("============================")
tibble(
  variable = c("qtl_perm", "proxy_dir", "out_file"),
  value    = c(qtl_perm, proxy_dir, out_file)) |> 
  knitr::kable(format = "simple", align = "l") |>
  print()
message("\n============================\n")

message("\n--- Generating SNP lookup file ---")
message("\nProcessing: ", qtl_perm)

# Load significant eQTL SNPs  -----
lead_variants <- read_tsv(qtl_perm) %>%
  filter(qval < 0.05) %>%
  pull(variant_id)

if (length(lead_variants) == 0) {
  message("No significant eQTLs found for ", cell_type)
  quit(save = "no", status = 0)
}

# Load proxies and filter -----
message('Loading proxies from files ...')
perfect_proxies <- list.files(proxy_dir, pattern = "_grch38\\.txt$", full.names = TRUE) %>%
  set_names(~ basename(.) %>% str_remove("_grch38\\.txt$")) %>%
  map_dfr(~ {
    suppressMessages(
      suppressWarnings(
        read_tsv(.x, show_col_types = FALSE)
      )
    ) %>%
      mutate(query_snp = basename(.x) %>% str_remove("_grch38\\.txt$"))
  }) %>%                       # no extra id column needed now
  filter(R2 == 1) %>%
  dplyr::select(query_snp, RS_Number, everything()) %>%
  distinct()

message('Loaded ', nrow(perfect_proxies), ' proxies')

message('Removing non-biallelic proxies ... ')
perfect_proxies_rs_ids <- perfect_proxies |>
  filter(RS_Number != '.') |> # Non-biallelic SNP have rsID == '.'
  pull(RS_Number)

message('  ', length(perfect_proxies_rs_ids), ' remain.\n')

message('Collating cis-eQTL and proxy rsIDs ... ')
all_rs_ids <- c(lead_variants, perfect_proxies_rs_ids) |>
  unique()
message('  ', length(all_rs_ids),' rsIDs in total.\n')


# Initialize Ensembl once
temp_cache <- tempdir()
Sys.setenv(BIOMART_CACHE = temp_cache)
biomartCacheClear() # As running on head node need to use temp cache
ensembl <- useEnsembl("ENSEMBL_MART_SNP", dataset = "hsapiens_snp")

# SNP lookup
message("Pulling out ", length(all_rs_ids), " SNPs with BiomaRt ...")
snp_lookup_hg38 <- getBM(
  attributes = c("refsnp_id", "chr_name", "chrom_start", "chrom_end"),
  filters = "snp_filter",
  values =  all_rs_ids,
  mart = ensembl,
  uniqueRows = TRUE
)
message('  ', nrow(snp_lookup_hg38), " SNPs annotated with coords.\n") 

message("Filtering SNPs wit chr 1:22, X, or Y ...")
snp_lookup_hg38_filt <- snp_lookup_hg38 %>%
  filter(chr_name %in% c(1:22, "X", "Y")) %>%
  mutate(chr_name = paste0("chr", chr_name)) %>%
  distinct()

message('. ', nrow(snp_lookup_hg38) - nrow(snp_lookup_hg38_filt), 
        ' eQTL removed due to non-standard chr annotation.\n')

message("Final number of annotated SNPs: ", nrow(snp_lookup_hg38_filt), '.')
message("Writing SNP list ...")
write_rds(snp_lookup_hg38_filt, out_file)
message('Done.')

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------