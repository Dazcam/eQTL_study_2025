#--------------------------------------------------------------------------------------
#
#    Extract sig. eGenes - Pass as a list to run_susie.R for each cell type
#
#--------------------------------------------------------------------------------------

## Load libraries and variables -------------------------------------------------------
library(tidyverse)

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

message('\n\nExtracting table of sig. eGenes to use as input for run_susie.R ...')
input_cs <- snakemake@input[['cred_set']]
input_bim <- snakemake@input[['bim']]
out_file <- snakemake@output[[1]]

# Extract cell type from input_cs
cell_type <- input_cs %>%
  basename() %>%
  str_extract("^[^.]+")  # Extracts "ExN-1"

# Extract chromosome and genome build from input_bim
bim_info <- input_bim %>%
  basename() %>%
  str_match(".*(hg\\d+).*\\.(\\d+)\\.bim$") %>%
  as_tibble() %>%
  transmute(build = V2, chr = V3)

message("Creating maxCPP ", bim_info$build, " annot file for ", 
        cell_type, " on chr ", bim_info$chr, "...")
  
susie_tbl <- read_tsv(input_cs)
bim_file <- read_tsv(input_bim, col_names = F)

message(nrow(susie_tbl), " variants in credible sets before MaxCPP calculation ...")

# For all variants in CS, maxCPP is the maximum CPP, or PIP here, of all eGenes
susie_tbl2 <- susie_tbl %>% 
  group_by(variant_id) %>% 
  mutate(maxCPP = max(pip)) %>% 
  select(chr, pos, variant_id, maxCPP)
susie_tbl2 <- susie_tbl2[!duplicated(susie_tbl2$variant_id),]
message(nrow(susie_tbl2), " variants in credible sets after MaxCPP calculation ...")

# Reference variants not in susie CS have maxCPP as 0
bim <- bim_file %>% left_join(susie1, by=c("X1" = "chr", "X4" = "pos"))
bim$maxCPP[ is.na( bim$maxCPP ) ] <- 0
bim <- bim %>% 
  select(CHR = X1, BP = X4, SNP = X2, CM = X3, maxCPP) 

message(nrow(bim |> filter(maxCPP > 0)), " variants with MaxCPP > 0 on chr ", 
        bim_info$chr, ' ...')

write.table(bim, gzfile(out_file), col.names = T, row.names = F, quote = F, sep = "\t")
message('All done.')

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
