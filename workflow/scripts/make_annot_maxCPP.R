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
  select(chr, pos, variant_id, maxCPP) %>%
  filter(chr == bim_info$chr)
susie_tbl2 <- susie_tbl2[!duplicated(susie_tbl2$variant_id),]
message(nrow(susie_tbl2), " variants in credible sets after MaxCPP calculation for chr ", bim_info$chr,"...")

# Ensure the annotation file matches the .bim file's SNPs and order
# Create a reference table from bim_file with necessary columns
bim_ref <- bim_file %>% 
  select(CHR = X1, BP = X4, SNP = X2, CM = X3) %>%
  mutate(row_index = row_number())  # Preserve original order

# Perform a left join to keep all bim SNPs, assigning maxCPP = 0 for non-CS SNPs
bim_final <- bim_ref %>%
  left_join(susie_tbl2, by = c("CHR" = "chr", "BP" = "pos")) %>%
  mutate(maxCPP = replace(maxCPP, is.na(maxCPP), 0)) %>%
  arrange(row_index)  

# Check that the SNPs in the annotation match the .bim file
stopifnot(
  "Annotation file does not contain the same SNPs as the .bim file" = 
    all.equal(bim$SNP, bim_ref$SNP, check.attributes = FALSE)
)

message(nrow(bim_final |> filter(maxCPP > 0)), " variants with MaxCPP > 0 on chr ", 
        bim_info$chr, ' ...')

write.table(bim_final, gzfile(out_file), col.names = T, row.names = F, quote = F, sep = "\t")
message('All done.')

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
