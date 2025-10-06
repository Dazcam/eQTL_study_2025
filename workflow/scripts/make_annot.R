#--------------------------------------------------------------------------------------
#
#    Make annotation files (maxCPP and CS95) for SLDSR
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
out_maxcpp <- snakemake@output[['maxcpp_annot']]
out_cs95 <- snakemake@output[['cs95_annot']]

cat("============================")
tibble(
  variable = c("input_cs", "input_bim", "out_maxcpp", "out_cs95"),
  value    = c(input_cs, input_bim, out_maxcpp, out_cs95)) |> 
  knitr::kable(format = "simple", align = "l") |>
  print()
message("\n============================\n")

# Extract cell type from input_cs
cell_type <- input_cs %>%
  basename() %>%
  str_extract("^[^.]+")  # Extracts "ExN-DL"

# Extract chromosome and genome build from input_bim
bim_info <- input_bim %>%
  basename() %>%
  str_match(".*(hg\\d+).*\\.(\\d+)\\.bim$") %>%
  as_tibble() %>%
  transmute(build = V2, chr = V3)

message("Creating maxCPP and CS95 ", bim_info$build, " annot files for ", 
        cell_type, " on chr ", bim_info$chr, "...")

susie_tbl <- read_tsv(input_cs)
bim_file <- read_tsv(input_bim, col_names = F)

message(nrow(susie_tbl), " variants in credible sets before annotation calculation ...")

# Filter to target chr early
susie_tbl_chr <- susie_tbl %>% filter(chr == bim_info$chr)

# Compute maxCPP: max PIP across eGenes per variant
susie_maxcpp <- susie_tbl_chr %>% 
  group_by(variant_id) %>% 
  mutate(maxCPP = max(pip)) %>% 
  select(chr, pos, maxCPP) %>%  # Drop variant_id; join on chr/pos
  distinct()  # Dedup if needed

# Compute CS95: binary indicator (1 if in any CS, via distinct chr/pos)
susie_cs95 <- susie_tbl_chr %>% 
  distinct(chr, pos) %>% 
  mutate(CS95 = 1L)

message(nrow(susie_maxcpp), " unique variants for maxCPP on chr ", bim_info$chr, "...")
message(nrow(susie_cs95), " unique variants for CS95 on chr ", bim_info$chr, "...")

# Ensure the annotation files match the .bim file's SNPs and order
# Create a reference table from bim_file with necessary columns
bim_ref <- bim_file %>% 
  select(CHR = X1, BP = X4, SNP = X2, CM = X3) %>%
  mutate(row_index = row_number())  # Preserve original order

# For maxCPP: Left join on chr/pos and set NA to 0
bim_maxcpp <- bim_ref %>%
  left_join(susie_maxcpp, by = c("CHR" = "chr", "BP" = "pos")) %>%
  mutate(maxCPP = replace(maxCPP, is.na(maxCPP), 0)) %>%
  arrange(row_index) |>
  select(CHR, BP, SNP, CM, maxCPP)

# For CS95: Left join on chr/pos and set NA to 0
bim_cs95 <- bim_ref %>%
  left_join(susie_cs95, by = c("CHR" = "chr", "BP" = "pos")) %>%
  mutate(CS95 = replace(CS95, is.na(CS95), 0L)) %>%
  arrange(row_index) |>
  select(CHR, BP, SNP, CM, CS95)

# Check that the SNPs in the annotations match the .bim file
stopifnot(
  "maxCPP annotation does not contain the same SNPs as the .bim file" = 
    all.equal(bim_maxcpp$SNP, bim_ref$SNP, check.attributes = FALSE),
  "CS95 annotation does not contain the same SNPs as the .bim file" = 
    all.equal(bim_cs95$SNP, bim_ref$SNP, check.attributes = FALSE)
)

message(nrow(bim_maxcpp |> filter(maxCPP > 0)), " variants with maxCPP > 0 on chr ", 
        bim_info$chr, ' ...')
message(nrow(bim_cs95 |> filter(CS95 > 0)), " variants in any 95% CS on chr ", 
        bim_info$chr, ' ...')

# Write outputs
write.table(bim_maxcpp, gzfile(out_maxcpp), col.names = T, row.names = F, quote = F, sep = "\t")
write.table(bim_cs95, gzfile(out_cs95), col.names = T, row.names = F, quote = F, sep = "\t")
message('All done.')

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------