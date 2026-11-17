#--------------------------------------------------------------------------------------
#
#    Significant SMR + HEIDI results per GWAS (Supp Table 9)
#
#--------------------------------------------------------------------------------------
#
# Pipeline:  16MANUSCRIPT-TABLES | Rule: smr_sig_tbl
#            Upstream:   12 SMR (.smr results per cell type and GWAS)
#            Downstream: smr_ctwas_venn_tbl (16); 17 smr_ctwas_venn_fig (Fig 6);
#                        12 smr_report.Rmd; Supp Table 9
#
# Purpose:   For each GWAS trait:
#            - Keep SMR + HEIDI significant results in every cell type
#              (Bonferroni over probes tested in that cell type, HEIDI P > p_heidi)
#            - Add gene symbols and write {gwas}_smr.tsv (pipeline labels)
#            - Add one sheet per trait to a combined xlsx (manuscript labels)
#
# Inputs:    in_dir            12 SMR results directory ({cell_type}/{cell_type}_{gwas}.smr)
#            gene_lookup_file  Ensembl ID to gene symbol lookup
#            p_smr, p_heidi    Significance thresholds (config)
#
# Outputs:   tsv               {gwas}_smr.tsv, one per trait, same order as traits
#            xlsx              Combined table, one sheet per disorder
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

message("\n\nBuilding significant SMR tables ...")

##  Load packages, functions and variables  -------------------------------------------
suppressMessages({
  library(tidyverse)
  library(openxlsx)
})

# Load functions (relabel_cell_types)
source(file.path(snakemake@scriptdir, "functions.R"))

# Inputs, outputs and settings
traits            <- snakemake@params[["traits"]]
cell_types        <- snakemake@params[["cell_types"]]
in_dir            <- snakemake@params[["in_dir"]]
gene_lookup_file  <- snakemake@params[["gene_lookup_file"]]
manuscript_labels <- snakemake@params[["labels"]]
p_smr             <- as.numeric(snakemake@params[["p_smr"]])
p_heidi           <- as.numeric(snakemake@params[["p_heidi"]])
genes_out_files   <- snakemake@output[["tsv"]]
xlsx_out          <- snakemake@output[["xlsx"]]

# Print what each variable is set to
message("\nVariables")
message("============================")
tibble(
  variable = c("traits", "cell_types", "in_dir", "gene_lookup_file", "labels",
               "p_smr", "p_heidi", "genes_out_files", "xlsx_out"),
  value    = c(paste(traits, collapse = ", "), paste(cell_types, collapse = ", "),
               in_dir, gene_lookup_file,
               paste(names(manuscript_labels), unlist(manuscript_labels), sep = " -> ", collapse = ", "),
               p_smr, p_heidi, paste(genes_out_files, collapse = ", "), xlsx_out)) |>
  knitr::kable(format = "simple") |>
  cat(sep = "\n")
message("============================\n")

if (length(genes_out_files) != length(traits)) {
  stop("Number of Snakemake outputs (", length(genes_out_files),
       ") doesn't match number of traits (", length(traits),
       ") -- check the expand(..., gwas = ...) list matches `traits` exactly.")
}

# gene_lookup_hg38.tsv has >1 row for some gene IDs; without distinct() the
# left_join below duplicates significant rows
gene_lookup <- read_tsv(gene_lookup_file, show_col_types = FALSE) |>
  dplyr::select(ensembl_gene_id, external_gene_name) |>
  distinct(ensembl_gene_id, .keep_all = TRUE)

##  Functions  ------------------------------------------------------------------------
# SMR + HEIDI significance for one (cell_type, gwas) pair: SMR p < p_smr / n_probes
# (Bonferroni across probes tested in that cell type) and HEIDI p > p_heidi
get_heidi <- function(cell_type, in_dir, gwas, p_smr, p_heidi) {
  file_path <- file.path(in_dir, cell_type, paste0(cell_type, "_", gwas, ".smr"))
  if (!file.exists(file_path)) {
    message("  [", cell_type, "] file not found: ", file_path)
    return(NULL)
  }
  smr <- read_tsv(file_path, show_col_types = FALSE)
  n_probes <- length(unique(smr$probeID))
  smr |>
    filter(p_SMR < (p_smr / n_probes), p_HEIDI > p_heidi & !is.na(p_HEIDI)) |>
    mutate(Cell_Type = cell_type)
}

# Columns written when no .smr file could be read for a trait, so every
# {gwas}_smr.tsv has the same header
empty_final_tbl <- tibble(
  Chr = character(), Symbol = character(), `Ensembl ID` = character(),
  `Cell Type` = character(), SNP = character(), A1 = character(), A2 = character(),
  Beta = double(), SE = double(), `P (SMR)` = double(), `P (HEIDI)` = double(),
  `NSNP HEIDI` = double()
)

# Full disorder names for the xlsx sheet names (as in manuscript_table_ctwas_sig.R)
disorder_labels <- c(
  scz  = "Schizophrenia",
  bpd  = "Bipolar disorder",
  mdd  = "Major depressive disorder",
  adhd = "ADHD",
  ocd  = "Obsessive-compulsive disorder"
)

##  Significant results per trait  ----------------------------------------------------
sig_sheets <- list()

for (i in seq_along(traits)) {
  gwas <- traits[i]
  message("--- ", gwas, " ---")

  all_sig_results <- map_dfr(cell_types, get_heidi, in_dir = in_dir, gwas = gwas,
                              p_smr = p_smr, p_heidi = p_heidi)

  if (is.null(all_sig_results) || ncol(all_sig_results) == 0) {
    # Every cell type's .smr file was missing for this GWAS: no schema to use
    message(gwas, ": no readable SMR files for any cell type")
    final_tbl <- empty_final_tbl
  } else if (nrow(all_sig_results) == 0) {
    # Files were read but nothing passed (e.g. ADHD): filter() keeps the full
    # column schema with 0 rows, so the normal formatting is used
    message(gwas, ": no significant SMR+HEIDI results for any cell type")
    final_tbl <- all_sig_results |>
      left_join(gene_lookup, by = join_by(probeID == ensembl_gene_id)) |>
      relocate(Cell_Type, ensembl_id = probeID, symbol = external_gene_name) |>
      distinct() |>
      rename(Chr = topSNP_chr, Symbol = symbol, `Ensembl ID` = ensembl_id,
             `Cell Type` = Cell_Type, SNP = topSNP, Beta = b_SMR, SE = se_SMR,
             `P (SMR)` = p_SMR, `P (HEIDI)` = p_HEIDI, `NSNP HEIDI` = nsnp_HEIDI)
  } else {
    final_tbl <- all_sig_results |>
      left_join(gene_lookup, by = join_by(probeID == ensembl_gene_id)) |>
      relocate(Cell_Type, ensembl_id = probeID, symbol = external_gene_name) |>
      distinct() |>
      rename(Chr = topSNP_chr, Symbol = symbol, `Ensembl ID` = ensembl_id,
             `Cell Type` = Cell_Type, SNP = topSNP, Beta = b_SMR, SE = se_SMR,
             `P (SMR)` = p_SMR, `P (HEIDI)` = p_HEIDI, `NSNP HEIDI` = nsnp_HEIDI)

    n_pairs <- final_tbl |> distinct(`Ensembl ID`, SNP) |> nrow()
    message(gwas, ": ", nrow(final_tbl), " significant (SNP, gene, cell type) rows across ",
            n_distinct(final_tbl$`Cell Type`), " cell types, ", n_pairs,
            " unique SNP-gene pairs")
  }

  write_tsv(final_tbl, genes_out_files[[i]])
  message("Wrote: ", genes_out_files[[i]])

  # xlsx sheet: same table with manuscript cell type labels; blank gene
  # symbols become NA, as in the cTWAS xlsx
  sheet_tbl <- final_tbl |>
    mutate(
      `Cell Type` = relabel_cell_types(`Cell Type`, manuscript_labels),
      Symbol      = na_if(Symbol, "")
    )

  sig_sheets[[disorder_labels[[gwas]]]] <- sheet_tbl
}

##  Write combined xlsx  --------------------------------------------------------------
message("Writing combined xlsx: ", xlsx_out)
write.xlsx(
  sig_sheets,
  file = xlsx_out,
  overwrite = TRUE,
  headerStyle = createStyle(textDecoration = "bold")
)

message("Done.")

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
