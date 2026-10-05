#--------------------------------------------------------------------------------------
#
#    Generate per-GWAS significant SMR+HEIDI tables for manuscript
#
#--------------------------------------------------------------------------------------

# For each GWAS trait: SMR+HEIDI-significant (SNP, gene, cell type) results,
# unioned across all cell types, joined to gene symbols, written as
# {gwas}_smr.tsv.
#
# Also writes a single combined xlsx (Supplementary Table 7) with one sheet
# per GWAS trait (sheet name = full disorder name), Glu-UL/Glu-DL relabelled
# to Glu-A/Glu-B and "-Q4-" stripped from pseudotime bin cell type names in
# the "Cell Type" column -- same pattern as manuscript_table_ctwas_sig.R.
#
# This was previously generated silently as a side effect of rendering
# smr_report.Rmd (one write_tsv() call per disorder tab). Moved here so it's
# a tracked, reproducible Snakemake rule that runs once in the tables
# pipeline rather than every time the report happens to render.
# smr_report.Rmd now reads these files instead of recomputing them.
#
# Feeds: smr_tbl (all_smr_tbl.xlsx), smr_ctwas_venn_tbl (SMR vs cTWAS gene
# tables), and smr_report.Rmd's "SMR and HEIDI sig" / "Unique SNP-gene
# pairs" tabs.

## Info  ------------------------------------------------------------------------------

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

message("\n\nGenerating per-GWAS significant SMR tables for the manuscript ...")

# -------------------------------------------------------------------------------------
suppressMessages({
  library(tidyverse)
  library(openxlsx)
})

# --- Set variables
traits           <- snakemake@params[["traits"]]
cell_types       <- snakemake@params[["cell_types"]]
in_dir           <- snakemake@params[["in_dir"]]
gene_lookup_file <- snakemake@params[["gene_lookup_file"]]
p_smr            <- as.numeric(snakemake@params[["p_smr"]])
p_heidi          <- as.numeric(snakemake@params[["p_heidi"]])

# Snakemake output is now a named output: `tsv` (list, positionally ordered
# to match `traits` -- both come from expand(..., gwas = GWAS_TRAITS) against
# the same list) and `xlsx` (single combined Supplementary Table 7 path).
genes_out_files <- snakemake@output[["tsv"]]
xlsx_out        <- snakemake@output[["xlsx"]]

# traits           <- c("scz", "bpd", "mdd", "adhd", "ocd")
# cell_types       <- c("Glu-UL", "Glu-DL", "NPC", ...)
# in_dir           <- "../results/10SMR/smr/"
# gene_lookup_file <- "../resources/sheets/gene_lookup_hg38.tsv"
# p_smr            <- 0.05
# p_heidi          <- 0.01

message("=== Per-GWAS significant SMR tables: ", paste(traits, collapse = ", "), " ===")

if (length(genes_out_files) != length(traits)) {
  stop("Number of Snakemake outputs (", length(genes_out_files),
       ") doesn't match number of traits (", length(traits),
       ") -- check the expand(..., gwas = ...) list matches `traits` exactly.")
}

gene_lookup <- suppressMessages(read_tsv(gene_lookup_file, show_col_types = FALSE)) |>
  select(ensembl_gene_id, external_gene_name) |>
  distinct(ensembl_gene_id, .keep_all = TRUE)  # gene_lookup_hg38.tsv has >1 row for some
                                                # gene IDs; without this the left_join below
                                                # fans out and duplicates significant rows

# --- SMR + HEIDI significance for one (cell_type, gwas) pair. Same filter as
# smr_report.Rmd's former get_heidi(): SMR p < (p_smr threshold / n_probes)
# [Bonferroni across probes tested in that cell type] AND HEIDI p > p_heidi.
get_heidi <- function(cell_type, in_dir, gwas, p_smr, p_heidi) {
  file_path <- file.path(in_dir, cell_type, paste0(cell_type, "_", gwas, ".smr"))
  if (!file.exists(file_path)) {
    message("  [", cell_type, "] file not found: ", file_path)
    return(NULL)
  }
  smr <- suppressMessages(read_tsv(file_path, show_col_types = FALSE))
  n_probes <- length(unique(smr$probeID))
  smr |>
    filter(p_SMR < (p_smr / n_probes), p_HEIDI > p_heidi & !is.na(p_HEIDI)) |>
    mutate(Cell_Type = cell_type)
}

# Columns the report expects, even when a trait has zero significant results
# (e.g. ADHD) -- keeps smr_report.Rmd's read_tsv() uniform across all traits
# with no per-trait special-casing needed.
empty_final_tbl <- tibble(
  Chr = character(), Symbol = character(), `Ensembl ID` = character(),
  `Cell Type` = character(), SNP = character(), A1 = character(), A2 = character(),
  Beta = double(), SE = double(), `P (SMR)` = double(), `P (HEIDI)` = double(),
  `NSNP HEIDI` = double()
)

# --- Relabel helper (display only, xlsx output): Glu-UL -> Glu-A,
# Glu-DL -> Glu-B, and strip "-Q4-" from pseudotime bin names, e.g.
# "NPC-to-Glu-DL-Q4-Bin1" -> "NPC-to-Glu-B-Bin1". Matches the convention used
# across the rest of the manuscript pipeline (bulk_egene_overlaps_plt.R,
# manuscript_plots_fig4.py, manuscript_table_ctwas_sig.R, etc.) The per-trait
# TSVs above are left untouched -- this only affects the new xlsx.
relabel_cell_type <- function(x) {
  x <- str_replace(x, "Glu-UL", "Glu-A")
  x <- str_replace(x, "Glu-DL", "Glu-B")
  x <- str_replace(x, "-Q4-", "-")
  x
}

# --- Full disorder names for xlsx sheet names (same mapping used in
# manuscript_table_ctwas_sig.R)
disorder_labels <- c(
  scz  = "Schizophrenia",
  bpd  = "Bipolar disorder",
  mdd  = "Major depressive disorder",
  adhd = "ADHD",
  ocd  = "Obsessive-compulsive disorder"
)

sig_sheets <- list()

for (i in seq_along(traits)) {
  gwas <- traits[i]
  message("--- ", gwas, " ---")

  all_sig_results <- map_dfr(cell_types, get_heidi, in_dir = in_dir, gwas = gwas,
                              p_smr = p_smr, p_heidi = p_heidi)

  if (is.null(all_sig_results) || ncol(all_sig_results) == 0) {
    # Genuinely no data: every cell type's .smr file was missing/unreadable
    # for this GWAS, so there's no schema to work with at all.
    message(gwas, ": no readable SMR files for any cell type")
    final_tbl <- empty_final_tbl
  } else if (nrow(all_sig_results) == 0) {
    # Files existed and were read fine, just nothing passed the
    # significance filter (e.g. ADHD) - filter() still preserves the full
    # column schema (A1, A2, etc.) even with 0 rows, so run the normal
    # pipeline rather than substituting the hand-typed fallback, which
    # would silently drop columns the real schema actually has.
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

  # --- xlsx sheet: same table, Cell Type relabelled for display
  # (Glu-UL/Glu-DL -> Glu-A/Glu-B, "-Q4-" stripped from pseudotime bin
  # names). Blank gene symbols become NA, consistent with the cTWAS xlsx.
  sheet_tbl <- final_tbl |>
    mutate(
      `Cell Type` = relabel_cell_type(`Cell Type`),
      Symbol      = na_if(Symbol, "")
    )

  sig_sheets[[disorder_labels[[gwas]]]] <- sheet_tbl
}

message("Writing combined xlsx (Supplementary Table 7): ", xlsx_out)
write.xlsx(
  sig_sheets,
  file = xlsx_out,
  overwrite = TRUE,
  headerStyle = createStyle(textDecoration = "bold")
)

message("Done: ", paste(traits, collapse = ", "))

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
