#!/bin/bash
#SBATCH --job-name=snp_lookup
#SBATCH --account=SCWF00021_n_bray_19
#SBATCH --partition=htc_genoa
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=01:00:00
#SBATCH --output=snp_lookup_%j.out
#SBATCH --error=snp_lookup_%j.err

# query_eqtl_nom_files.sh
# Looks up a fixed list of SNPs across all *_cis_eQTL_nominal.tsv.gz files
# for the 19 non-trajectory cell types in
# 13MANUSCRIPT_PLOTS_TABLES/data_sharing/, and writes matching rows to a
# single combined TSV (snp_lookup_combined.tsv) with a cell_type column,
# retaining all original columns. Also writes a summary log
# (snp_lookup_summary.log) with per-SNP and per-cell-type match counts.
#
# Run from: workflow/  (i.e. same directory you'd normally run snakemake from)
# Usage:    sbatch query_eqtl_nom_files.sh

set -euo pipefail

SNPS="rs9274407 rs11749126 rs35904745 rs3129305"
BASE="../results/13MANUSCRIPT_PLOTS_TABLES/data_sharing"
OUT="snp_lookup_combined.tsv"
LOG="snp_lookup_summary.log"

CELL_TYPES="GABA GABA-0 GABA-1 GABA-2 Glu-DL Glu-DL-0 Glu-DL-1 Glu-DL-2 \
Glu-UL Glu-UL-0 Glu-UL-1 Glu-UL-2 MG NPC NPC-0 NPC-1 NPC-2 OPC Endo-Peri"

# Build an awk SNP-matching condition (SNP id is column 4 in these files)
SNP_COND=$(echo "$SNPS" | awk '{for(i=1;i<=NF;i++) printf "%s$4==\"%s\"", (i>1?" || ":""), $i}')

# Write header once (prefix with cell_type)
FIRST_FILE="$BASE/GABA_cis_eQTL_nominal.tsv.gz"
echo -e "cell_type\t$(zcat "$FIRST_FILE" | head -n1)" > "$OUT"

MISSING_FILES=""
for ct in $CELL_TYPES; do
  f="$BASE/${ct}_cis_eQTL_nominal.tsv.gz"
  if [[ -f "$f" ]]; then
    zcat "$f" | awk -F'\t' -v ct="$ct" -v cond="$SNP_COND" '
      NR>1 && ('"$SNP_COND"') { print ct"\t"$0 }
    ' >> "$OUT"
  else
    echo "WARNING: missing file for $ct: $f" >&2
    MISSING_FILES="$MISSING_FILES $ct"
  fi
done

echo "Done. Output: $OUT"

# --- Summary log ---
{
  echo "SNP lookup summary"
  echo "Run date: $(date)"
  echo "SNPs queried: $SNPS"
  echo "Cell types queried (${#CELL_TYPES}): $CELL_TYPES"
  echo ""

  echo "Cell types with missing input file:"
  if [[ -z "$MISSING_FILES" ]]; then
    echo "  none"
  else
    for ct in $MISSING_FILES; do echo "  $ct"; done
  fi
  echo ""

  echo "Total matching rows: $(($(wc -l < "$OUT") - 1))"
  echo ""

  echo "Matches per SNP (across all cell types):"
  for snp in $SNPS; do
    n=$(awk -F'\t' -v s="$snp" 'NR>1 && $5==s' "$OUT" | wc -l)
    echo "  $snp: $n"
  done
  echo ""

  echo "Matches per cell type:"
  for ct in $CELL_TYPES; do
    n=$(awk -F'\t' -v c="$ct" 'NR>1 && $1==c' "$OUT" | wc -l)
    echo "  $ct: $n"
  done
} > "$LOG"

echo "Summary log written to: $LOG"
