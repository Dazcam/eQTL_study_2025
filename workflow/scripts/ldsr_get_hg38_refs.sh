#!/bin/bash
#--------------------------------------------------------------------------------------
#
#    Download LDSC hg38 reference files
#
#--------------------------------------------------------------------------------------
#
# Pipeline:  09SLDSR | Rule: get_hg38_refs
#            Upstream:   none (downloads from Zenodo)
#            Downstream: make_annot, lift_hapmap3_snps, ldsr_ld_scores_hg38,
#                        ldsr_strat_hg38_bl_v12
#
# Purpose:   - Download the GRCh38 LDSC reference bundle (Zenodo 7768714) and the
#              HapMap3 no-MHC SNP list (Zenodo 10515792)
#            - Unpack baseline_v1.2, plink_files and weights into the output dir
#            - Remove the downloaded archives
#
# Usage:     ldsr_get_hg38_refs.sh [output_dir]
#
# Inputs:    output_dir  Directory for the reference files (default: LDSC_hg38)
#
# Outputs:   <output_dir>/baseline_v1.2/  baseline v1.2 LD scores per chr
#            <output_dir>/plink_files/    1000G EUR hg38 plink and .frq files per chr
#            <output_dir>/weights/        HapMap3 no-MHC regression weights per chr
#            <output_dir>/hm3_no_MHC.list.txt
#
# Notes:     The SNP list is kept so the hg38 refs are complete on their own; the
#            liftover rule currently uses the hg19 copy.
#
#--------------------------------------------------------------------------------------

set -euo pipefail

OUTPUT_DIR="${1:-LDSC_hg38}"

mkdir -p "${OUTPUT_DIR}"
cd "${OUTPUT_DIR}"

##  Download  -------------------------------------------------------------------------
wget -O GRCh38.tgz "https://zenodo.org/records/7768714/files/GRCh38.tgz?download=1"
wget -O hm3_no_MHC.list.txt "https://zenodo.org/records/10515792/files/hm3_no_MHC.list.txt?download=1"

##  Extract  --------------------------------------------------------------------------
tar -xvzf GRCh38.tgz || { echo "Error: Failed to extract GRCh38.tgz"; exit 1; }

cd GRCh38
for file in baseline_v1.2.tgz plink_files.tgz weights.tgz; do
  if [[ -f "${file}" ]]; then
    tar -xvzf "${file}" || { echo "Error: Failed to extract ${file}"; exit 1; }
  else
    echo "Warning: ${file} not found, skipping extraction"
  fi
done

mv baseline_v1.2 plink_files weights ../ ||
  { echo "Error: Failed to move files to ${OUTPUT_DIR}"; exit 1; }
cd ..

##  Tidy  -----------------------------------------------------------------------------
rm -rf GRCh38 GRCh38.tgz
rm -f baseline_v1.2.tgz plink_files.tgz weights.tgz

echo "Download and extraction complete. Files are in $(pwd)"

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
