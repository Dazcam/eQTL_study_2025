#!/bin/bash
#--------------------------------------------------------------------------------------
#
#    Download LDSC hg19 reference files
#
#--------------------------------------------------------------------------------------
#
# Pipeline:  09SLDSR | Rule: get_hg19_refs
#            Upstream:   none (downloads from Zenodo)
#            Downstream: lift_hapmap3_snps
#
# Purpose:   - Download the 1000G Phase 3 hg19 LDSC reference files (Zenodo 10515792)
#            - Unpack and rename to baseline_v1.2, ldscores, plink_files and weights
#            - Move the .frq files into plink_files
#            - Download the HapMap3 no-MHC SNP list
#
# Usage:     ldsr_get_hg19_refs.sh [output_dir]   (run from workflow/)
#
# Inputs:    output_dir  Directory for the reference files
#                        (default: ../resources/ldsr/ldsr_hg19_refs)
#
# Outputs:   <output_dir>/baseline_v1.2/  baseline v1.2 LD scores per chr
#            <output_dir>/ldscores/       1000G Phase 3 LD scores per chr
#            <output_dir>/plink_files/    1000G EUR plink and .frq files per chr
#            <output_dir>/weights/        HapMap3 no-MHC regression weights per chr
#            <output_dir>/hm3_no_MHC.list.txt
#
# Notes:     ldsc itself is not downloaded here. It is cloned once by hand:
#              git clone git@github.com:bulik/ldsc.git ../resources/ldsr/ldsc
#
#--------------------------------------------------------------------------------------

set -euo pipefail

LDSR_DIR="${1:-../resources/ldsr/ldsr_hg19_refs}"
LDSR_DIR="${LDSR_DIR%/}"
WEB_ROOT="https://zenodo.org/records/10515792/files/"

mkdir -p "${LDSR_DIR}"

# Archive name -> target directory name
declare -A extract_dirs=(
    ["1000G_Phase3_baseline_v1.2_ldscores"]="baseline_v1.2"
    ["1000G_Phase3_ldscores"]="ldscores"
    ["1000G_Phase3_frq"]="frq"
    ["1000G_Phase3_plinkfiles"]="plink_files"
    ["1000G_Phase3_weights_hm3_no_MHC"]="weights"
)

##  Download and extract  -------------------------------------------------------------
for file in "${!extract_dirs[@]}"; do
    tgz_file="${file}.tgz"
    url="${WEB_ROOT}${tgz_file}?download=1"
    output_path="${LDSR_DIR}/${tgz_file}"

    echo "Downloading ${url} -> ${output_path}"
    wget "${url}" -O "${output_path}"

    echo "Extracting ${output_path} to ${LDSR_DIR}"
    tar -zxf "${output_path}" -C "${LDSR_DIR}"

    # The plinkfiles archive extracts to 1000G_EUR_Phase3_plink
    extracted_path="${LDSR_DIR}/${file}"
    if [[ ! -d "${extracted_path}" && "${file}" == "1000G_Phase3_plinkfiles" ]]; then
        extracted_path="${LDSR_DIR}/1000G_EUR_Phase3_plink"
    fi

    target_path="${LDSR_DIR}/${extract_dirs[$file]}"
    if [[ "${extracted_path}" != "${target_path}" ]]; then
        echo "Moving ${extracted_path} -> ${target_path}"
        mv "${extracted_path}" "${target_path}"
    fi

    echo "Removing ${output_path}"
    rm -rf "${output_path}"
done

echo "Moving frequency files into plink directory"
mv "${LDSR_DIR}/frq/"* "${LDSR_DIR}/plink_files/"

##  Download SNP list  ----------------------------------------------------------------
hm3_url="${WEB_ROOT}hm3_no_MHC.list.txt?download=1"
hm3_out="${LDSR_DIR}/hm3_no_MHC.list.txt"
echo "Downloading ${hm3_url} -> ${hm3_out}"
wget "${hm3_url}" -O "${hm3_out}"

##  Tidy  -----------------------------------------------------------------------------
rm -rf "${LDSR_DIR}/frq/"
rm -rf "${LDSR_DIR}/baseline_v1.2/1000G_Phase3_baseline_v1.2_ldscores.tgz"
rm -rf "${LDSR_DIR}/weights/1000G_Phase3_weights_hm3_no_MHC"

echo "Download and extraction complete. Files are in ${LDSR_DIR}"

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
