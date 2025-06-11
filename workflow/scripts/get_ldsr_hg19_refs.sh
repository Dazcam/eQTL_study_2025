#!/bin/bash
# Run from workflow dir 

# Set variables
LDSR_ROOT=../resources/ldsr/
LDSR_DIR=../resources/ldsr/reference_files
WEB_ROOT="https://zenodo.org/records/10515792/files/"

# Clone repo
git clone git@github.com:bulik/ldsc.git ${LDSR_ROOT}

# Create dir for LDSR support files
mkdir -p ${LDSR_DIR}

# Declare associative array for expected extraction names
declare -A extract_dirs=(
    ["1000G_Phase3_baseline_v1.2_ldscores"]="baseline_v1.2"
    ["1000G_Phase3_frq"]="1000G_Phase3_frq"
    ["1000G_Phase3_ldscores"]="LDscore"
    ["1000G_Phase3_plinkfiles"]="1000G_EUR_Phase3_plink"
    ["1000G_Phase3_weights_hm3_no_MHC"]="1000G_Phase3_weights_hm3_no_MHC"
)

# Download and extract files in a loop
for file in "${!extract_dirs[@]}"; do
    tgz_file="${file}.tgz"
    wget "${WEB_ROOT}${tgz_file}?download=1" -O "${LDSR_DIR}/${tgz_file}"
  
    # Extract into LDSR_DIR
    tar -zxf "${LDSR_DIR}/${tgz_file}" -C "${LDSR_DIR}"
    extracted_path="${LDSR_DIR}/${extract_dirs[$file]}"
    target_path="${LDSR_DIR}/${file}"

    # Check if extracted folder is different from the target path before moving
    if [[ "$extracted_path" != "$target_path" ]]; then
        mv "$extracted_path" "$target_path"
    fi
    
    # Remove the downloaded .tgz fil
    rm -rf "${LDSR_DIR}/${tgz_file}"
done

# Download hm3_no_MHC.list.txt separately
wget "${WEB_ROOT}hm3_no_MHC.list.txt?download=1" -O "${LDSR_DIR}/hm3_no_MHC.list.txt"
