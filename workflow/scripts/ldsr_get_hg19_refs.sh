#!/bin/bash
# Run from workflow dir

# Set variables
LDSR_ROOT=../resources/ldsr/
LDSR_DIR=../resources/ldsr/ldsr_hg19_refs
WEB_ROOT="https://zenodo.org/records/10515792/files/"

# Clone repo
#git clone git@github.com:bulik/ldsc.git ${LDSR_ROOT}

# Create dir for LDSR support files
mkdir -p ${LDSR_DIR}

# Declare associative array for expected extraction names
declare -A extract_dirs=(
    ["1000G_Phase3_baseline_v1.2_ldscores"]="baseline_v1.2"
    ["1000G_Phase3_ldscores"]="ldscores"
    ["1000G_Phase3_frq"]="frq"
    ["1000G_Phase3_plinkfiles"]="plink_files"
    ["1000G_Phase3_weights_hm3_no_MHC"]="weights"
)

# Download and extract files in a loop
for file in "${!extract_dirs[@]}"; do
    tgz_file="${file}.tgz"
    url="${WEB_ROOT}${tgz_file}?download=1"
    output_path="${LDSR_DIR}/${tgz_file}"

    echo "Downloading ${url} -> ${output_path}"
    wget "${url}" -O "${output_path}"

    # Extract into LDSR_DIR
    echo "Extracting ${output_path} to ${LDSR_DIR}"
    tar -zxf "${output_path}" -C "${LDSR_DIR}"

    # Set extracted_path, with a fallback for the plinkfiles case
    extracted_path="${LDSR_DIR}/${file}"

    if [[ ! -d "$extracted_path" ]]; then
        # Special case: for plinkfiles, folder is named without 'files'
        if [[ "$file" == "1000G_Phase3_plinkfiles" ]]; then
            extracted_path="${LDSR_DIR}/1000G_EUR_Phase3_plink"
        fi
    fi

    target_path="${LDSR_DIR}/${extract_dirs[$file]}"

    # Rename extracted folder to target name
    if [[ "$extracted_path" != "$target_path" ]]; then
        echo "Moving ${extracted_path} -> ${target_path}"
        mv "$extracted_path" "$target_path"
    fi

    echo "Removing ${output_path}"
    rm -rf "${output_path}"
done

# Move frq files into plink dir
echo "Moving frequency files into plink directory"
mv "${LDSR_DIR}/frq/"* "${LDSR_DIR}/plink_files/"

# Download hm3_no_MHC.list.txt separately
hm3_url="${WEB_ROOT}hm3_no_MHC.list.txt?download=1"
hm3_out="${LDSR_DIR}/hm3_no_MHC.list.txt"
echo "Downloading ${hm3_url} -> ${hm3_out}"
wget "${hm3_url}" -O "${hm3_out}"

#Tidy
rm -rf "${LDSR_DIR}/frq/"
rm -rf "${LDSR_DIR}/baseline_v1.2/1000G_Phase3_baseline_v1.2_ldscores.tgz"
rm -rf "${LDSR_DIR}/weights/1000G_Phase3_weights_hm3_no_MHC"
