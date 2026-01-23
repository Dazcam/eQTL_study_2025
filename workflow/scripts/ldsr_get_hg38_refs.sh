#!/bin/bash

# Script to download and unpack LDSC hg38 reference files from Zenodo
# Allows specifying output directory as a parameter (defaults to LDSC_hg38)
# Moves all files to the parent directory and removes GRCh38 directory and GRCh38.tgz

# Set output directory (default: LDSC_hg38 if no parameter provided)
OUTPUT_DIR="${1:-LDSC_hg38}"

# Create output directory
mkdir -p "$OUTPUT_DIR"
cd "$OUTPUT_DIR" || exit 1

# Download files from Zenodo
wget -O GRCh38.tgz "https://zenodo.org/records/7768714/files/GRCh38.tgz?download=1"
wget -O hm3_no_MHC.list.txt "https://zenodo.org/records/10515792/files/hm3_no_MHC.list.txt?download=1"

# Check if downloads were successful
if [[ ! -f GRCh38.tgz || ! -f hm3_no_MHC.list.txt ]]; then
  echo "Error: One or more downloads failed."
  exit 1
fi

# Extract GRCh38.tgz
tar -xvzf GRCh38.tgz
if [[ $? -ne 0 ]]; then
  echo "Error: Failed to extract GRCh38.tgz"
  exit 1
fi

# Move into GRCh38 directory
cd GRCh38 || exit 1

# Extract nested tar.gz files
for file in baseline_v1.2.tgz plink_files.tgz weights.tgz; do
  if [[ -f "$file" ]]; then
    tar -xvzf "$file"
    if [[ $? -ne 0 ]]; then
      echo "Error: Failed to extract $file"
      exit 1
    fi
  else
    echo "Warning: $file not found, skipping extraction"
  fi
done

# Move all files from GRCh38 to parent directory (OUTPUT_DIR)
mv baseline_v1.2 plink_files weights ../
if [[ $? -ne 0 ]]; then
  echo "Error: Failed to move files to $OUTPUT_DIR"
  exit 1
fi

# Move back to parent directory
cd ..

# Remove GRCh38 directory and GRCh38.tgz
rm -rf GRCh38 GRCh38.tgz
if [[ $? -ne 0 ]]; then
  echo "Error: Failed to remove GRCh38 directory or GRCh38.tgz"
  exit 1
fi

# Clean up nested tar.gz files (if any remain)
rm -f baseline_v1.2.tgz plink_files.tgz weights.tgz

echo "Download and extraction complete. Files are in $(pwd)"

