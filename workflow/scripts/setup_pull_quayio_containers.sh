#!/usr/bin/env bash
#
# Pull multiple Apptainer/Singularity containers (.sif files)
# - Skips if file already exists
# - Supports oras:// and docker:// (auto-converts to .sif)
# - Saves in ./containers/

set -euo pipefail

# Check for apptainer or singularity
if command -v apptainer &> /dev/null; then
    CMD="apptainer"
elif command -v singularity &> /dev/null; then
    CMD="singularity"
else
    echo "ERROR: Neither 'apptainer' nor 'singularity' command found."
    echo "       Please install Apptainer or Singularity."
    exit 1
fi

# Check if apptainer is available
if ! command -v apptainer &> /dev/null; then
    echo "ERROR: 'apptainer' command not found. Install Apptainer or replace with 'singularity'."
    exit 1
fi

OUTPUT_DIR="./containers"
mkdir -p "$OUTPUT_DIR"

# Your Quay.io ORAS containers
declare -A QUAY_CONTAINERS=(
    ["genotype-qc2hrc"]="oras://quay.io/dazcam1/genotype-qc2hrc:latest"
    ["r_eqtl"]="oras://quay.io/dazcam1/r_eqtl:latest"
    ["twas"]="oras://quay.io/dazcam1/twas:latest"
    ["susier"]="oras://quay.io/dazcam1/susier:latest"
)

# External Docker containers (converted to .sif)
declare -A EXTERNAL_CONTAINERS=(
    ["ubuntu-22.04"]="docker://ubuntu:22.04"
    ["tensorqtl"]="docker://francois4/tensorqtl:latest"
)

# Combine both
declare -A ALL_CONTAINERS
for k in "${!QUAY_CONTAINERS[@]}"; do ALL_CONTAINERS["$k"]="${QUAY_CONTAINERS[$k]}"; done
for k in "${!EXTERNAL_CONTAINERS[@]}"; do ALL_CONTAINERS["$k"]="${EXTERNAL_CONTAINERS[$k]}"; done

echo "Checking/pulling ${#ALL_CONTAINERS[@]} containers → $OUTPUT_DIR"
echo ""

for name in "${!ALL_CONTAINERS[@]}"; do
    uri="${ALL_CONTAINERS[$name]}"
    filename="${name}.sif"
    target="$OUTPUT_DIR/$filename"

    if [ -f "$target" ]; then
        echo "SKIP: $filename already exists"
        ls -lh "$target" | awk '{print "  Size: " $5}'
        continue
    fi

    echo "────────────────────────────────────────"
    echo "Pulling: $uri"
    echo "Saving as: $target"

    if apptainer pull "$target" "$uri"; then
        echo "Success"
        ls -lh "$target" | awk '{print "Size: " $5}'
    else
        echo "ERROR pulling $uri"
    fi
    echo ""
done

echo "Done. Files in $OUTPUT_DIR:"
ls -lh "$OUTPUT_DIR/" | grep ".sif" || echo "No .sif files found."
