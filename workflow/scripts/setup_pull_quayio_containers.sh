#!/usr/bin/env bash
#
# Pull public Quay.io Apptainer containers
# Run in repo root dir. Saves .sif files to 
# Requires: apptainer (or singularity) installed

set -euo pipefail  # Exit on error, undefined vars, pipe failures

# Directory to save the pulled containers
OUTPUT_DIR="resources/containers"
mkdir -p "$OUTPUT_DIR"

# List of Quay.io repositories
REPOS=(
    "dazcam1/genotype-qc2hrc"
    "dazcam1/ctwas"
    "dazcam1/r_eqtl"        
    "dazcam1/twas"
    "dazcam1/susier"
)

# Tag to pull (change if you used something other than latest)
TAG="latest"

echo "Starting pull of ${#REPOS[@]} containers from quay.io"
echo "Saving to: $OUTPUT_DIR"
echo ""

for repo in "${REPOS[@]}"; do
    # Convert repo name to safe filename 
    filename="${repo#*/}"          # remove dazcam1/
    filename="${filename}_$TAG.sif"

    full_uri="oras://quay.io/$repo:$TAG"

    echo "────────────────────────────────────────"
    echo "Pulling: $full_uri"
    echo "Saving as: $OUTPUT_DIR/$filename"

    if apptainer pull "$OUTPUT_DIR/$filename" "$full_uri"; then
        echo "Success: $filename pulled"
        ls -lh "$OUTPUT_DIR/$filename" | awk '{print "Size: " $5}'
    else
        echo "Error pulling $full_uri — check connection, auth (for public should be fine), or tag"
    fi
    echo ""
done

echo "────────────────────────────────────────"
echo "All pulls completed."
echo "Check files in $OUTPUT_DIR/"
ls -lh "$OUTPUT_DIR/" | grep ".sif"