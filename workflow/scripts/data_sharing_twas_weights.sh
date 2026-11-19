#!/bin/bash
#--------------------------------------------------------------------------------------
#
#    Tarball of FUSION TWAS weights for data sharing
#
#--------------------------------------------------------------------------------------
#
# Pipeline:  18DATA-SHARING | Rule: data_weights
#            Upstream:   14CTWAS copy_fusion_weights (weights_dir/<pipeline label>/)
#            Downstream: final output (figshare)
#
# Purpose:   - Archive weights_dir/<pipeline label>/ for each cell type, stored in the
#              archive as <archive_root>/<manuscript label>/
#            - Compress with pigz, then check the gzip stream and the tar listing
#
# Usage:     data_sharing_twas_weights.sh <out_file> <threads> <weights_dir>
#                <archive_root> <pipeline_label>=<manuscript_label> [...]
#
# Outputs:   out_file     .tar.gz with one directory per cell type
#
#--------------------------------------------------------------------------------------

set -euo pipefail

OUTFILE=$1
THREADS=$2
WEIGHTS_DIR=${3%/}
ARCHIVE_ROOT=$4
shift 4
PAIRS=("$@")
COMPRESSION=6

if [[ ! "$THREADS" =~ ^[0-9]+$ ]]; then
    echo "[ERROR] THREADS must be a positive integer, got: '$THREADS'" >&2
    exit 2
fi

if [[ ${#PAIRS[@]} -eq 0 ]]; then
    echo "[ERROR] No cell types provided" >&2
    exit 2
fi

BASE_DIR=$(dirname "$WEIGHTS_DIR")
SOURCE_NAME=$(basename "$WEIGHTS_DIR")

##  Build tar arguments: one directory and one name transform per cell type  ----------
tar_dirs=()
transforms=()
for pair in "${PAIRS[@]}"; do
    old=${pair%%=*}
    new=${pair#*=}
    if [[ ! "$old" =~ ^[A-Za-z0-9_-]+$ || ! "$new" =~ ^[A-Za-z0-9_-]+$ ]]; then
        echo "[ERROR] Unexpected characters in cell type pair: $pair" >&2
        exit 2
    fi
    if [[ ! -d "${WEIGHTS_DIR}/${old}" ]]; then
        echo "[ERROR] No weights directory: ${WEIGHTS_DIR}/${old}" >&2
        exit 2
    fi
    tar_dirs+=("${SOURCE_NAME}/${old}")
    transforms+=("--transform=s|^${SOURCE_NAME}/${old}\(/\|$\)|${ARCHIVE_ROOT}/${new}\1|")
    echo "[INFO] ${WEIGHTS_DIR}/${old} -> ${ARCHIVE_ROOT}/${new}"
done

mkdir -p "$(dirname "$OUTFILE")"

##  Create tarball  -------------------------------------------------------------------
echo "[INFO] Creating ${OUTFILE} from ${#tar_dirs[@]} cell types (${THREADS} threads)"
tar -C "$BASE_DIR" "${transforms[@]}" -cf - "${tar_dirs[@]}" \
    | pigz -${COMPRESSION} -p "$THREADS" > "$OUTFILE"

##  Verify  ---------------------------------------------------------------------------
echo "[INFO] Verifying gzip integrity ..."
pigz -t "$OUTFILE"

echo "[INFO] Verifying tar structure ..."
n_files=$(tar -tzf "$OUTFILE" | wc -l)
echo "[SUCCESS] Archive created and verified (${n_files} entries)"

echo "[INFO] Done."

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
