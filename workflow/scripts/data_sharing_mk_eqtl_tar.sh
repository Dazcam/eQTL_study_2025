#!/bin/bash

set -o pipefail

OUTDIR=$1
THREADS=$2
CELL_TYPES_STR=$3
COMPRESSION=6

echo "[INFO] Creating output directory: ${OUTDIR}"
mkdir -p "${OUTDIR}"
OUTFILE=$(realpath "${OUTDIR%/}/nbray_2026_dev_brain_eqtl.tar.gz")

if [[ -z "$THREADS" || ! "$THREADS" =~ ^[0-9]+$ ]]; then
    echo "[ERROR] THREADS must be a positive integer, got: '$THREADS'" >&2
    exit 2
fi

if [[ -z "$CELL_TYPES_STR" ]]; then
    echo "[ERROR] No cell types provided" >&2
    exit 2
fi

read -r -a CELL_TYPES <<< "$CELL_TYPES_STR"
echo "[INFO] Parsed cell types (${#CELL_TYPES[@]}): ${CELL_TYPES[*]}" >&2

# Create tarball
echo "[INFO] Creating tarball..."
cd "${OUTDIR}" || exit 1

files=()
for ct in "${CELL_TYPES[@]}"; do
  files+=("${ct}_cis_eQTL_perm.tsv.gz" "${ct}_cis_eQTL_nominal.tsv.gz")
done

expected_count=${#files[@]}

echo "[INFO] Creating tarball and compressing with pigz ..."

if ! tar -cf - "${files[@]}" | pigz -${COMPRESSION} -p "${THREADS}" > "${OUTFILE}"; then
    echo "[ERROR] Tar or Pigz failed during creation!"
    exit 1
fi

# Verification step
echo "[INFO] Verifying archive ..."

# Gzip integrity
if ! pigz -t "${OUTFILE}"; then
    echo "[ERROR] Gzip stream is corrupt!"
    exit 1
fi

# Count verification
actual_count=$(tar -tf "${OUTFILE}" | wc -l)

if [ "$actual_count" -eq "$expected_count" ]; then
    echo "[SUCCESS] Integrity check passed ($actual_count/$expected_count files present)."
else
    echo "[ERROR] File count mismatch! Expected $expected_count, found $actual_count."
    echo "[DEBUG] Files actually in archive:"
    tar -tf "${OUTFILE}"
    exit 1
fi

echo "[INFO] Job completed successfully."
