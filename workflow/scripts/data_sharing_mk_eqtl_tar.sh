#!/bin/bash
#--------------------------------------------------------------------------------------
#
#    Tarball of eQTL permutation and nominal files for data sharing
#
#--------------------------------------------------------------------------------------
#
# Pipeline:  18DATA-SHARING | Rule: mk_eqtl_tar
#            Upstream:   data_perm, data_nominal
#            Downstream: final output (figshare)
#
# Purpose:   - Archive the given files (all in one directory) under their base names,
#              in the order given
#            - Compress with pigz, then check the gzip stream and the file count
#
# Usage:     data_sharing_mk_eqtl_tar.sh <out_file> <threads> <file> [<file> ...]
#
# Outputs:   out_file     .tar.gz of the files, no leading directory
#
#--------------------------------------------------------------------------------------

set -euo pipefail

OUTFILE=$1
THREADS=$2
shift 2
FILES=("$@")
COMPRESSION=6

if [[ ! "$THREADS" =~ ^[0-9]+$ ]]; then
    echo "[ERROR] THREADS must be a positive integer, got: '$THREADS'" >&2
    exit 2
fi

if [[ ${#FILES[@]} -eq 0 ]]; then
    echo "[ERROR] No files provided" >&2
    exit 2
fi

IN_DIR=$(dirname "${FILES[0]}")
names=()
for f in "${FILES[@]}"; do
    if [[ "$(dirname "$f")" != "$IN_DIR" ]]; then
        echo "[ERROR] All files must be in ${IN_DIR}: $f" >&2
        exit 2
    fi
    names+=("$(basename "$f")")
done

mkdir -p "$(dirname "$OUTFILE")"
OUTFILE=$(realpath "$OUTFILE")
expected_count=${#names[@]}

echo "[INFO] Archiving ${expected_count} files from ${IN_DIR} to ${OUTFILE} (${THREADS} threads)"
echo "[INFO] Files: ${names[*]}"

##  Create tarball  -------------------------------------------------------------------
if ! tar -C "$IN_DIR" -cf - "${names[@]}" | pigz -${COMPRESSION} -p "$THREADS" > "$OUTFILE"; then
    echo "[ERROR] tar or pigz failed"
    exit 1
fi

##  Verify  ---------------------------------------------------------------------------
echo "[INFO] Verifying archive ..."
if ! pigz -t "$OUTFILE"; then
    echo "[ERROR] Gzip stream is corrupt"
    exit 1
fi

actual_count=$(tar -tf "$OUTFILE" | wc -l)
if [[ "$actual_count" -eq "$expected_count" ]]; then
    echo "[SUCCESS] Integrity check passed (${actual_count}/${expected_count} files present)"
else
    echo "[ERROR] File count mismatch: expected ${expected_count}, found ${actual_count}"
    tar -tf "$OUTFILE"
    exit 1
fi

echo "[INFO] Done."

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
