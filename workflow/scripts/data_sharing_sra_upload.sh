#!/bin/bash

# Script to upload FASTQ files to SRA via FTP
#
# Usage:
#   ./upload_sra.sh <CONFIG_FILE> <PLATE>
#
# Arguments:
#   CONFIG_FILE   Path to a Bash-compatible config file (see format below)
#   PLATE         Plate identifier (used to construct paths)
#
# ------------------------------------------------------------------
# CONFIG FILE FORMAT (IMPORTANT)
# ------------------------------------------------------------------
# The config file MUST be a plain text file containing valid Bash
# variable assignments. It is sourced at runtime using `source`.
# <USER>, <PASS> and <REMOTE_BASE> are provided by SRA.
#
# Example (../config/passwords/sra_upload.sh):
#
#   USER="username"
#   PASS="your_password_here"
#   REMOTE_BASE="uploads/your_upload_dir"
#   ROOT_BASE="/path/to/local/data"
#   LOG="/path/to/logfile.log"
#   MAX_JOBS=3
#
# ------------------------------------------------------------------

set -euo pipefail

# --- Input arguments ---
CONFIG_FILE="${1:-}"
PLATE="${2:-}"

if [[ -z "$CONFIG_FILE" || -z "$PLATE" ]]; then
    echo "Usage: $0 <CONFIG_FILE> <PLATE>"
    exit 1
fi

# --- Validate config file ---
if [[ ! -f "$CONFIG_FILE" ]]; then
    echo "ERROR: Config file not found: $CONFIG_FILE"
    exit 1
fi

# --- Load config ---
source "$CONFIG_FILE"

# --- Validate required vars ---
for var in USER PASS REMOTE_BASE ROOT_BASE LOG MAX_JOBS; do
    if [[ -z "${!var:-}" ]]; then
        echo "ERROR: Variable $var is not set in config"
        exit 1
    fi
done

# --- Derived paths ---
REMOTE_PATH="${REMOTE_BASE}/${PLATE}"
LOCAL_DIR="${ROOT_BASE}/${PLATE}"

# --- Execution ---
cd "$LOCAL_DIR" || { echo "ERROR: Cannot cd to $LOCAL_DIR"; exit 1; }

echo "--- Starting Upload $(date) ---" >> "$LOG"

for f in *.fastq.gz; do
    [[ -e "$f" ]] || continue

    (
        echo "[$(date +%T)] START: $f" >> "$LOG"

        curl -s -u "$USER:$PASS" \
             -T "$f" \
             -C - \
             --ftp-pasv \
             "ftp://ftp-private.ncbi.nlm.nih.gov/$REMOTE_PATH/"

        if [[ $? -eq 0 ]]; then
            echo "[$(date +%T)] SUCCESS: $f" >> "$LOG"
        else
            echo "[$(date +%T)] FAIL: $f" >> "$LOG"
        fi
    ) &

    # --- Limit parallel jobs ---
    while [[ $(jobs -r -p | wc -l) -ge "$MAX_JOBS" ]]; do
        sleep 10
    done
done

wait
echo "--- Finished Upload $(date) ---" >> "$LOG"
