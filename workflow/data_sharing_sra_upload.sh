#!/bin/bash

# --- Config ---
USER="subftp"
PASS="REMOVED_SECRET"
REMOTE_PATH="uploads/BrayN3_cardiff.ac.uk_OymvBIkH/plate3"
LOCAL_DIR="/scratch/SCWF00021/eQTL_study_2025/results/14DATA_SHARING/plate3"
LOG="/scratch/SCWF00021/eQTL_study_2025/workflow/logs/sra_upload_progress.log"
MAX_JOBS=3

# --- Execution ---
cd "$LOCAL_DIR" || exit 1

echo "--- Starting Upload $(date) ---" >> "$LOG"

for f in *.fastq.gz; do
    # 1. Background the upload
    (
        echo "[$(date +%T)] START: $f" >> "$LOG"
        
        # -C - is vital for 4TB; it resumes if the connection drops
        curl -s -u "$USER:$PASS" \
             -T "$f" \
             -C - \
             --ftp-pasv \
             "ftp://ftp-private.ncbi.nlm.nih.gov/$REMOTE_PATH/"
        
        if [ $? -eq 0 ]; then
            echo "[$(date +%T)] SUCCESS: $f" >> "$LOG"
        else
            echo "[$(date +%T)] FAIL: $f" >> "$LOG"
        fi
    ) &

    # 2. Limit the number of parallel curls
    while [[ $(jobs -r -p | wc -l) -ge $MAX_JOBS ]]; do
        sleep 10
    done
done

wait
echo "--- Finished Upload $(date) ---" >> "$LOG"
