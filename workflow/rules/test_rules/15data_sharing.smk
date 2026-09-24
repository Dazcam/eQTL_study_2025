# Extra smk pipeline for SRA upload
configfile: "../config/config.yaml"
localrules: extract_fq

import json
import os

PLATES = config["plate"]

rule all:
    input:
        expand("../results/14DATA_SHARING/{plate}/local_checksums_{plate}.txt", plate=PLATES)
#        [f"../results/14DATA_SHARING/{plate}/.copy_complete" for plate in PLATES]

rule extract_fq:
    """
    This rule reads the JSON for a plate, extracts all file paths, 
    and uses rsync to copy them while preserving the run-folder name
    to prevent filename collisions.
    """
    input:  "../config/samples_{plate}.json"
    output: "../results/14DATA_SHARING/{plate}/.copy_complete"
    params: out_dir = "../results/14DATA_SHARING/{plate}/"
    log:    "logs/{plate}.rsync.log"
    shell:
        r"""
        echo "--- Starting Plate {wildcards.plate} Transfer $(date) ---" > {log}
        module load parallel
        
        # 1. Get the list of files (the sed part creates the ./ relative point)
        FILE_LIST=$(jq -r '.. | .R1?, .R2? | arrays[]' {input} | sed 's|\([^/]*\)/[^/]*$|./\0|')

        # 2. Run parallel
        echo "$FILE_LIST" | parallel --no-notice -j 3 \
            "echo '[$(date +%T)] START: {{}}' >> {log} && \
             rsync -aqR {{}} {params.out_dir} && \
             echo '[$(date +%T)] SUCCESS: {{}}' >> {log} || \
             echo '[$(date +%T)] FAIL: {{}}' >> {log}"

        # 3. Summary logic stays the same...
        touch {output}
        """

rule create_json:
    output: "../config/rename_fastq.json"
    params: json_dir = "../config"
    log:    "logs/create_json.log"
    shell:  """
	    python scripts/data_sharing_rename_fqs_for_sra.py \
              --json-dir {params.json_dir} \
              --plates plate1 plate2 plate3 \
              --output {output} >> {log}
            """

rule rename_fastqs:
    """
    Rename fastq files for SRA submission using a dedicated python script.
    """
    input:  json = "../config/rename_fastq.json",
            copy_done = "../results/14DATA_SHARING/{plate}/.copy_complete"
    output: touch("../results/14DATA_SHARING/{plate}/.rename_complete")
    params: plate_dir = "../results/14DATA_SHARING/{plate}/"
    log:    "logs/{plate}.rename.log"
    script: "../scripts/data_sharing_rename_fqs_for_sra.py"

rule create_checksum:
    input:  "../results/14DATA_SHARING/{plate}/.rename_complete"
    output: "../results/14DATA_SHARING/{plate}/local_checksums_{plate}.txt"
    params: outdir = "../results/14DATA_SHARING/{plate}/"
    log:    "logs/{plate}.rename.log"
    shell:  "md5sum {params.outdir}*.fastq.gz > {output}"
