configfile: "../config/config.yaml"

from pathlib import Path
from snakemake.io import glob_wildcards

FASTQ_DIR = "/nfs/neurocluster/sequencers2/Nick_Bray/250716_A00748_0722_AHC52WDSXF_fastq"


# Use all 4 wildcards to find the files
files = glob_wildcards(f"{FASTQ_DIR}/{{sample}}_S{{sample_id}}_{{lane}}_R{{read}}_001.fastq.gz")

# Filter out Undetermined and just keep sample + lane
valid_samples = [
    (s, l) for s, l, r in zip(files.sample, files.lane, files.read)
    if s != "Undetermined" and r == "1"
]

# De-duplicate
ALL_SAMPLE_LANES = sorted(set(valid_samples))
print(ALL_SAMPLE_LANES)

# Separate into lists if needed
ALL_SAMPLES = sorted(set(s for s, l in ALL_SAMPLE_LANES))
LANES = sorted(set(l for s, l in ALL_SAMPLE_LANES))


#LANES = ["L001", "L002", "L003", "L004"]

localrules: copy_fqs

rule copy_fqs:
    input:
        r1 = lambda wildcards: sorted(Path(FASTQ_DIR).glob(f"{wildcards.sample}_S*_{wildcards.lane}_R1_001.fastq.gz"))[0],
        r2 = lambda wildcards: sorted(Path(FASTQ_DIR).glob(f"{wildcards.sample}_S*_{wildcards.lane}_R2_001.fastq.gz"))[0]
    output:
        r1 = temp("fastqs/250716/{sample}_{lane}_R1_001.fastq.gz"),
        r2 = temp("fastqs/250716/{sample}_{lane}_R2_001.fastq.gz")
    shell:
        """
        cp {input.r1} {output.r1}
        cp {input.r2} {output.r2}
        """

rule check_read_counts:
    input:
        r1 = "fastqs/250716/{sample}_{lane}_R1_001.fastq.gz",
        r2 = "fastqs/250716/{sample}_{lane}_R2_001.fastq.gz"
    output:
        "reports/fq_check/read_counts/{sample}_{lane}_read_counts.txt"
    priority: 50
    resources:
        mem_mb = 8000,
        time = "02:00:00",
        threads = 8
    shell:
        """
        echo "Sample: {wildcards.sample}" > {output}
        echo "Lane\tR1_File\tR1_Count\tR2_File\tR2_Count\tStatus" >> {output}

        r1={input.r1}
        r2={input.r2}
        lane={wildcards.lane}
        r1_count=$(zcat $r1 | grep -c "^@")
        r2_count=$(zcat $r2 | grep -c "^@")

        if [ "$r1_count" -eq "$r2_count" ]; then
            status="Pass"
        else
            status="Fail"
        fi

        echo "$lane\t$(basename $r1)\t$r1_count\t$(basename $r2)\t$r2_count\t$status" >> {output}
        """

rule check_fastq_integrity:
    input:
        r1 = "fastqs/250716/{sample}_{lane}_R1_001.fastq.gz",
        r2 = "fastqs/250716/{sample}_{lane}_R2_001.fastq.gz"
    output:
        "reports/fq_check/integrity/{sample}_{lane}_integrity.txt"
    priority: 50
    shell:
        """
        echo -n "{input.r1}: " > {output}
        gunzip -t {input.r1} && echo "OK" || echo "Corrupted" >> {output}
        echo -n "{input.r2}: " >> {output}
        gunzip -t {input.r2} && echo "OK" || echo "Corrupted" >> {output}
        """



rule aggregate_read_counts:
    input:
        integrity=expand("reports/fq_check/integrity/{sample}_{lane}_integrity.txt", sample=ALL_SAMPLES, lane=["L001", "L002", "L003", "L004"]),
        reports=expand("reports/fq_check/read_counts/{sample}_{lane}_read_counts.txt", sample=ALL_SAMPLES, lane=["L001", "L002", "L003", "L004"])
    output:
        summary="reports/read_counts/summary_read_counts.txt"
    shell:
        """
        mkdir -p reports/read_counts
        echo "Summary of read counts for 250716_A00748_0722_AHC52WDSXF_fastq" > {output.summary}
        echo "Sample\tLane\tR1_File\tR1_Count\tR2_File\tR2_Count\tStatus" >> {output.summary}
        for r in {input.reports}; do
            tail -n +3 $r >> {output.summary}
        done
        """
