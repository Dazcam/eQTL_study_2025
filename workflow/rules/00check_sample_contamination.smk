import json
configfile: '../config/config.yaml'

SUBLIBS = json.load(open(config['BAM_FILES']))
CHROM = "chr22"  # Adjust to "22" if BAMs lack "chr" prefix
COMMON_VCF = "../results/04GENOTYPES/GenotypeQCtoHRC/dragondata_extra/imputation-sites/ALL.TOPMed_freeze5_hg38_dbSNP.vcf.gz"
REF_VCF = "../results/04GENOTYPES/GenotypeQCtoHRC/eqtlgenotypeshg38pretm/imputeMe/eqtl_genotypes_hg19.gh.topmed.chr22.vcf.gz"  # Your reference genotypes for chr22

# Note that I'm using preimputed genotypes here after running GenotypeQCtoHRC

localrules: copy_bam, make_whitelist

SCRATCH_BAM_DIR = "../results/00CHECK_SMPL_CONTAMINATION/bam_files"

rule all:
    input:
        expand(os.path.join(SCRATCH_BAM_DIR, "{sublib}_vireo/donor_ids.tsv"), sublib=SUBLIBS.keys())

# Copy BAMs from /nfs to local scratch
rule copy_bam:
    input:  lambda wc: SUBLIBS[wc.sublib]
    output: os.path.join(SCRATCH_BAM_DIR, "{sublib}.bam")
    benchmark: "reports/benchmarks/check_sample_contamination_copy_bam.{sublib}.benchmark.txt"
    log:    "../results/00LOG/00CHECK_SMPL_CONTAMINATION/copy_{sublib}.log"
    shell:  "cp {input} {output} > {log} 2>&1"

rule make_whitelist:
    input:  lambda wildcards: f"/nfs/neurocluster/projects/fetal_brain_single_cell_eQTL_NBray/{wildcards.sublib.split('_')[1]}/{wildcards.sublib}/all-sample/DGE_filtered/cell_metadata.csv"
    output: "../results/00CHECK_SMPL_CONTAMINATION/refs/{sublib}_whitelist.txt"
    shell:
        # Extract first column (bc_wells), skip header
        "cut -d',' -f1 {input} | tail -n +2 > {output}"


rule subset_vcf_chr22:
    input: COMMON_VCF
    output: "../results/00CHECK_SMPL_CONTAMINATION/refs/ALL.TOPMed_freeze5_hg38_dbSNP.chr22.vcf.gz"
    conda: "../envs/cellsnp_lite.yml"
    log: "../results/00LOG/00CHECK_SMPL_CONTAMINATION/subset_vcf_chr22.log"
    shell:
        """
        bcftools view -r chr22 {input} \
        | bcftools view -i 'MAF>=0.1 && INFO/AN>=20' \
        -Oz -o {output} 2> {log}
        bcftools index -t {output}
        """

# Sort full BAMs with sambamba
rule sort_bam:
    input:  os.path.join(SCRATCH_BAM_DIR, "{sublib}.bam")
    output: os.path.join(SCRATCH_BAM_DIR, "{sublib}.sorted.bam")
    conda: "../envs/cellsnp_lite.yml"
    resources: threads=16, mem_mb=80000, time="2-0:00:00"
    priority: 10
    benchmark: "reports/benchmarks/check_sample_contamination_sort_bam.{sublib}.benchmark.txt"
    log:    "../results/00LOG/00CHECK_SMPL_CONTAMINATION/sort_{sublib}.log"
    shell:  "sambamba sort -t {resources.threads} -m 56G -o {output} {input} > {log} 2>&1"

# Index full sorted BAMs with sambamba
rule index_bam:
    input:  os.path.join(SCRATCH_BAM_DIR, "{sublib}.sorted.bam")
    output: os.path.join(SCRATCH_BAM_DIR, "{sublib}.sorted.bam.bai")
    conda: "../envs/cellsnp_lite.yml"
    resources: threads=8, mem_mb=56000, time="2-0:00:00"
    benchmark: "reports/benchmarks/check_sample_contamination_index_bam.{sublib}.benchmark.txt"
    priority: 30
    log:    "../results/00LOG/00CHECK_SMPL_CONTAMINATION/index_{sublib}.log"
    shell:  "sambamba index -t {resources.threads} {input} > {log} 2>&1"


rule extract_and_recode_chr22:
    input:  bam=os.path.join(SCRATCH_BAM_DIR, "{sublib}.sorted.bam"),
            bai=os.path.join(SCRATCH_BAM_DIR, "{sublib}.sorted.bam.bai")
    output: bam=os.path.join(SCRATCH_BAM_DIR, "{sublib}.chr22.filtered.bam"),
            bai=os.path.join(SCRATCH_BAM_DIR, "{sublib}.chr22.filtered.bam.bai")
    conda:  "../envs/cellsnp_lite.yml"
    resources: threads=8, mem_mb=32000, time="2-0:00:00"
    benchmark: "reports/benchmarks/extract_recode_chr22.{sublib}.benchmark.txt"
    log:    "../results/00LOG/00CHECK_SMPL_CONTAMINATION/extract_recode_chr22_{sublib}.log"
    shell:
        """
        # Extract only chr22 (with hg38_22 prefix)
        sambamba view -t {resources.threads} -f bam {input.bam} hg38_22 -o {output.bam} > {log} 2>&1

        # Recode header to remove hg38_ prefix
        samtools view -H {output.bam} | \
        sed 's/^@SQ.*SN:hg38_22/@SQ\tSN:chr22/' > {input}_new_header.sam

        samtools reheader {input}_new_header.sam {output.bam} > {output.bam}.tmp && mv {output.bam}.tmp {output.bam}

        # Index the final BAM
        samtools index {output.bam}

        # Clean up
        rm {input}_new_header.sam
        """



## Extract hg38_22 from sorted BAM files with sambamba
#rule extract_chr:
#    input: 
#        bam=os.path.join(SCRATCH_BAM_DIR, "{sublib}.sorted.bam"),
#        bai=os.path.join(SCRATCH_BAM_DIR, "{sublib}.sorted.bam.bai")
#    output: os.path.join(SCRATCH_BAM_DIR, "{sublib}.hg38_22.sorted.bam")
#    conda: "../envs/cellsnp_lite.yml"
#    resources: threads=8, mem_mb=32000, time="2-0:00:00"
#    priority: 50
#    benchmark: "reports/benchmarks/check_sample_contamination_extract_chr.{sublib}.benchmark.txt"
#    log:    "../results/00LOG/00CHECK_SMPL_CONTAMINATION/extract_chr_{sublib}.log"
#    shell:  "sambamba view -t {resources.threads} -f bam {input.bam} hg38_22 -o {output} > {log} 2>&1"

# Rm hg38_ from bam chr encoding
#rule recode_bam_chr:
#    input: os.path.join(SCRATCH_BAM_DIR, "{sublib}.hg38_22.sorted.bam")
#    output:
#        bam=os.path.join(SCRATCH_BAM_DIR, "{sublib}.chr22.sorted.bam"),
#    conda: "../envs/cellsnp_lite.yml"
#    envmodules: "samtools/1.9", "bcftools/1.16.0"
#    resources: threads=4, mem_mb=16000, time="2-0:00:00"
#    priority: 60    
#    benchmark: "reports/benchmarks/check_sample_contamination_recode_bam_chr.{sublib}.benchmark.txt"
#    log: "../results/00LOG/00CHECK_SMPL_CONTAMINATION/recode_chr_{sublib}.log"
#    shell: """
#           # Generate rename map from BAM header
#           samtools view -H {input} | grep '^@SQ' | awk '{{split($2,a,":"); print a[2] "\t" (a[2]=="chr22" ? "chr22" : "chr"substr(a[2],6))}}' > chrom_recode.txt
#
#           # Apply new header
#           samtools view -H {input} | \
#             sed -f <(awk '{{print "s/"$1"/"$2"/"}}' chrom_recode.txt) \
#             > new_header.sam
#
#             samtools reheader new_header.sam {input} > {output.bam}
#
#             # Index the new BAM
#             samtools index {output.bam}
#
#             # Clean up
#             rm new_header.sam chrom_recode.txt
#             """

# Rule to run cellsnp-lite
rule cellsnp_lite:
    input:  bam=os.path.join(SCRATCH_BAM_DIR, "{sublib}.chr22.filtered.bam"),
            vcf=rules.subset_vcf_chr22.output,
            whitelist=rules.make_whitelist.output
    output: "../results/00CHECK_SMPL_CONTAMINATION/cellsnp_lite/{sublib}_pileup/cellSNP.cells.vcf.gz"
    conda: "../envs/cellsnp_lite.yml"
    resources: threads=8, mem_mb=80000
    priority: 70
    benchmark: "reports/benchmarks/check_sample_contamination_cellsnp_lite.{sublib}.benchmark.txt"
    params: "../results/00CHECK_SMPL_CONTAMINATION/cellsnp_lite/{sublib}_pileup/" 
    log:    "../results/00LOG/00CHECK_SMPL_CONTAMINATION/cellsnp_lite_{sublib}.log"
    shell:
        """
        cellsnp-lite -s {input.bam} \
                     -O {params} \
                     --cellTAG CB --UMItag UB \
                     -R {input.vcf} \
                     --barcodeFile {input.whitelist} \
                     --genotype --gzip \
                     --chrom chr22 2>&1 \
                     --minMAF 0.05 --minCOUNT 1 \
                     -p {threads} 2>&1 | tee -a {log}
        """

# Rule to run Vireo
rule vireo:
    input:  pileup_vcf = "../results/00CHECK_SMPL_CONTAMINATION/cellsnp_lite/{sublib}_pileup/cellSNP.cells.vcf.gz",
            ref_vcf = REF_VCF
    output: os.path.join(SCRATCH_BAM_DIR, "{sublib}_vireo/donor_ids.tsv")
    conda: "../envs/vireo.yml" 
    resources: threads = 2, mem_mb = 8000
    priority: 80
    benchmark: "reports/benchmarks/check_sample_contamination_vireo.{sublib}.benchmark.txt"
    params: os.path.join(SCRATCH_BAM_DIR, "{sublib}_vireo/") 
    log:    "../results/00LOG/00CHECK_SMPL_CONTAMINATION/vireo_{sublib}.log"
    shell:
        """
        vireo -c {input.pileup_vcf} \
              -d {input.ref_vcf} \
              -o {params} 2>&1 | tee -a {log}
        """

# Optional: Clean up extracted BAMs to save space
#rule cleanup:
#    input:  donor_ids="{scratch_dir}/{sublib}_vireo/donor_ids.tsv"
#    output: touch="{scratch_dir}/{sublib}_cleanup.done"
#    shell:
#        """
#        rm {wildcards.scratch_dir}/{wildcards.sublib}_{wildcards.chrom}.bam
#        echo "==== Cleaned up extracted BAM for {wildcards.sublib} ====" | tee -a logs/cleanup_{wildcards.sublib}.log
#        """
