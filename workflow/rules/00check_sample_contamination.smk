# Define input BAM paths and sublibraries
BAM_PATHS = {
    "13_plate1": "/nfs/neurocluster/projects/fetal_brain_single_cell_eQTL_NBray/plate1/13_plate1/process/barcode_headAligned_anno.bam",
    "5_plate1": "/nfs/neurocluster/projects/fetal_brain_single_cell_eQTL_NBray/plate1/5_plate1/process/barcode_headAligned_anno.bam",
    "9_plate1": "/nfs/neurocluster/projects/fetal_brain_single_cell_eQTL_NBray/plate1/9_plate1/process/barcode_headAligned_anno.bam",
    "15_plate1": "/nfs/neurocluster/projects/fetal_brain_single_cell_eQTL_NBray/plate1/15_plate1/process/barcode_headAligned_anno.bam",
    "14_plate1": "/nfs/neurocluster/projects/fetal_brain_single_cell_eQTL_NBray/plate1/14_plate1/process/barcode_headAligned_anno.bam",
    "3_plate1": "/nfs/neurocluster/projects/fetal_brain_single_cell_eQTL_NBray/plate1/3_plate1/process/barcode_headAligned_anno.bam",
    "7_plate1": "/nfs/neurocluster/projects/fetal_brain_single_cell_eQTL_NBray/plate1/7_plate1/process/barcode_headAligned_anno.bam",
    "2_plate1": "/nfs/neurocluster/projects/fetal_brain_single_cell_eQTL_NBray/plate1/2_plate1/process/barcode_headAligned_anno.bam",
    "4_plate1": "/nfs/neurocluster/projects/fetal_brain_single_cell_eQTL_NBray/plate1/4_plate1/process/barcode_headAligned_anno.bam",
    "10_plate1": "/nfs/neurocluster/projects/fetal_brain_single_cell_eQTL_NBray/plate1/10_plate1/process/barcode_headAligned_anno.bam",
    "11_plate1": "/nfs/neurocluster/projects/fetal_brain_single_cell_eQTL_NBray/plate1/11_plate1/process/barcode_headAligned_anno.bam",
    "12_plate1": "/nfs/neurocluster/projects/fetal_brain_single_cell_eQTL_NBray/plate1/12_plate1/process/barcode_headAligned_anno.bam",
    "16_plate1": "/nfs/neurocluster/projects/fetal_brain_single_cell_eQTL_NBray/plate1/16_plate1/process/barcode_headAligned_anno.bam",
    "6_plate1": "/nfs/neurocluster/projects/fetal_brain_single_cell_eQTL_NBray/plate1/6_plate1/process/barcode_headAligned_anno.bam",
    "8_plate1": "/nfs/neurocluster/projects/fetal_brain_single_cell_eQTL_NBray/plate1/8_plate1/process/barcode_headAligned_anno.bam",
    "14_plate2": "/nfs/neurocluster/projects/fetal_brain_single_cell_eQTL_NBray/plate2/14/process/barcode_headAligned_anno.bam",
    "16_plate2": "/nfs/neurocluster/projects/fetal_brain_single_cell_eQTL_NBray/plate2/16/process/barcode_headAligned_anno.bam",
    "1_plate2": "/nfs/neurocluster/projects/fetal_brain_single_cell_eQTL_NBray/plate2/1/process/barcode_headAligned_anno.bam",
    "10_plate2": "/nfs/neurocluster/projects/fetal_brain_single_cell_eQTL_NBray/plate2/10/process/barcode_headAligned_anno.bam",
    "11_plate2": "/nfs/neurocluster/projects/fetal_brain_single_cell_eQTL_NBray/plate2/11/process/barcode_headAligned_anno.bam",
    "12_plate2": "/nfs/neurocluster/projects/fetal_brain_single_cell_eQTL_NBray/plate2/12/process/barcode_headAligned_anno.bam",
    "13_plate2": "/nfs/neurocluster/projects/fetal_brain_single_cell_eQTL_NBray/plate2/13/process/barcode_headAligned_anno.bam",
    "15_plate2": "/nfs/neurocluster/projects/fetal_brain_single_cell_eQTL_NBray/plate2/15/process/barcode_headAligned_anno.bam",
    "2_plate2": "/nfs/neurocluster/projects/fetal_brain_single_cell_eQTL_NBray/plate2/2/process/barcode_headAligned_anno.bam",
    "4_plate2": "/nfs/neurocluster/projects/fetal_brain_single_cell_eQTL_NBray/plate2/4/process/barcode_headAligned_anno.bam",
    "5_plate2": "/nfs/neurocluster/projects/fetal_brain_single_cell_eQTL_NBray/plate2/5/process/barcode_headAligned_anno.bam",
    "6_plate2": "/nfs/neurocluster/projects/fetal_brain_single_cell_eQTL_NBray/plate2/6/process/barcode_headAligned_anno.bam",
    "7_plate2": "/nfs/neurocluster/projects/fetal_brain_single_cell_eQTL_NBray/plate2/7/process/barcode_headAligned_anno.bam",
    "8_plate2": "/nfs/neurocluster/projects/fetal_brain_single_cell_eQTL_NBray/plate2/8/process/barcode_headAligned_anno.bam",
    "9_plate2": "/nfs/neurocluster/projects/fetal_brain_single_cell_eQTL_NBray/plate2/9/process/barcode_headAligned_anno.bam",
    "5_plate3": "../results/01PARSE/5_plate3/process/barcode_headAligned_anno.bam",
    "15_plate3": "../results/01PARSE/15_plate3/process/barcode_headAligned_anno.bam",
    "4_plate3": "../results/01PARSE/4_plate3/process/barcode_headAligned_anno.bam",
    "14_plate3": "../results/01PARSE/14_plate3/process/barcode_headAligned_anno.bam",
    "10_plate3": "../results/01PARSE/10_plate3/process/barcode_headAligned_anno.bam",
    "3_plate3": "../results/01PARSE/3_plate3/process/barcode_headAligned_anno.bam",
    "12_plate3": "../results/01PARSE/12_plate3/process/barcode_headAligned_anno.bam",
    "11_plate3": "../results/01PARSE/11_plate3/process/barcode_headAligned_anno.bam",
    "1_plate3": "../results/01PARSE/1_plate3/process/barcode_headAligned_anno.bam",
    "2_plate3": "../results/01PARSE/2_plate3/process/barcode_headAligned_anno.bam",
    "8_plate3": "../results/01PARSE/8_plate3/process/barcode_headAligned_anno.bam",
    "6_plate3": "../results/01PARSE/6_plate3/process/barcode_headAligned_anno.bam"
}

SUBLIBS = list(BAM_PATHS.keys())
CHROM = "chr22"  # Adjust to "22" if BAMs lack "chr" prefix
SCRATCH_DIR = "/path/to/scratch/area"  # Replace with your scratch path
COMMON_VCF = "/path/to/common_snps_chr22.vcf.gz"  # Filtered 1000 Genomes SNPs
REF_VCF = "/path/to/ref_chr22.vcf.gz"  # Your reference genotypes for chr22

# Define local rules for head node execution
localrules: index_bam, extract_chr

# Rule to generate all final outputs
rule all:
    input:
        expand("{scratch_dir}/{sublib}_vireo/donor_ids.tsv", scratch_dir=SCRATCH_DIR, sublib=SUBLIBS)

# Rule to index full BAM files (on head node)
rule index_bam:
    input:  bam=lambda wildcards: BAM_PATHS[wildcards.sublib]
    output: bai=lambda wildcards: BAM_PATHS[wildcards.sublib] + ".bai"
    log:    "logs/index_bam_{sublib}.log"
    shell:
        """
        echo "==== Indexing {input.bam} ====" | tee -a {log}
        samtools index {input.bam} 2>&1 | tee -a {log}
        echo "==== Completed indexing {input.bam} ====" | tee -a {log}
        echo | tee -a {log}
        """

# Rule to extract chr22 and copy to scratch (on head node)
rule extract_chr:
    input: bam=lambda wildcards: BAM_PATHS[wildcards.sublib],
           bai=lambda wildcards: BAM_PATHS[wildcards.sublib] + ".bai"
    output: extracted_bam="{scratch_dir}/{sublib}_{chrom}.bam"
    log:  "logs/extract_chr_{sublib}.log"
    shell:
        """
        echo "==== Extracting {wildcards.chrom} from {input.bam} ====" | tee -a {log}
        samtools view -b {input.bam} {wildcards.chrom} > {output.extracted_bam} 2>&1 | tee -a {log}
        echo "==== Completed extraction for {wildcards.sublib} ====" | tee -a {log}
        echo | tee -a {log}
        """

# Rule to run cellsnp-lite (on compute nodes, no index needed)
rule cellsnp_lite:
    input:  bam="{scratch_dir}/{sublib}_{chrom}.bam",
            vcf=COMMON_VCF
    output: vcf="{scratch_dir}/{sublib}_pileup/cellSNP.vcf.gz"
    log:    "logs/cellsnp_lite_{sublib}.log"
    threads: 4
    resources: mem_mb=16000
    shell:
        """
        echo "==== Starting pileup for {wildcards.sublib} ====" | tee -a {log}
        cellsnp-lite -s {input.bam} \
                     -O {wildcards.scratch_dir}/{wildcards.sublib}_pileup \
                     --cellTAG CB --UMItag UB \
                     -R {input.vcf} \
                     --genotype --gzip \
                     --chrom {wildcards.chrom} 2>&1 | tee -a {log}
        echo "==== Completed pileup for {wildcards.sublib} ====" | tee -a {log}
        echo | tee -a {log}
        """

# Rule to run Vireo (on compute nodes)
rule vireo:
    input:  pileup_vcf="{scratch_dir}/{sublib}_pileup/cellSNP.vcf.gz",
            ref_vcf=REF_VCF
    output: donor_ids="{scratch_dir}/{sublib}_vireo/donor_ids.tsv"
    log:    "logs/vireo_{sublib}.log"
    threads: 2
    resources:mem_mb=8000 
    shell:
        """
        echo "==== Starting Vireo for {wildcards.sublib} ====" | tee -a {log}
        vireo -c {input.pileup_vcf} \
              -d {input.ref_vcf} \
              -o {wildcards.scratch_dir}/{wildcards.sublib}_vireo 2>&1 | tee -a {log}
        echo "==== Completed Vireo for {wildcards.sublib} ====" | tee -a {log}
        echo | tee -a {log}
        """

# Optional: Clean up extracted BAMs to save space
rule cleanup:
    input:  donor_ids="{scratch_dir}/{sublib}_vireo/donor_ids.tsv"
    output: touch="{scratch_dir}/{sublib}_cleanup.done"
    shell:
        """
        rm {wildcards.scratch_dir}/{wildcards.sublib}_{wildcards.chrom}.bam
        echo "==== Cleaned up extracted BAM for {wildcards.sublib} ====" | tee -a logs/cleanup_{wildcards.sublib}.log
        """
