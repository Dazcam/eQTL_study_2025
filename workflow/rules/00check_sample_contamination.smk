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
        expand("../results/00CHECK_SMPL_CONTAMINATION/vireo/{sublib}_vireo/donor_ids.tsv", sublib=SUBLIBS.keys())

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
    shell:  "samtools sort -@ {resources.threads} -m 4G -o {output} {input} > {log} 2>&1"

# Index full sorted BAMs with sambamba
rule index_bam:
    input:  os.path.join(SCRATCH_BAM_DIR, "{sublib}.sorted.bam")
    output: os.path.join(SCRATCH_BAM_DIR, "{sublib}.sorted.bam.bai")
    conda: "../envs/cellsnp_lite.yml"
    resources: threads=8, mem_mb=56000, time="2-0:00:00"
    benchmark: "reports/benchmarks/check_sample_contamination_index_bam.{sublib}.benchmark.txt"
    priority: 30
    log:    "../results/00LOG/00CHECK_SMPL_CONTAMINATION/index_{sublib}.log"
    shell:  "samtools index -@ {resources.threads} {input} > {log} 2>&1"


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
         samtools view -@ {resources.threads} -b {input.bam} hg38_22 -o {output.bam} > {log} 2>&1

        # Recode header to remove hg38_ prefix
        samtools view -H {output.bam} | \
        sed 's/^@SQ.*SN:hg38_22/@SQ\tSN:chr22/' > {input.bam}_new_header.sam

        samtools reheader {input.bam}_new_header.sam {output.bam} > {output.bam}.tmp && mv {output.bam}.tmp {output.bam}

        # Index the final BAM
        samtools index {output.bam}
        
        # Check for presence of chr22 reads
        READ_COUNT=$(samtools view -c {output.bam} chr22)
        if [ "$READ_COUNT" -eq 0 ]; then
            echo "ERROR: No reads found on chr22 in {output.bam}" >&2
            exit 1
        fi

        # Print basic stats for logging/debugging
        samtools idxstats {output.bam}

        # Clean up
        rm {input.bam}_new_header.sam
        """

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
                     --UMItag None \
                     --cellTAG CB \
                     -R {input.vcf} \
                     --barcodeFile {input.whitelist} \
                     --genotype \
                     --gzip \
                     --chrom chr22 2>&1 \
                     -p {threads} 2>&1 | tee -a {log}
        """

# Rule to run Vireo
rule vireo:
    input:  pileup_vcf = "../results/00CHECK_SMPL_CONTAMINATION/cellsnp_lite/{sublib}_pileup/cellSNP.cells.vcf.gz",
            geno_vcf = REF_VCF
    output: "../results/00CHECK_SMPL_CONTAMINATION/vireo/{sublib}_vireo/donor_ids.tsv"
    conda: "../envs/vireo.yml" 
    resources: threads = 16, mem_mb = 160000
    priority: 80
    benchmark: "reports/benchmarks/check_sample_contamination_vireo.{sublib}.benchmark.txt"
    params: outdir = "../results/00CHECK_SMPL_CONTAMINATION/vireo/{sublib}_vireo",
            tmp = "../results/00CHECK_SMPL_CONTAMINATION/vireo/tmp_{sublib}" 
    log:    "../results/00LOG/00CHECK_SMPL_CONTAMINATION/vireo_{sublib}.log"
    shell:
        """
        # Genotypes VCF checks 
        zcat {input.geno_vcf} | grep -v '^#' | cut -f 1 | sort | uniq >> {log}   # Check chrs
        zcat {input.geno_vcf} | grep -v '^##' | wc -l                 >> {log}   # Count SNPs
        zcat {input.geno_vcf} | grep '^#CHROM' | cut -f 10- | wc -w   >> {log}   # Count samples

        # Check overlaps
        zcat ${input.pileup_vcf} | grep -v '^##' | cut -f 1-2 > {params.tmp}_cellsnp_positions.txt
        zcat ${input.geno_vcf} | grep -v '^##' | cut -f 1-2 > {params.tmp}_donor_positions.txt
        wc -l {params.tmp}_cellsnp_positions.txt >> {log}
        wc -l {params.tmp}_donor_positions.txt >> {log}
        comm -12 <(sort {params.tmp}_cellsnp_positions.txt) <(sort {params.tmp}_donor_positions.txt) | wc -l >> {log}
        rm {params.tmp}*         

        vireo -c {input.pileup_vcf} \
              -d {input.geno_vcf} \
              -t GT \
              -o {params.outdir} 2>&1 | tee -a {log}
        """

