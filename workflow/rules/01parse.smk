import json

configfile: "../config/config.yaml"

localrules: cat_fqs

MERGE_FQ = json.load(open(config['MERGE_FQ_JSON']))
ALL_SAMPLES = sorted(MERGE_FQ.keys())

rule get_refs:
    output:  fa = "../resources/refs/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz",
             gtf = "../resources/refs/Homo_sapiens.GRCh38.109.gtf.gz"
    params:  outdir = "../resources/refs/"
    message: "Downloading genome reference files"
    log:     "../results/00LOG/01PARSE/get_refs.log"
    shell:
             r"""
             wget https://ftp.ensembl.org/pub/release-109/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz -P {params.outdir}
             wget https://ftp.ensembl.org/pub/release-109/gtf/homo_sapiens/Homo_sapiens.GRCh38.109.gtf.gz -P {params.outdir}
             """

rule mk_ref:
    input:   fa = "../resources/refs/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz",
             gtf = "../resources/refs/Homo_sapiens.GRCh38.109.gtf.gz"
    output:  "../resources/refs/ref.done" 
    conda:   "../envs/parse.yml"
    params:  outdir = "../resources/refs/hg38", name = "hg38"
    resources: threads = 6, mem_mb = 64000, time="0-1:00:00"
    message: "Creating processed reference sequence index files"
    log:     "../results/00LOG/01parse/mk_ref.log"
    shell:
       	     r"""
             split-pipe \
             --mode mkref \
             --genome_name {params.name} \
             --threads 
             --fasta {input.fa} \
             --genes {input.gtf} \
             --output_dir {params.outdir} 

             touch {params.outdir}/{output}
       	     """

rule cat_fqs:
    input:  r1 = lambda wildcards: MERGE_FQ[wildcards.sample]['R1'],
            r2 = lambda wildcards: MERGE_FQ[wildcards.sample]['R2']
    output: r1 = "../results/01MRGD_fqs/{sample}_R1.fastq.gz",
            r2 = "../results/01MRGD_fqs/{sample}_R2.fastq.gz"
    log:    r1 = "../results/00LOG/01MRGD_fqs/{sample}_R1.log",
            r2 = "../results/00LOG/01MRGD_fqs/{sample}_R2.log"
    params: indir = config['fq_dir'], 
            outdir = "../results/01MRGD_fqs/"
#    priority: 10
    shell:
        """
        cat {input.r1} > {params.outdir}{wildcards.sample}_R1.fastq.gz 2> {log.r1} && \
        cat {input.r2} > {params.outdir}{wildcards.sample}_R2.fastq.gz 2> {log.r2}
        """


rule run_parse:
    input:  r1 = "../results/01MRGD_fqs/{sample}_R1.fastq.gz",
            r2 = "../results/01MRGD_fqs/{sample}_R2.fastq.gz"
    output:  "../results/02PARSE/{sample}/run.done"
    conda:   "../envs/parse.yml"
    params:  refdir = "../resources/refs/hg38"
    resources: threads = 16, mem_mb = 128000, time="1-0:00:00"
    message: "Running Parse combine"
    log:     "../results/00LOG/02parse/run_parse_{sample}.log"
    shell:
             """
             split-pipe \
               --mode all \
               --kit WT_mega \
               --chemistry v1 \
               --genome_dir {params.refdir} \
               --fq1 {input.r1} \
               --fq2 {input.r1} \
               --samp_list ../config/sample-list.txt \
               --output_dir ../results/02PARSE/{wildcards.sample} 2> {log}
              """

rule run_parse_combine:
    input:  expand("../results/02PARSE/{sample}/run.done", sample = ALL_SAMPLES)
    output: "../results/02PARSE/combine/run.done"
    conda: "../envs/Parse.yml"
    resources: threads = 16, mem_mb = 128000, time="1-0:00:00"
    message: "Combining Parse for fastq files"
    log:     "../results/00LOG/02parse/run_parse_combine.log"
    shell:
        """
        split-pipe \
        --mode comb \
        ../results/02PARSE/2 \
        ../results/02PARSE/3 \
        ../results/02PARSE/4 \
        ../results/02PARSE/5 \
        ../results/02PARSE/6 \
        ../results/02PARSE/7 \
        ../results/02PARSE/8 \
        ../results/02PARSE/9 \
        ../results/02PARSE/10 \
        ../results/02PARSE/11 \
        ../results/02PARSE/12 \
        ../results/02PARSE/13 \
        ../results/02PARSE/14 \
        ../results/02PARSE/15 \
        ../results/02PARSE/16 \
        --output_dir ../results/02PARSE/combine 2> {log}
        touch {output}
	"""


