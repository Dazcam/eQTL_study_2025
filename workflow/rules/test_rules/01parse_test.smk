# It's possible to split spipe process into 3 rules to save time and resources
# This would be a far more efficent way to run things, but needs testing 

rule run_parse_pre:
    input:  r1 = "../results/01MRGD_fqs/{sample}_R1.fastq.gz",
            r2 = "../results/01MRGD_fqs/{sample}_R2.fastq.gz",
            ref = "../resources/refs/ref.done"
    output:  "../results/02PARSE/{sample}/pre_run.done"
    conda:   "../envs/spipe.yml"
    params:  refdir = "../resources/refs/hg38", 
             sample_list = "../config/sample-list_plate1.txt"
    priority: 50
    benchmark: "reports/benchmarks/{sample}.run_parse_pre.benchmark.txt"
    resources: threads = 1, mem_mb = 32000, time="1-0:00:00"
    message: "Running Parse combine"
    log:     "../results/00LOG/02PARSE/run_parse_pre_{sample}.log"
    shell:
             """
             source activate spipe-1.3.1
             split-pipe \
               --mode inqc \
               --until_step pre \
               --kit WT_mega \
               --nthreads {threads} \
               --chemistry v2 \
               --genome_dir {params.refdir} \
               --fq1 {input.r1} \
               --fq2 {input.r2} \
               --samp_list {params.sample_list} \
               --output_dir ../results/02PARSE/{wildcards.sample} 2> {log}
             touch {output} 
             """

rule run_parse_align:
    input:  r1 = "../results/01MRGD_fqs/{sample}_R1.fastq.gz",
            r2 = "../results/01MRGD_fqs/{sample}_R2.fastq.gz",
            ref = "../resources/refs/ref.done",
            parse_pre = "../results/02PARSE/{sample}/pre_run.done"
    output:  "../results/02PARSE/{sample}/align_run.done"
    conda:   "../envs/spipe.yml"
    params:  refdir = "../resources/refs/hg38",
             sample_list = "../config/sample-list_plate1.txt"
    priority: 60
    benchmark: "reports/benchmarks/{sample}.run_parse_align.benchmark.txt"
    resources: threads = 32, mem_mb = 256000, time="3-0:00:00"
    message: "Running Parse combine"
    log:     "../results/00LOG/02PARSE/run_parse_align_{sample}.log"
    shell:
             """
             source activate spipe-1.3.1
             split-pipe \
               --mode align \
               --until_step dge	\
               --kit WT_mega \
               --nthreads {threads} \
               --chemistry v2 \
               --genome_dir {params.refdir} \
               --fq1 {input.r1} \
               --fq2 {input.r2} \
               --samp_list {params.sample_list} \
               --output_dir ../results/02PARSE/{wildcards.sample} 2> {log}
             touch {output}
             """

rule run_parse_ana:
    input:  r1 = "../results/01MRGD_fqs/{sample}_R1.fastq.gz",
            r2 = "../results/01MRGD_fqs/{sample}_R2.fastq.gz",
            ref = "../resources/refs/ref.done",
            parse_align = "../results/02PARSE/{sample}/align_run.done"
    output:  "../results/02PARSE/{sample}/ana_run.done"
    conda:   "../envs/spipe.yml"
    params:  refdir = "../resources/refs/hg38",
             sample_list = "../config/sample-list_plate1.txt"
    priority: 70
    benchmark: "reports/benchmarks/{sample}.run_parse_ana.benchmark.txt"
    resources: threads = 32, mem_mb = 256000, time="3-0:00:00"
    message: "Running Parse combine"
    log:     "../results/00LOG/02PARSE/run_parse_ana_{sample}.log"
    shell:
             """
             source activate spipe-1.3.1
             split-pipe \
               --mode ana \
               --kit WT_mega \
               --nthreads {threads} \
               --chemistry v2 \
               --genome_dir {params.refdir} \
               --fq1 {input.r1} \
               --fq2 {input.r2} \
               --samp_list {params.sample_list} \
               --reuse \
               --output_dir ../results/02PARSE/{wildcards.sample} 2> {log}
             touch {output}
             """

rule run_parse_combine:
    input:  expand("../results/02PARSE/{sample}/ana_run.done", sample = ALL_SAMPLES)
    output: "../results/02PARSE/combine_plate1/run.done"
    conda: "../envs/spipe.yml"
    resources: threads = 32, mem_mb = 256000, time="3-0:00:00"
    benchmark: "reports/benchmarks/parse_combine_plate1.benchmark.txt"
    message: "Combining Parse for fastq files"
    log:     "../results/00LOG/02PARSE/run_parse_combine_plate1.log"
    shell:
        """
        source activate spipe-1.3.1
        split-pipe \
        --mode comb \
        --nthreads {threads} \
       --sublibraries ../results/02PARSE/2_plate1 \
       ../results/02PARSE/3_plate1 \
        ../results/02PARSE/4_plate1 \
        ../results/02PARSE/5_plate1 \
        ../results/02PARSE/6_plate1 \
        ../results/02PARSE/7_plate1 \
        ../results/02PARSE/8_plate1 \
        ../results/02PARSE/9_plate1 \
        ../results/02PARSE/10_plate1 \
        ../results/02PARSE/11_plate1 \
        ../results/02PARSE/12_plate1 \
        ../results/02PARSE/13_plate1 \
        ../results/02PARSE/14_plate1 \
        ../results/02PARSE/15_plate1 \
        ../results/02PARSE/16_plate1 \
        --output_dir ../results/02PARSE/combine_plate1 2> {log}
        touch {output}
	"""
