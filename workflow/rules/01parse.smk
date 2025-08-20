import json

configfile: "../config/config.yaml"

localrules: get_refs, cat_fqs

## Note: there is a clash between snakemake env and spipe env (using yaml and dir) that I can't resolve
## As a workaround I just activated spipe env within the shell script

MERGE_FQ = json.load(open(config['MERGE_FQ_JSON']))
ALL_SAMPLES = sorted(MERGE_FQ.keys())

rule get_refs:
    output:  fa = config["parse"]["get_refs"]["fa_out"],
             gtf = config["parse"]["get_refs"]["gtf_out"]
    params:  outdir = config["parse"]["get_refs"]["outdir"],
             fa_link = config["parse"]["get_refs"]["fa_link"],
       	     gtf_link = config["parse"]["get_refs"]["gtf_link"]
    message: "Downloading genome reference files"
    log:     config["parse"]["get_refs"]["log"]
    shell:
             r"""
             wget {params.fa_link} -P {params.outdir}
             wget {params.gtf_link} -P {params.outdir}
             """

rule mk_ref:
    input:   fa = rules.get_refs.output.fa,
             gtf = rules.get_refs.output.gtf
    output:  config["parse"]["mk_ref"]["outfile"]
    priority: 50
    params:  outdir = config["parse"]["mk_ref"]["outdir"],
             build = config["parse"]["mk_ref"]["build"]
    resources: threads = 6, mem_mb = 64000, time="1-0:00:00"
    message: "Creating processed reference sequence index files"
    log:     config["parse"]["mk_ref"]["log"]
    shell:
       	     r"""
             source activate spipe-1.3.1
             split-pipe \
             --mode mkref \
             --genome_name {params.build} \
             --nthreads {threads} \
             --fasta {input.fa} \
             --genes {input.gtf} \
             --output_dir {params.outdir} 2> {log} 

             touch {output}
       	     """

rule cat_fqs:
    input:  r1 = lambda wildcards: MERGE_FQ[wildcards.sample]['R1'],
            r2 = lambda wildcards: MERGE_FQ[wildcards.sample]['R2']
    output: r1 = temp(config["parse"]["cat_fqs"]["out_r1"]),
            r2 = temp(config["parse"]["cat_fqs"]["out_r2"])
    log:    r1 = config["parse"]["cat_fqs"]["log_r1"],
            r2 = config["parse"]["cat_fqs"]["log_r2"]
    benchmark: "reports/benchmarks/{sample}.cat_fq.benchmark.txt"
    params: outdir = config["parse"]["cat_fqs"]["outdir"],
            fq_size = "reports/benchmarks/input_fq_sizes.tsv"
    shell:
        """
        cat {input.r1} > {params.outdir}{wildcards.sample}_R1.fastq.gz 2> {log.r1} && \
        cat {input.r2} > {params.outdir}{wildcards.sample}_R2.fastq.gz 2> {log.r2} 
        """

rule run_parse:
    input:  r1 = rules.cat_fqs.output.r1,
            r2 = rules.cat_fqs.output.r2
    output: config["parse"]["run_parse"]["outfile"]
    params:  refdir = config["parse"]["mk_ref"]["outdir"],
    priority: 50
    benchmark: "reports/benchmarks/{sample}.run_parse.benchmark.txt"
    resources: threads = 32, mem_mb = 360000, time="10-0:00:00"
    message: "Running Parse alignment"
    log:     config["parse"]["run_parse"]["log"]
    shell:
             """
             source activate spipe-1.3.1
             split-pipe \
               --mode all \
               --kit WT_mega \
               --chemistry v2 \
               --genome_dir {params.refdir} \
               --fq1 {input.r1} \
               --fq2 {input.r2} \
               --samp_list ../config/sample-list_plate3.txt \
               --output_dir ../results/01PARSE/{wildcards.sample} 2> {log}
             touch {output}
             """

rule run_parse_combine:
    input:  expand("../results/01PARSE/{sample}/run.done", sample = ALL_SAMPLES)
    output: "../results/01PARSE/combine_plate3/run.done" 
    resources: threads = 32, mem_mb = 360000, time="3-0:00:00"
    benchmark: "reports/benchmarks/run_parse_combine.benchmark.txt"
    message: "Combining Parse for fastq files"
    log:     "../results/00LOG/01PARSE/run_parse_combine.log"
    shell:
        """
        source activate spipe-1.3.1
        split-pipe \
        --mode comb \
        --sublib_list ../config/sublib_lst_plate3.txt \
        --output_dir ../results/01PARSE/combine_plate3 2> {log}
        touch {output}
	"""

# It's possible to split spipe process into 3 rules to save time and resources

#rule run_parse_pre:
#    input:  r1 = "../results/01MRGD_fqs/{sample}_R1.fastq.gz",
#            r2 = "../results/01MRGD_fqs/{sample}_R2.fastq.gz",
#            ref = "../resources/refs/ref.done"
#    output:  "../results/02PARSE/{sample}/pre_run.done"
#    conda:   "../envs/spipe.yml"
#    params:  refdir = "../resources/refs/hg38", 
#             sample_list = "../config/sample-list_plate1.txt"
#    priority: 50
#    benchmark: "reports/benchmarks/{sample}.run_parse_pre.benchmark.txt"
#    resources: threads = 1, mem_mb = 32000, time="1-0:00:00"
#    message: "Running Parse combine"
#    log:     "../results/00LOG/02PARSE/run_parse_pre_{sample}.log"
#    shell:
#             """
#             source activate spipe-1.3.1
#             split-pipe \
#               --mode inqc \
#               --until_step pre \
#               --kit WT_mega \
#               --nthreads {threads} \
#               --chemistry v2 \
#               --genome_dir {params.refdir} \
#               --fq1 {input.r1} \
#               --fq2 {input.r2} \
#               --samp_list {params.sample_list} \
#               --output_dir ../results/02PARSE/{wildcards.sample} 2> {log}
#             touch {output} 
#             """

#rule run_parse_align:
#    input:  r1 = "../results/01MRGD_fqs/{sample}_R1.fastq.gz",
#            r2 = "../results/01MRGD_fqs/{sample}_R2.fastq.gz",
#            ref = "../resources/refs/ref.done",
#            parse_pre = "../results/02PARSE/{sample}/pre_run.done"
##    output:  "../results/02PARSE/{sample}/align_run.done"
#    conda:   "../envs/spipe.yml"
#    params:  refdir = "../resources/refs/hg38",
#             sample_list = "../config/sample-list_plate1.txt"
#    priority: 60
#    benchmark: "reports/benchmarks/{sample}.run_parse_align.benchmark.txt"
#    resources: threads = 32, mem_mb = 256000, time="3-0:00:00"
#    message: "Running Parse combine"
#    log:     "../results/00LOG/02PARSE/run_parse_align_{sample}.log"
#    shell:
#             """
#             source activate spipe-1.3.1
#             split-pipe \
#               --mode align \
#               --until_step dge	\
#               --kit WT_mega \
#               --nthreads {threads} \
#               --chemistry v2 \
#               --genome_dir {params.refdir} \
#               --fq1 {input.r1} \
#               --fq2 {input.r2} \
#               --samp_list {params.sample_list} \
#               --output_dir ../results/02PARSE/{wildcards.sample} 2> {log}
#             touch {output}
#             """

#rule run_parse_ana:
#    input:  r1 = "../results/01MRGD_fqs/{sample}_R1.fastq.gz",
#            r2 = "../results/01MRGD_fqs/{sample}_R2.fastq.gz",
#            ref = "../resources/refs/ref.done",
#            parse_align = "../results/02PARSE/{sample}/align_run.done"
#    output:  "../results/02PARSE/{sample}/ana_run.done"
##    conda:   "../envs/spipe.yml"
#    params:  refdir = "../resources/refs/hg38",
#             sample_list = "../config/sample-list_plate1.txt"
#    priority: 70
#    benchmark: "reports/benchmarks/{sample}.run_parse_ana.benchmark.txt"
#    resources: threads = 32, mem_mb = 256000, time="3-0:00:00"
#    message: "Running Parse combine"
#    log:     "../results/00LOG/02PARSE/run_parse_ana_{sample}.log"
#    shell:
#             """
#             source activate spipe-1.3.1
#             split-pipe \
#               --mode ana \
#               --kit WT_mega \
#               --nthreads {threads} \
#               --chemistry v2 \
#               --genome_dir {params.refdir} \
#               --fq1 {input.r1} \
#               --fq2 {input.r2} \
#               --samp_list {params.sample_list} \
#               --reuse \
#               --output_dir ../results/02PARSE/{wildcards.sample} 2> {log}
#             touch {output}
#             """

#rule run_parse_combine:
#    input:  expand("../results/02PARSE/{sample}/ana_run.done", sample = ALL_SAMPLES)
#    output: "../results/02PARSE/combine_plate1/run.done"
#    conda: "../envs/spipe.yml"
#    resources: threads = 32, mem_mb = 256000, time="3-0:00:00"
#    benchmark: "reports/benchmarks/parse_combine_plate1.benchmark.txt"
#    message: "Combining Parse for fastq files"
#    log:     "../results/00LOG/02PARSE/run_parse_combine_plate1.log"
#    shell:
#        """
#        source activate spipe-1.3.1
#        split-pipe \
#        --mode comb \
#        --nthreads {threads} \
#       --sublibraries ../results/02PARSE/2_plate1 \
#       ../results/02PARSE/3_plate1 \
#        ../results/02PARSE/4_plate1 \
#        ../results/02PARSE/5_plate1 \
#        ../results/02PARSE/6_plate1 \
#        ../results/02PARSE/7_plate1 \
#        ../results/02PARSE/8_plate1 \
#        ../results/02PARSE/9_plate1 \
#        ../results/02PARSE/10_plate1 \
#        ../results/02PARSE/11_plate1 \
#        ../results/02PARSE/12_plate1 \
#        ../results/02PARSE/13_plate1 \
#        ../results/02PARSE/14_plate1 \
#        ../results/02PARSE/15_plate1 \
#        ../results/02PARSE/16_plate1 \
#        --output_dir ../results/02PARSE/combine_plate1 2> {log}
#        touch {output}
#	"""


