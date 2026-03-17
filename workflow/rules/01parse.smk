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
             source activate spipe
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
             source activate spipe
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
        source activate spipe
        split-pipe \
        --mode comb \
        --sublib_list ../config/sublib_lst_plate3.txt \
        --output_dir ../results/01PARSE/combine_plate3 2> {log}
        touch {output}
	"""

