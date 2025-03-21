configfile: "../config/config.yaml"

rule zip_cnts:
    input:   config["input_files"]["counts"]
    output:  config["output_files"]["counts_gz"]
    singularity: config["containers"]["fastqtl"] # Use bgzip and tabix in container
    log:     config["log_files"]["fastqtl_zip"]
    shell:
             """
             bgzip -c {input} > {output}
             tabix -p bed {output}
             """

rule fast_qtl:
    input:  counts = config["output_files"]["counts_gz"],
            genotypes = config["input_files"]["genotypes"],
            covariates = config["input_files"]["covariates"]
    output: config["output_files"]["fastqtl_output"]
    singularity: config["containers"]["fastqtl"]
    resources: threads = 1, mem_mb = 6000, time="5:00:00"
    params: chunk = lambda wc: {wc.chunk},
            num_chunks = config["eQTL"]["num_chunks"]
    log:    config["log_files"]["fastqtl"]
    shell:
            """
            fastQTL --vcf {input.genotypes} \
                    --bed {input.counts} \
                    --chunk {params.chunk} {params.num_chunks} \
                    --cov {input.covariates} \
                    --out {output} \
                    --log {log} \
                    --normal >> {log} 2>&1
            """

rule fast_qtl_perms:
    input:  counts = config["output_files"]["counts_gz"],
            genotypes = config["input_files"]["genotypes"],
            covariates = config["input_files"]["covariates"]
    output: config["output_files"]["fastqtl_perm_output"]
    singularity: config["containers"]["fastqtl"]
    resources: threads = 1, mem_mb = 6000, time="5:00:00"
    params: chunk = lambda wc: {wc.chunk},
            num_chunks = config["eQTL"]["num_chunks"],
            min = config["eQTL"]["perm_min"],
            max = config["eQTL"]["perm_max"]
    log:    config["log_files"]["fastqtl_perm"]
    shell:
            """
            fastQTL --vcf {input.genotypes} \
                    --bed {input.counts} \
                    --chunk {params.chunk} {params.num_chunks} \
                    --cov {input.covariates} \
                    --permute {params.min} {params.max} \
                    --out {output} \
                    --log {log} \
                    --normal >> {log} 2>&1
            """

rule convert_genotypes:
    input:  config["input_files"]["genotypes"]
    output: config["output_files"]["genotypes_plink"]
    params: config["eQTL"]["out_prefix"]
    envmodules: "plink/2.0"
    shell: "plink2 --vcf {input} --make-pgen --out {params}" 
     
rule tensorqtl:
    input:  genotypes = config["output_files"]["genotypes_plink"],
            counts = config["output_files"]["counts_gz"],
            covariates = config["input_files"]["covariates"]
    output: config["output_files"]["tensorqtl_output"]
    singularity: config["containers"]["tensorqtl"]
    resources: threads = 10, mem_mb = 100000, time="5:00:00"
    params: plink_prefix = config["eQTL"]["out_prefix"],
            out_prefix = config["output_files"]["tensorqtl_out_prefix"],
            window = config["eQTL"]["window"]
    log:    config["log_files"]["tensorqtl"]
    shell:
            """
            python3 -m tensorqtl {params.plink_prefix} {input.counts} {params.out_prefix} \
               --covariates {input.covariates} \
               --window {params.window} \
               --mode cis_nominal >> {log} 2>&1
             """

rule tensorqtl_perm:
    input:  genotypes = config["output_files"]["genotypes_plink"],
            counts = config["output_files"]["counts_gz"],
            covariates = config["input_files"]["covariates"]
    output: out = config["output_files"]["tensorqtl_perm_output"],
            log = config["log_files"]["tensorqtl_perm"]
    singularity: config["containers"]["tensorqtl"]
    resources: threads = 10, mem_mb = 100000, time="5:00:00"
    params: plink_prefix = config["eQTL"]["out_prefix"],
            out_prefix = config["output_files"]["tensorqtl_perm_out_prefix"],
            window = config["eQTL"]["window"]
    log:    config["log_files"]["tensorqtl_perm"]
    shell:
            """
            python3 -m tensorqtl {params.plink_prefix} {input.counts} {params.out_prefix} \
               --covariates {input.covariates} \
               --window {params.window} \
               --mode cis >> {log} 2>&1 
            """

rule tensorqtl_cat_log:
    input:  expand(config["log_files"]["tensorqtl_perm"], cell_type = config['cell_types'])
    output: config["output_files"]["tensorqtl_cat_log_output"],
    shell:  """python scripts/cat_tensorqtl_logs.py -i "{input}" -o {output}"""

rule plot_qtl:
    input:  genotypes = config["input_files"]["genotypes"],
            pairs_file = config["plot_qtl"]["pairs_file"]
    output: config["plot_qtl"]["out_file"]
    singularity: config["containers"]["tensorqtl"]
    params: expression_dir = config["plot_qtl"]["expression_dir"],
            output_dir = config["plot_qtl"]["out_dir"]
    log:    config["plot_qtl"]["log"]
    shell:
            """
            python3 scripts/eqtl_plot.py \
               --pairs_file {input.pairs_file} \
               --genotype_file {input.genotypes} \
               --expression_dir {params.expression_dir} \
               --output_dir {params.output_dir} >> {log}
             """




