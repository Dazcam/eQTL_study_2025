configfile: "../config/config.yaml"

rule prep_tensorQTL_input:
    input:  cov_file = config["tensorQTL"]["prep_tensorQTL_input"]["cov_file"],
            sex_file = config["tensorQTL"]["prep_tensorQTL_input"]["sex_file"],
            gene_lookup = config["tensorQTL"]["prep_tensorQTL_input"]["gene_lookup"] 
    output: cov_out = config["tensorQTL"]["prep_tensorQTL_input"]["cov_out"],
            exp_out = config["tensorQTL"]["prep_tensorQTL_input"]["exp_out"]
    params: pseudoblk_dir = config["tensorQTL"]["prep_tensorQTL_input"]["pseudoblk_dir"],
            report_dir = config["tensorQTL"]["prep_tensorQTL_input"]["report_dir"],
            out_dir = config["tensorQTL"]["prep_tensorQTL_input"]["out_dir"],
            batch_var = config["tensorQTL"]["prep_tensorQTL_input"]["batch_var"],
            norm_method = "{norm_method}"
    singularity: config["containers"]["r_eqtl"]
    resources: threads = 1, mem_mb = 6000, time="5:00:00"
    message: "Prep pseudoblk GeX and covariate matricies for tensorQTL with norm: {wildcards.norm_method}"
    benchmark: "reports/benchmarks/05tensorQTL.prep_tensorQTL_input_{cell_type}_{norm_method}.txt"
    log:    config["tensorQTL"]["prep_tensorQTL_input"]["log"]
    script: "../scripts/prep_tensorQTL_input_files.R"

rule zip_pblk_cnts:
    input:   rules.prep_tensorQTL_input.output.exp_out
    output:  config["tensorQTL"]["zip_pblk_cnts"]["output"]
    singularity: config["containers"]["fastqtl"] # Use bgzip and tabix in container
    message: "bgzip and index pseudoblk counts for norm: {wildcards.norm_method}"
    benchmark: "reports/benchmarks/05tensorQTL.zip_pbulk_cnts_{cell_type}_{norm_method}.txt"
    log:     config["tensorQTL"]["zip_pblk_cnts"]["log"]
    shell:
             """
             bgzip -c {input} > {output}
             tabix -p bed {output}
             """

rule convert_genotypes:
    input:  config["tensorQTL"]["convert_genotypes"]["input"]
    output: config["tensorQTL"]["convert_genotypes"]["output"]
    params: config["tensorQTL"]["convert_genotypes"]["prefix_out"]
    envmodules: "plink/2.0"
    message: "Convert genotypes to plink2 format for tensorQTL"
    benchmark: "reports/benchmarks/05tensorQTL.convert_genotypes.txt"
    shell: "plink2 --vcf {input} --make-pgen --out {params}" 

rule split_covariates:
    input:  rules.prep_tensorQTL_input.output.cov_out
    output:  config["tensorQTL"]["split_covariates"]["output"],
    singularity: config["containers"]["R"]
    resources: threads = 1, mem_mb = 6000, time="5:00:00"
    message: "Divide covariate file to test different PC thresholds in tensorQTL for norm: {wildcards.norm_method}"
    benchmark: "reports/benchmarks/05tensorQTL.split_covariates_{cell_type}_{norm_method}_genPC_{geno_pc}_expPC_{exp_pc}.txt"   
    log:    config["tensorQTL"]["split_covariates"]["log"]
    script: "../scripts/split_covariates_for_tensorQTL.R"

rule tensorqtl_nom:
    input:  genotypes = rules.convert_genotypes.output,
            counts = rules.zip_pblk_cnts.output,
            covariates = rules.split_covariates.output
    output: config["tensorQTL"]["tensorqtl_nom"]["output"]
    params: prefix_in = config["tensorQTL"]["tensorqtl_nom"]["prefix_in"],
            prefix_out = config["tensorQTL"]["tensorqtl_nom"]["prefix_out"],
            window = config["tensorQTL"]["window"]
    singularity: config["containers"]["tensorqtl"]
    resources: threads = 10, mem_mb = 100000, time="5:00:00"
    message: "Run tensorQTL nominal for norm: {wildcards.norm_method}"
    benchmark: "reports/benchmarks/05tensorQTL.nom_{cell_type}_{norm_method}_genPC_{geno_pc}_expPC_{exp_pc}.txt"
    log:    config["tensorQTL"]["tensorqtl_nom"]["log"]
    shell:
            """
            python3 -m tensorqtl {params.prefix_in} {input.counts} {params.prefix_out} \
               --covariates {input.covariates} \
               --window {params.window} \
               --mode cis_nominal >> {log} 2>&1
            """

rule tensorqtl_perm:
    input:  genotypes = rules.convert_genotypes.output,
            counts = rules.zip_pblk_cnts.output,
            covariates = rules.split_covariates.output
    output: config["tensorQTL"]["tensorqtl_perm"]["output"]
    params: prefix_in = config["tensorQTL"]["tensorqtl_perm"]["prefix_in"],
            prefix_out = config["tensorQTL"]["tensorqtl_perm"]["prefix_out"],
            window = config["tensorQTL"]["window"]
    singularity: config["containers"]["tensorqtl"]
    resources: threads = 10, mem_mb = 100000, time="5:00:00"
    message: "Run tensorQTL permutation for norm: {wildcards.norm_method}"
    benchmark: "reports/benchmarks/05tensorQTL.perm_{cell_type}_{norm_method}_genPC_{geno_pc}_expPC_{exp_pc}.txt"
    log:    config["tensorQTL"]["tensorqtl_perm"]["log"]
    shell:
            """
            python3 -m tensorqtl {params.prefix_in} {input.counts} {params.prefix_out} \
               --covariates {input.covariates} \
               --window {params.window} \
               --mode cis >> {log} 2>&1 
            """

rule tensotqtl_report:
    # Note diff paths for output and out_file; Rmarkdown needs outfile to be relative to Rmd file
    input:  qtl = expand(rules.tensorqtl_perm.output, cell_type=config["cell_types"],geno_pc=config["tensorQTL"]["geno_pcs"],exp_pc=config["tensorQTL"]["exp_pcs"],norm_method=config["tensorQTL"]["norm_methods"]),
            rmd_script = "scripts/tensorQTL_report.Rmd"
    output: "reports/05TENSORQTL/tensotqtl_report.html"
    params: in_dir = "../../results/05TENSORQTL/tensorqtl_perm/",
            bmark_dir = "../reports/benchmarks/",
            output_file = "../reports/05TENSORQTL/tensotqtl_report.html",
    singularity: config["containers"]["r_eqtl"]
    message: "Generate tensorQTL report"
    benchmark: "reports/benchmarks/05tensorQTL.tensorqtl_report.benchmark.txt"
    log:     "../results/00LOG/05TENSORQTL/tensorqtl_report.log"
    shell:
        """
        Rscript -e "rmarkdown::render('{input.rmd_script}', \
            output_file = '{params.output_file}', \
            params = list(in_dir = '{params.in_dir}', bmark_dir = '{params.bmark_dir}'))" > {log} 2>&1
        """

#rule tensorqtl_tss_and_sumstats:
#    input:  expand(config["output_files"]["tensorqtl_perm_log"], cell_type = config['cell_types'])
#    output: config["output_files"]["tensorqtl_tss_plt"],
#            config["output_files"]["tensorqtl_tss_tbl"]
#    singularity: config["containers"]["R"]
#    params: root_dir = config["root_dir"],
#            cell_types = config["cell_types"]
#    log:    config["log_files"]["tensorqtl_tss"]    
#    script: "scripts/tensorqtl_tss_and_sumstats.R"

#rule plot_qtl:
#    input:  genotypes = config["input_files"]["genotypes"],
#            pairs_file = config["plot_qtl"]["pairs_file"]
#    output: config["plot_qtl"]["out_file"]
#    singularity: config["containers"]["tensorqtl"]
#    params: expression_dir = config["plot_qtl"]["expression_dir"],
#            output_dir = config["plot_qtl"]["out_dir"]
#    log:    config["plot_qtl"]["log"]
#    shell:
#            """
#            python3 scripts/eqtl_plot.py \
#               --pairs_file {input.pairs_file} \
#               --genotype_file {input.genotypes} \
#               --expression_dir {params.expression_dir} \
#               --output_dir {params.output_dir} >> {log}
#             """
