configfile: '../config/config.yaml'

localrules: ctwas_report

rule create_ld_matrices:
    output: config["ctwas"]["create_ld_matrices"]["output"]
    params: ld_ref_dir = config["ctwas"]["create_ld_matrices"]["ld_ref_dir"]
            out_dir = config["ctwas"]["create_ld_matrices"]["outdir"]
    resources: threads = 5, mem_mb = 20000, time="3-0:00:00"
    singularity: config["containers"]["twas"]
    message:  "Creating LD matrices and variant info files for causal TWAS"
    benchmark: "reports/benchmarks/12ctwas.create_ctwas_ld_matrices.benchmark.txt"    
    log:    config["ctwas"]["create_ld_matrices"]["log"]
    script: "../scripts/ctwas_create_ld_matrices.R"


rule copy_fusion_weights:
    """
    Copy the real .wgt.RDat files (column 2 of the .pos file) into a clean
    directory that will be used by ctwas.  This removes the stub/skipped
    files that make preprocess_weights() choke.
    Note: Keep stub files in separate directory from weights or cTWAS chokes 
    """
    input:   config["ctwas"]["copy_fusion_weights"]["input"]
    output:  touch(config["ctwas"]["copy_fusion_weights"]["output"])
    params:  src_dir = config["ctwas"]["copy_fusion_weights"]["scr_dir"],
             dest_dir  = config["ctwas"]["copy_fusion_weights"]["dest_dir"]
    message: "Copying genuine FUSION weights for {wildcards.cell_type} (listed in .pos) to {params.dest_dir}"
    log:     config["ctwas"]["copy_fusion_weights"]["log"]
    benchmark: "reports/benchmarks/12ctwas.copy_weights_{cell_type}.benchmark.txt"
    shell:  """
            mkdir -p {params.dest_dir} && \\
            awk 'NR>1 {{print $2}}' {input} | \\
            while read wgt; do \\
            cp -v {params.src_dir}/"$wgt" {params.dest_dir}/ ; \\
            done > {log} 2>&1
            """

rule run_ctwas:
    input:  ld_mat = rules.create_ld_matrices.output,
            gwas = config["ctwas"]["run_ctwas"]["gwas"]
            weights = config["ctwas"]["run_ctwas"]["weights"]
            bim_file = config["ctwas"]["run_ctwas"]["bim_file"]
    output: config["ctwas"]["run_ctwas"]["output"]
    params: ld_dir = config["ctwas"]["run_ctwas"]["ld_dir"]
            weights_dir = config["ctwas"]["run_ctwas"]["weights_dir"]
    resources: threads = 16, mem_mb = 380000, time="1-0:00:00"
#    resources: threads = 1, mem_mb = 40000, time="1-0:00:00"
#    resources: threads = 6, mem_mb = 96000, time="1-0:00:00" 
#    resources: threads = 6, mem_mb = 148000, time="1-0:00:00"  # InN-0 / adhd
    singularity: config["containers"]["twas"]
    message:  "Running cTWAS"
    benchmark: "reports/benchmarks/12ctwas.run_ctwas_{cell_type}_{gwas}.benchmark.txt"
    log:    config["ctwas"]["run_ctwas"]["log"]
    script: "../scripts/ctwas_run.R"

rule run_ctwas_multi:
    input:  ld_mat = rules.create_ld_matrices.output,
            gwas = config["ctwas"]["run_ctwas"]["gwas"]
            weights = config["ctwas"]["run_ctwas"]["weights"]
            bim_file = config["ctwas"]["run_ctwas"]["bim_file"]
    output: config["ctwas"]["run_ctwas_multi"]["output"]
    params: ld_dir = config["ctwas"]["run_ctwas"]["ld_dir"],
            weights_dir = config["ctwas"]["run_ctwas"]["weights_dir"],
            cell_types = config["cell_types"]
#    resources: threads = 16, mem_mb = 380000, time="3-0:00:00"
    resources: threads = 6, mem_mb = 96000, time="1-0:00:00"
    singularity: config["containers"]["twas"]
    message:  "Running multi-group cTWAS"
    benchmark: "reports/benchmarks/12ctwas.run_ctwas_multi_{gwas}.benchmark.txt"
    log:    config["ctwas"]["run_ctwas_multi"]["log"]
    script: "../scripts/ctwas_run_multi.R"

rule ctwas_report:
    # Note diff paths for output and out_file; Rmarkdown needs outfile to be relative to Rmd file
    input:  
#            ctwas_res = expand(rules.run_ctwas.output, cell_type = config['cell_types'], gwas = config['gwas']),
            ctwas_multi = expand(rules.run_ctwas_multi.output, gwas = config['gwas']),
            rmd_script = config["ctwas"]["ctwas_report"]["rmd_script"]
    output: config["ctwas"]["ctwas_report"]["output"]
    params: in_dir = config["ctwas"]["ctwas_report"]["in_dir"],
            bmark_dir = config["ctwas"]["ctwas_report"]["bmark_dir"],
            lookup_dir = config["ctwas"]["ctwas_report"]["lookup_dir"],
            tbl_dir = config["ctwas"]["ctwas_report"]["tbl_dir"],
            output_file = config["ctwas"]["ctwas_report"]["output_file"]
    singularity: config["containers"]["r_eqtl"] # Need to add ctwas to r_eqtl conatiner to print locus plot
    message: "Generate cTWAS report"
    benchmark: "reports/benchmarks/12ctwas.ctwas_report.benchmark.txt"
    log:     config["ctwas"]["ctwas_report"]["log"]
    shell:
        """
        Rscript -e "rmarkdown::render('{input.rmd_script}', \
            output_file = '{params.output_file}', \
            params = list(in_dir = '{params.in_dir}', \
            bmark_dir = '{params.bmark_dir}', \
            lookup_dir = '{params.lookup_dir}', \
            tbl_dir = '{params.tbl_dir}'))" > {log} 2>&1
        """

