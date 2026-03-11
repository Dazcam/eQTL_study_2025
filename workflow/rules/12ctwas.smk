xconfigfile: '../config/config.yaml'

localrules: ctwas_report

rule create_ld_matrices:
    output: "../results/12CTWAS/ld_mats/create_ctwas_ld_matrices.done"
    params: ld_ref_dir = "../resources/ldsr/ldsr_hg38_refs/plink_files/",
            out_dir = "../results/12CTWAS/ld_mats/"
    resources: threads = 5, mem_mb = 20000, time="3-0:00:00"
    singularity: config["containers"]["twas"]
    message:  "Creating LD matrices and variant info files for causal TWAS"
    benchmark: "reports/benchmarks/12ctwas.create_ctwas_ld_matrices.benchmark.txt"    
    log:    "../results/00LOG/12CTWAS/create_ld_matrices.log"
    script: "../scripts/ctwas_create_ld_matrices.R"


rule copy_fusion_weights:
    """
    Copy the real .wgt.RDat files (column 2 of the .pos file) into a clean
    directory that will be used by ctwas.  This removes the stub/skipped
    files that make preprocess_weights() choke.
    Note: Keep stub files in separate directory from weights or cTWAS chokes 
    """
    input:   "../results/11TWAS/weights/{cell_type}/{cell_type}.pos"
    output:  touch("../results/12CTWAS/weights/stubs/{cell_type}_copy_weights.done")
    params:  src_dir = "../results/11TWAS/weights/{cell_type}",
             dest_dir  = "../results/12CTWAS/weights/{cell_type}"
    message: "Copying genuine FUSION weights for {wildcards.cell_type} (listed in .pos) → {params.dest_dir}"
    log:     "../results/00LOG/12CTWAS/copy_weights_{cell_type}.log"
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
            gwas = "../results/07PREP-GWAS/{gwas}_hg38.tsv",
            weights = "../results/12CTWAS/weights/stubs/{cell_type}_copy_weights.done",
            bim_file = "../results/10SMR/smr_input/1000G.EUR.hg38.bim"
    output: "../results/12CTWAS/output/ctwas_{cell_type}_{gwas}_ctwas.rds"
    params: ld_dir = "../results/12CTWAS/ld_mats/",
            weights_dir = "../results/12CTWAS/weights/{cell_type}/" # cel-specific
    resources: threads = 16, mem_mb = 380000, time="1-0:00:00"
#    resources: threads = 1, mem_mb = 40000, time="1-0:00:00"
#    resources: threads = 6, mem_mb = 96000, time="1-0:00:00" 
#    resources: threads = 6, mem_mb = 148000, time="1-0:00:00"  # InN-0 / adhd
    singularity: config["containers"]["twas"]
    message:  "Running cTWAS"
    benchmark: "reports/benchmarks/12ctwas.run_ctwas_{cell_type}_{gwas}.benchmark.txt"
    log:    "../results/00LOG/12CTWAS/ctwas_run_{cell_type}_{gwas}.log"
    script: "../scripts/ctwas_run.R"

rule run_ctwas_multi:
    input:  ld_mat = rules.create_ld_matrices.output,
            gwas = "../results/07PREP-GWAS/{gwas}_hg38.tsv",
            weights = expand("../results/12CTWAS/weights/stubs/{cell_type}_copy_weights.done", cell_type = config['cell_types']),
            bim_file = "../results/10SMR/smr_input/1000G.EUR.hg38.bim"
    output: "../results/12CTWAS/multi/ctwas_multi_{gwas}_ctwas.rds"
    params: ld_dir = "../results/12CTWAS/ld_mats/",
            weights_dir = "../results/12CTWAS/weights/", # general
            cell_types = config["cell_types"]
#    resources: threads = 16, mem_mb = 380000, time="3-0:00:00"
    resources: threads = 6, mem_mb = 96000, time="1-0:00:00"
    singularity: config["containers"]["twas"]
    message:  "Running multi-group cTWAS"
    benchmark: "reports/benchmarks/12ctwas.run_ctwas_multi_{gwas}.benchmark.txt"
    log:    "../results/00LOG/12CTWAS/ctwas_run_multi_{gwas}.log"
    script: "../scripts/ctwas_run_multi.R"

rule ctwas_report:
    # Note diff paths for output and out_file; Rmarkdown needs outfile to be relative to Rmd file
    input:  
#            ctwas_res = expand(rules.run_ctwas.output, cell_type = config['cell_types'], gwas = config['gwas']),
            ctwas_multi = expand(rules.run_ctwas_multi.output, gwas = config['gwas']),
            rmd_script = "scripts/ctwas_report.Rmd"
    output: "reports/12CTWAS/12ctwas_report.html"
    params: in_dir = "../../results/12CTWAS/",
            bmark_dir = "../reports/benchmarks/",
            lookup_dir = "../../resources/sheets/",
            tbl_dir = "../../results/13MANUSCRIPT_PLOTS_TABLES/tables/",
            output_file = "../reports/12CTWAS/12ctwas_report.html"
    singularity: config["containers"]["r_eqtl"] # Need to add ctwas to r_eqtl conatiner to print locus plot
    message: "Generate cTWAS report"
    benchmark: "reports/benchmarks/12ctwas.ctwas_report.benchmark.txt"
    log:     "../results/00LOG/12CTWAS/ctwas_report.log"
    shell:
        """
        Rscript -e "rmarkdown::render('{input.rmd_script}', \
            output_file = '{params.output_file}', \
            params = list(in_dir = '{params.in_dir}', \
            bmark_dir = '{params.bmark_dir}', \
            lookup_dir = '{params.lookup_dir}', \
            tbl_dir = '{params.tbl_dir}'))" > {log} 2>&1
        """

