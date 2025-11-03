configfile: '../config/config.yaml'

rule all:
    input:
        expand("../results/12CTWAS/output/ctwas_{cell_type}_{gwas}_pip.tsv", cell_type = config['cell_types'], gwas = config['gwas'])

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
            gwas = "../results/07PREP-GWAS/{gwas}_hg38_ldsr_ready.sumstats.gz",
            weights = "../results/12CTWAS/stubs/{cell_type}_copy_weights.done"
    output: "../results/12CTWAS/output/ctwas_{cell_type}_{gwas}_pip.tsv"
    params: ld_dir = "../results/12CTWAS/ld_mats/",
            weights_dir = "../results/12CTWAS/weights/{cell_type}/" 
    resources: threads = 6, mem_mb = 24000, time="3-0:00:00"
    singularity: config["containers"]["twas"]
    message:  "Creating LD matrices and variant info files for causal TWAS"
    benchmark: "reports/benchmarks/12ctwas.run_ctwas_{cell_type}_{gwas}.benchmark.txt"
    log:    "../results/00LOG/12CTWAS/ctwas_run_{cell_type}_{gwas}.log"
    script: "../scripts/ctwas_run.R"
