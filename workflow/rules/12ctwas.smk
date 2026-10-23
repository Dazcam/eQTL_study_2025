configfile: '../config/config.yaml'

localrules: ctwas_report

rule create_ld_matrices:
    output: config["ctwas"]["create_ld_matrices"]["output"]
    params: ld_ref_dir = config["ctwas"]["create_ld_matrices"]["ld_ref_dir"],
            out_dir    = config["ctwas"]["create_ld_matrices"]["out_dir"]
    resources: threads = 1, mem_mb = 10000, time="5:30:00"
    singularity: config["containers"]["twas"]
    message:  "Creating LD matrices and variant info files for causal TWAS"
    benchmark: config["ctwas"]["create_ld_matrices"]["bmrk"]
    log:       config["ctwas"]["create_ld_matrices"]["log"]
    script:    config["ctwas"]["create_ld_matrices"]["script"]


rule copy_fusion_weights:
    """
    Copy the real .wgt.RDat files (column 2 of the .pos file) into a clean
    directory that will be used by ctwas.  This removes the stub/skipped
    files that make preprocess_weights() choke.
    Note: Keep stub files in separate directory from weights or cTWAS chokes 
    """
    input:   config["ctwas"]["copy_fusion_weights"]["input"]
    output:  touch(config["ctwas"]["copy_fusion_weights"]["output"])
    params:  src_dir   = config["ctwas"]["copy_fusion_weights"]["src_dir"],
             dest_dir  = config["ctwas"]["copy_fusion_weights"]["dest_dir"]
    message:  "Copying genuine FUSION weights for {wildcards.cell_type} (listed in .pos) to {params.dest_dir}"
    log:       config["ctwas"]["copy_fusion_weights"]["log"]
    benchmark: config["ctwas"]["copy_fusion_weights"]["bmrk"]
    shell:  """
            mkdir -p {params.dest_dir} && \\
            awk 'NR>1 {{print $2}}' {input} | \\
            while read wgt; do \\
            cp -v {params.src_dir}/"$wgt" {params.dest_dir}/ ; \\
            done > {log} 2>&1
            """

rule run_ctwas:
    input:  ld_mat   = rules.create_ld_matrices.output,
            gwas     = config["ctwas"]["run_ctwas"]["gwas"],
            weights  = rules.copy_fusion_weights.output,
            bim_file = config["ctwas"]["run_ctwas"]["bim_file"]
    output: config["ctwas"]["run_ctwas"]["output"]
    params: ld_dir        = config["ctwas"]["create_ld_matrices"]["out_dir"],
            weights_dir   = config["ctwas"]["run_ctwas"]["weights_dir"],
            run_weight_qc = config["ctwas"]["run_ctwas"]["run_weight_qc"]
    resources: threads = 6, mem_mb = 96000, time="1:00:00" # mem_mb can vary 
    threads: 6
    singularity: config["containers"]["twas"]
    message:     "Running cTWAS"
    benchmark:   config["ctwas"]["run_ctwas"]["bmrk"]
    log:         config["ctwas"]["run_ctwas"]["log"]
    script:      config["ctwas"]["run_ctwas"]["script"]

rule ctwas_report:
    # Note diff paths for output and out_file; Rmarkdown needs outfile to be relative to Rmd file
    input:
            ctwas_res  = expand(rules.run_ctwas.output, cell_type = config['cell_types'], gwas = config['gwas']),
            rmd_script = config["ctwas"]["ctwas_report"]["rmd_script"]
    output: config["ctwas"]["ctwas_report"]["output"]
    params: in_dir      = config["ctwas"]["ctwas_report"]["in_dir"],
            bmark_dir   = config["ctwas"]["ctwas_report"]["bmark_dir"],
            lookup_dir  = config["ctwas"]["ctwas_report"]["lookup_dir"],
            cell_types  = ",".join(config["cell_types"]),
            gwas        = ",".join(config["gwas"].keys()),
            output_file = config["ctwas"]["ctwas_report"]["output_file"]
    singularity: config["containers"]["r_eqtl"]
    message:     "Generate cTWAS report"
    benchmark:   config["ctwas"]["ctwas_report"]["bmrk"]
    log:         config["ctwas"]["ctwas_report"]["log"]
    shell:
        """
        Rscript -e "rmarkdown::render('{input.rmd_script}', \
            output_file = '{params.output_file}', \
            params = list(in_dir = '{params.in_dir}', \
                          bmark_dir = '{params.bmark_dir}', \
                          lookup_dir = '{params.lookup_dir}', \
                          cell_types = '{params.cell_types}', \
                          gwas = '{params.gwas}'))" > {log} 2>&1
        """


