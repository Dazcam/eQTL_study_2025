configfile: "../config/config.yaml" 

rule all:
    input: 
#        config["review"]["singler_report"]["html"],
#        config["review"]["prep_singler_query"]["h5ad"],
#        config["review"]["prep_singler_query"]["umap_csv"],
#        config["review"]["singler_classify"]["csv"]
        config["review"]["batch_correction_lisi"]["csv"],
        config["review"]["batch_correction_lisi"]["plot"],
#        config["review"]["batch_correction_kbet"]["csv"],
#        config["review"]["batch_correction_kbet"]["plot"]

# Move to set up / or ref download rule
rule dwnld_qian: 
    output: config["review"]["dwnld_qian"]["output"]
    params: web_link = config["review"]["dwnld_qian"]["web_link"],
    message: "Download snRNAseq spatial data from Qian 2025, PMID:40369074"
    benchmark: "reports/benchmarks/15review.dwnld_qian.txt"
    log:    config["review"]["dwnld_qian"]["log"]
    shell:  """
            wget -O {output} {params.web_link} &>> {log}
            """

rule prep_qian_ref:
    input:  config["review"]["dwnld_qian"]["output"]
    output: config["review"]["prep_qian_ref"]["output"]
    conda:  config["scanpy"]["env"]
    resources: threads = 4, mem_mb = 64000, time = "2:00:00"
    params:    n_cells_per_cluster = 500,
               random_seed         = 1234
    log:    config["review"]["prep_qian_ref"]["log"]
    message: "Preparing Qian 2025 MERFISH reference: filtering V1/V2, subsampling H2 clusters"
    script: config["review"]["prep_qian_ref"]["script"]

rule prep_singler_query:
    input:   query_h5ad = config["review"]["prep_singler_query"]["query_h5ad"],
             ref_h5ad   = config["review"]["prep_qian_ref"]["output"]
    output:  h5ad = config["review"]["prep_singler_query"]["h5ad"],
             umap_csv = config["review"]["prep_singler_query"]["umap_csv"]
    conda:   config["scanpy"]["env"]
    resources: threads = 10, mem_mb = 200000, time = "2:00:00"
    log:     config["review"]["prep_singler_query"]["log"]
    message: "Preparing snRNA-seq query: assigning L1 labels, subsetting to MERFISH panel"
    script:  config["review"]["prep_singler_query"]["script"]

rule singler_classify:
    input:   ref_h5ad   = config["review"]["prep_qian_ref"]["output"],
             query_h5ad = config["review"]["prep_singler_query"]["h5ad"]
    output:  config["review"]["singler_classify"]["csv"]
    singularity: config["containers"]["singler"]
    resources: threads  = 8, mem_mb = 128000, time = "2:00:00"
    threads: 8
    log:     config["review"]["singler_classify"]["log"]
    message: "Running SingleR classification: Qian H2 reference vs snRNA-seq query"
    script:  config["review"]["singler_classify"]["script"]

rule batch_correction_lisi:
    input:   h5ad = config["review"]["prep_singler_query"]["query_h5ad"]
    output:  csv  = config["review"]["batch_correction_lisi"]["csv"],
             plot = config["review"]["batch_correction_lisi"]["plot"]
    conda:   config["scanpy"]["env"]
    resources: threads = 10, mem_mb = 200000, time = "2:00:00"
    threads: 10
    params:  batch_col     = "sample",
             cell_type_col = "cell_type",
             random_seed   = config["review"]["batch_correction_lisi"].get("random_seed", 1234),
             perplexity    = config["review"]["batch_correction_lisi"].get("perplexity", 30)
    log:     config["review"]["batch_correction_lisi"]["log"]
    message: "Running LISI (iLISI / cLISI) batch correction assessment pre- and post-Harmony"
    script:  config["review"]["batch_correction_lisi"]["script"]

rule batch_correction_kbet:
    input:   h5ad = config["review"]["prep_singler_query"]["query_h5ad"]
    output:  csv  = config["review"]["batch_correction_kbet"]["csv"],
             plot = config["review"]["batch_correction_kbet"]["plot"]
    conda:   config["scanpy"]["env"]
    resources: threads = 10, mem_mb = 200000, time = "2:00:00"
    threads: 10
    params:  batch_col      = "sample",
             cell_type_col  = "cell_type",
             n_subsample    = config["review"]["batch_correction_kbet"].get("n_subsample", 10000),
             random_seed    = config["review"]["batch_correction_kbet"].get("random_seed", 1234),
             k0             = config["review"]["batch_correction_kbet"].get("k0", 30)
    log:     config["review"]["batch_correction_kbet"]["log"]
    message: "Running kBET batch correction assessment (pre- and post-Harmony) on subsampled cells"
    script:  config["review"]["batch_correction_kbet"]["script"]

rule singler_report:
    input:   singler_csv = config["review"]["singler_classify"]["csv"],
             umap_csv    = config["review"]["prep_singler_query"]["umap_csv"],
             rmd_script  = "scripts/review_singler_report.Rmd"
    output:  config["review"]["singler_report"]["html"]
    params:  singler_csv = "../../results/15REVIEW/singler_results.csv", # Paths relative to the Rmd script location (scripts/)
             umap_csv    = "../../results/15REVIEW/umap_coords.csv",
             output_file = "../reports/15REVIEW/15singler_report.html"
    resources:   threads = 4, mem_mb = 32000, time = "1:00:00"
    singularity: config["containers"]["singler"]
    log:     config["review"]["singler_report"]["log"]
    message: "Generating SingleR spatial validation report"
    shell:
        """
        Rscript -e "rmarkdown::render('{input.rmd_script}', \
            output_file = '{params.output_file}', \
            params = list(
                singler_csv = '{params.singler_csv}',
                umap_csv    = '{params.umap_csv}'))" > {log} 2>&1
        """
