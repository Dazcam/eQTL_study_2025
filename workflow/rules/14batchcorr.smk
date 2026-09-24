configfile: "../config/config.yaml" 

rule all:
    input: 
        config["review"]["batch_corr_report"]["html"],

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

rule batch_corr_report:
    input:   kbet_csv  = config["review"]["batch_correction_kbet"]["csv"],
             lisi_csv  = config["review"]["batch_correction_lisi"]["csv"].replace('.csv', '_summary.csv'),
             kbet_plot = config["review"]["batch_correction_kbet"]["plot"],
             lisi_plot = config["review"]["batch_correction_lisi"]["plot"]
    output:  html = config["review"]["batch_corr_report"]["html"]
    conda:   config["scanpy"]["env"]
    resources:  mem_mb  = 8000,
                time    = "0:30:00"
    params:     n_cells      = 679738,
                n_batches    = 136,
                n_cell_types = 7,
                harmony_vars = "sample"
    log:     config["review"]["batch_corr_report"]["log"]
    message: "Generating batch correction QC report"
    script:  config["review"]["batch_corr_report"]["script"]
