localrules: get_sig_eGenes

rule get_sig_eGenes:
    input:  config["output_files"]["tensorqtl_perm_output"]
    output: config["slsdr"]["get_sig_eGenes"]["output"]
    singularity: config["containers"]["R"]
    log:    config["slsdr"]["get_sig_eGenes"]["log"]
    script: "../scripts/get_sig_egenes_with_cis_window.R"

rule prep_susie_input:
    input:  sig_eGenes = config["slsdr"]["get_sig_eGenes"]["output"],
            pseudobulk = config["input_files"]["counts"],
            genotypes = config["output_files"]["genotypes_plink"],
            covariates = config["input_files"]["covariates"] 
    output: config["slsdr"]["prep_susie_input"]["output"]
    envmodules: "plink/2.0"
    resources: threads = 8, mem_mb = 64000, time="1:00:00"
    message: "Generating input files for SuSIE for all sig eGenes"
    log:    config["slsdr"]["prep_susie_input"]["log"]
    script: "../scripts/prep_susie_input_files.py"
