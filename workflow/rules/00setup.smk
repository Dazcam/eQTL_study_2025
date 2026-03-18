configfile: "../config/config.yaml"

localrules: create_parse_json

rule get_containers:
    output:  config["setup"]["get_containers"]["output"],
    message: "Downloading singulairty / apptainer containers"
    log:     config["setup"]["get_containers"]["log"]
    shell:
             r"""
             scripts/setup_pull_quayio_containers.sh
             """

rule create_parse_json:
    output:  config["setup"]["create_parse_json"]["output"],
    params:  dirs=lambda wildcards: config["fastq"][wildcards.plate]
    message:  "Create json to map fastq files to sample IDs for Parse alignment"
    log:     config["setup"]["create_parse_json"]["log"]
    shell:
             """
             python scripts/setup_create_parse_json.py \
             --plate {wildcards.plate} \
             --fastq_dirs {params.dirs} > {log} 2>&1
             """


