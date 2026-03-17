configfile: "../config/config.yaml"

rule get_containers:
    output:  config["setup"]["get_containers"]["output"],
    message: "Downloading singulairty / apptainer containers"
    log:     config["setup"]["get_containers"]["log"]
    shell:
             r"""
             scripts/setup_pull_quayio_containers.sh
             """
  
