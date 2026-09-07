configfile: "../config/config.yaml"

localrules: create_parse_json, setup_liftover

PLATES = list(config["fastq"].keys())

rule all:
    input:
#        config["setup"]["get_containers"]["output"],
#        expand(config["setup"]["create_parse_json"]["output"], plate=PLATES),
#        config["setup"]["liftover"]["output"],
#        config["setup"]["get_fugita_data"]["output"],
#        config["setup"]["get_jang_singlebrain"]["output"]
        config["setup"]["get_obrien_supp"]["output"]

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

rule setup_liftover:
    output:  config["setup"]["liftover"]["output"]
    params:  outdir = config["setup"]["liftover"]["outdir"]
    log:     config["setup"]["liftover"]["log"]
    shell:
             """
             wget -nc -P {params.outdir} http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/liftOver
             wget -nc -P {params.outdir} http://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz
             wget -nc -P {params.outdir} http://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz
        
             chmod +x {params.outdir}/liftOver
             touch {output}
             """

# Requires scanpy env so need to set up envs first
rule get_fugita_data:
    output:  config["setup"]["get_fugita_data"]["output"]
    params:  outdir=config["setup"]["get_fugita_data"]["outdir"],
             token=config["setup"]["get_fugita_data"]["token"]
    conda:   config["scanpy"]["env"]
    message: "Downloading Fugita et al. cell cell type and subtype data from Synapse"
    log:     config["setup"]["get_fugita_data"]["log"]
    shell:
             """
             python scripts/setup_get_fugita_data.py \
             --token {params.token} \
             --outdir {params.outdir} > {log} 2>&1
             touch {output}
             """

rule get_jang_singlebrain:
    output:  config["setup"]["get_jang_singlebrain"]["output"]
    params:  outdir=config["setup"]["get_jang_singlebrain"]["outdir"]
    conda:   config["scanpy"]["env"]
    message: "Downloading Jang et al. SingleBrain eQTL summary statistics from Zenodo"
    log:     config["setup"]["get_jang_singlebrain"]["log"]
    shell:
             """
             mkdir -p {params.outdir}
             zenodo_get 10.5281/zenodo.14908182 -o {params.outdir} > {log} 2>&1
             touch {output}
             """

#rule get_obrien_data:
#    output: all_qtl = config["setup"]["get_obrien_data"]["all_qtl"],
#            top_qtl = config["setup"]["get_obrien_data"]["top_qtl"]
#    params: out_dir = config["setup"]["get_obrien_data"]["out_dir"],
#            web_link = config["setup"]["get_obrien_data"]["web_link"],
#            zip_out = config["setup"]["get_obrien_data"]["zip_out"]
#    message: "Download bulk brain gene eQTL files from O'Brien 2018, PMID:30419947"
#    log:    config["setup"]["get_obrien_data"]["log"]
#    shell:  """
#            curl -L -b "cookies.txt" -A "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/91.0.4472.124 Safari/537.36" -o {params.zip_out} {params.web_link} &&
#            unzip -j -o {params.zip_out} all_eqtls_gene.txt.gz -d {params.out_dir} &&
#            unzip -j -o {params.zip_out} top_eqtls_gene.txt.gz -d {params.out_dir} &&
#            rm -f {params.zip_out} 2>> {log}
#            """

rule get_obrien_supp:
    output: config["setup"]["get_obrien_supp"]["output"],
    params: web_link = config["setup"]["get_obrien_supp"]["web_link"],
    message: "Download bulk eQTL supp tables for from O'Brien 2018, PMID:30419947"
    log:    config["setup"]["get_obrien_supp"]["log"]
    shell:   "wget -O {output} {params.web_link}"  
