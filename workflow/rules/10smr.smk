# Get SNP position from genotype files
rule get_smr_binary:
   output: config["smr"]["get_smr_binary"]["smr"]
   params: config["smr"]["get_smr_binary"]["prefix"]
   shell:
      """
      https://yanglab.westlake.edu.cn/software/smr/download/smr-1.4.0-linux-x86_64.zip
      unzip smr-1.4.0-linux-x86_64.zip 'smr-1.4.0-linux-x86_64/smr' 
      mv smr-1.4.0-linux-x86_64/smr {params}
      rm -rf smr-1.4.0-linux-x86_64.zip smr-1.4.0-linux-x86_64
      """

rule get_snp_positions:
    input: config["input_files"]["genotypes"]
    output: config["smr"]["get_snp_positions"]["snp_pos"]
    envmodules: "bcftools"
    log: config["smr"]["get_snp_positions"]["log"]
    shell: "bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\n' {input} > {output}"    

rule create_query:
    input:  eqtl = config["output_files"]["tensorqtl_perm_output"],
            snps = rules.get_snp_positions.output,
            gene = config["slsdr"]["prep_susie_gene_meta"]["output"],
            r_script = config["smr"]["create_query"]["r_script"]
    output: query = config["smr"]["create_query"]["query"],
            gene_lst = config["smr"]["create_query"]["gene_lst"]
    params: qval_thresh = config["eqtl_fdr"],
    singularity: config["containers"]["R"]
    log: config["smr"]["query"]["log"]
    script: "{input.r_script}"

rule create_besd:
    input: query = lambda wildcards: rules.create_query.output.query.format(cell_type=wildcards.cell_type),
           smr = rules.get_smr_binary.output
    output: config["smr"]["create_besd"]["besd"]
    log: config["smr"]["create_besd"]["log"]
    shell:
        """
        {input.smr} --qfile {input.query} --make-besd --out {output} > {log} 2>&1
        """



