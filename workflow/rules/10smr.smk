# Get SNP position from genotype files

rule get_smr_binary:
   output: config["smr"]["get_smr_binary"]["smr"]
   params: config["smr"]["get_smr_binary"]["prefix"]
   shell:
      """
      wget https://yanglab.westlake.edu.cn/software/smr/download/smr-1.4.0-linux-x86_64.zip
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
            genes = config["slsdr"]["prep_susie_gene_meta"]["output"],
            r_script = config["smr"]["create_query"]["r_script"]
    output: query = config["smr"]["create_query"]["query"],
            gene_lst = config["smr"]["create_query"]["gene_lst"]
    params: qval_thresh = config["eqtl_fdr"],
    singularity: config["containers"]["R"]
    log: config["smr"]["create_query"]["log"]
    script: "../scripts/prep_eQTL_for_smr.R"

rule create_besd:
    input: query = config["smr"]["create_query"]["query"],
           smr = rules.get_smr_binary.output
    output: config["smr"]["create_besd"]["besd"]
    params: config["smr"]["create_besd"]["prefix"]
    envmodules: "compiler/gnu/5/5.0"
    log: config["smr"]["create_besd"]["log"]
    shell:
        """
        {input.smr} --qfile {input.query} --make-besd --out {params} > {log} 2>&1
        """

CHR = list(range(1, 23))

# Rule to concatenate frequency files
rule cat_freq:
    input:  expand(["smr"]["cat_freq"]["frq"], chr = CHR)
    output: ["smr"]["cat_freq"]["cat_frq"]
    shell:
        """
        head -n 1 {input.frq_files[0]} > {output}
        for file in {input.frq_files}; do
            tail -n +2 $file >> {output}
        done
        """

# Rule to format GWAS data
rule format_gwas:
    input:  gwas = config["smr"]["format_gwas"]["gwas"],
            frq  = results.cat_freq.output,
            script = "scripts/format_gwas_for_smr.R"
    output: config["smr"]["format_gwas"]["ma"]
    params: config["smr"]["format_gwas"]["prefix"]
    singularity: config["containers"]["R"]
    log: config["smr"]["format_gwas"]["log"]
    shell:
        """
        Rscript {input.script} {input.gwas} {input.ref_af} {output.ma} > {log} 2>&1
        """


