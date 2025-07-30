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

# Cat ref files; just tracking frq as prefix of other only required in smr rule
rule cat_refs:
    # This rule uses the LDSR hg38 1000G frq files, but smr throws an allele frq error
    input:  frq = expand(config["smr"]["cat_refs"]["frq_in"], chr = range(1, 23)),
            bed = expand(config["smr"]["cat_refs"]["bed_in"], chr = range(1, 23)),
            bim = expand(config["smr"]["cat_refs"]["bim_in"], chr = range(1, 23)),
            fam = config["smr"]["cat_refs"]["fam_in"]
    output: frq = config["smr"]["cat_refs"]["cat_frq"],
            bed = config["smr"]["cat_refs"]["cat_bed"],
            bim = config["smr"]["cat_refs"]["cat_bim"],
            fam = config["smr"]["cat_refs"]["cat_fam"]
    log: config["smr"]["cat_refs"]["log"]
    shell: 
        """
        head -n 1 {input.frq[0]} > {output.frq}
        for file in {input.frq}; do
            tail -n +2 $file >> {output.frq}
        done

        cat {input.bed} > {output.bed}
        cat {input.bim} > {output.bim}
        cp {input.fam} {output.fam}

        echo "Concatenated frq, bed, bim, fam files" > {log}
        """

rule get_snp_positions:
    input: rules.cat_refs.output.bim
    output: config["smr"]["get_snp_positions"]["snp_pos"]
    envmodules: "bcftools"
    log: config["smr"]["get_snp_positions"]["log"]
#    shell: "bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\n' {input} > {output}" # if input is vcf
    shell: """awk '{{print $1 "\t" $4 "\t" $2 "\t" $6 "\t" $5}}' {input} > {output}"""

# Move this to tensorflow rules when debugged
rule cat_tensorqtl_nom_snps:
    input: config["output_files"]["tensorqtl_output"]
    output: cat_snps = config["smr"]["cat_tensorqtl_nom_snps"]["cat_snps"],
            summary = config["smr"]["cat_tensorqtl_nom_snps"]["summary"]
    params: config["smr"]["cat_tensorqtl_nom_snps"]["dir"] 
    log: config["smr"]["cat_tensorqtl_nom_snps"]["log"]
    shell:
        """
        python scripts/cat_tensorqtl_nom_snps.py \
          --nom_dir {params} \
          --cell_type {wildcards.cell_type} \
          --concat_out {output.cat_snps} \
          --summary_out {output.summary}  >> {log} 2>&1
        """

rule create_query:
    input:  eqtl = rules.cat_tensorqtl_nom_snps.output.cat_snps,
#            eqtl = config["output_files"]["tensorqtl_perm_output"], # All probes failed HEIDI when using perm
            snps = rules.get_snp_positions.output,
            genes = config["slsdr"]["prep_susie_gene_meta"]["output"],
            frq = rules.cat_refs.output.frq,
            r_script = config["smr"]["create_query"]["r_script"]
    output: query = config["smr"]["create_query"]["query"],
            gene_lst = config["smr"]["create_query"]["gene_lst"]
#    params: qval_thresh = config["eqtl_fdr"]
    singularity: config["containers"]["R"]
    resources: threads = 4, mem_mb = 20000
    log: config["smr"]["create_query"]["log"]
    script: "../scripts/prep_eQTL_for_smr.R"

rule create_besd:
    input: query = config["smr"]["create_query"]["query"],
           smr = rules.get_smr_binary.output
    output: config["smr"]["create_besd"]["besd"]
    params: config["smr"]["create_besd"]["prefix"]
    resources: threads = 4, mem_mb = 20000
    envmodules: "compiler/gnu/5/5.0"
    log: config["smr"]["create_besd"]["log"]
    shell:
        """
        {input.smr} --qfile {input.query} --make-besd --out {params} > {log} 2>&1
        """

# Generate frq files from my samples; fails smr due to allele freq discrepancies
#rule generate_topmed_freq:
#    input:  bed = config["smr"]["generate_topmed_freq"]["bed"],
#            bim = config["smr"]["generate_topmed_freq"]["bim"],
#            fam = config["smr"]["generate_topmed_freq"]["fam"]
#    output: config["smr"]["generate_topmed_freq"]["out"]
#    params: prefix_in = config["smr"]["generate_topmed_freq"]["prefix_in"],
#            prefix_out = config["smr"]["generate_topmed_freq"]["prefix_out"]
#    envmodules: "plink/1.9"  
#    shell:
#        """
#        plink --bfile {params.prefix_in} \
#              --freq \
#              --out {params.prefix_out}
#        """

# Rule to format GWAS data; my script my be causing the error
rule format_gwas:
    input:  gwas = config["smr"]["format_gwas"]["gwas"],
            frq  = rules.cat_refs.output.frq,
            script = "scripts/prep_gwas_for_smr.R"
    output: config["smr"]["format_gwas"]["ma"]
    params: config["smr"]["format_gwas"]["prefix"]
    resources: threads = 4, mem_mb = 20000
    singularity: config["containers"]["R"]
    log: config["smr"]["format_gwas"]["log"]
    script: "../scripts/prep_gwas_for_smr.R"    

rule smr:
    input:  bin = rules.get_smr_binary.output,
            gwas = rules.format_gwas.output,
            besd = rules.create_besd.output
    output: config["smr"]["smr"]["smr"]
    params: geno_prefix = config["smr"]["smr"]["geno_prefix"],
            besd_prefix = config["smr"]["create_besd"]["prefix"],
            out_prefix = config["smr"]["smr"]["smr_prefix"]
    resources: threads = 4, mem_mb = 20000
    envmodules: "compiler/gnu/5/5.0"
    log:    config["smr"]["smr"]["log"]
    shell:  """
            {input.bin} --bfile {params.geno_prefix} \
                --gwas-summary {input.gwas} \
                --beqtl-summary {params.besd_prefix} \
                --out {params.out_prefix} \
                --peqtl-smr 0.01 >> {log} 2>&1 
            """

rule smr_report:
    # Note diff paths for output and out_file; Rmarkdown needs outfile to be relative to Rmd file
    input:  smr_files = expand(rules.smr.output, cell_type = config["cell_types"], gwas = config["gwas"]),
            rmd_script = "scripts/smr_report.Rmd"
    output: config["smr"]["smr_report"]["html"]
    params: cell_types = ','.join(['\'{}\''.format(x) for x in config["cell_types"]]),
            in_dir = config["smr"]["smr_report"]["in_dir"],
            output_file = "../reports/07smr_report.html"
    singularity: config["containers"]["R"]
    log: config["smr"]["smr_report"]["log"]
    shell:
        """
        Rscript -e "rmarkdown::render('{input.rmd_script}', \
            output_file = '{params.output_file}', \
            params = list(cell_types = c({params.cell_types}), in_dir = '{params.in_dir}'))" > {log} 2>&1
        """
