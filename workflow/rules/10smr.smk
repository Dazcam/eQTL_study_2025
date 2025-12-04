configfile: "../config/config.yaml"

rule all:
    input:
        "reports/10SMR/10smr_report.html"
#        expand("../results/10SMR/smr/{cell_type}/{cell_type}_{gwas}.smr", cell_type=config["cell_types"], gwas=config["gwas"])
     

rule get_smr_binary:
   output: config["smr"]["get_smr_binary"]["smr"]
   params: config["smr"]["get_smr_binary"]["prefix"]
   message: "Download the smr binary file"
   benchmark: "reports/benchmarks/10smr.get_smr_binary.txt"
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
    message: "Cat hg38 refs into a single file"
    benchmark: "reports/benchmarks/10smr.cat_refs.txt"
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

# Move this to tensorflow rules when debugged
rule cat_tensorqtl_nom_snps:
    input: lambda w, norm_method=config['tensorQTL']['norm_methods'][0],
                geno_pc=config['tensorQTL']['geno_pcs'],
                exp_pc=config['tensorQTL']['exp_pcs'][0]:
                f"../results/05TENSORQTL/tensorqtl_nom/{w.cell_type}_{norm_method}_genPC_{geno_pc}_expPC_{exp_pc}/{w.cell_type}_{norm_method}_nom.cis_qtl_pairs.4.parquet"
    output: cat_snps = config["smr"]["cat_tensorqtl_nom_snps"]["cat_snps"],
            summary = config["smr"]["cat_tensorqtl_nom_snps"]["summary"]
    params: config["smr"]["cat_tensorqtl_nom_snps"]["dir"]
    message: "Cat TensorQTL parquet files into single file for SMR"
    benchmark: "reports/benchmarks/10smr.cat_tensorqtl_nom_snps_{cell_type}.txt"
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
            snps = rules.cat_refs.output.bim,
            genes = config["susie"]["prep_susie_gene_meta"]["output"],
            frq = rules.cat_refs.output.frq,
            r_script = config["smr"]["create_query"]["r_script"]
    output: query = config["smr"]["create_query"]["query"],
            gene_lst = config["smr"]["create_query"]["gene_lst"]
#    params: qval_thresh = config["eqtl_fdr"]
    singularity: config["containers"]["R"]
    resources: threads = 4, mem_mb = 20000
    message: "Create query file for SMR"
    benchmark: "reports/benchmarks/10smr.create_query_{cell_type}.txt"
    log: config["smr"]["create_query"]["log"]
    script: "../scripts/smr_create_query.R"

rule create_besd:
    input: query = config["smr"]["create_query"]["query"],
           smr = rules.get_smr_binary.output
    output: config["smr"]["create_besd"]["besd"]
    params: config["smr"]["create_besd"]["prefix"]
    resources: threads = 4, mem_mb = 20000
    envmodules: "compiler/gnu/5/5.0"
    message: "Create besd file for SMR"
    benchmark: "reports/benchmarks/10smr.create_besd_{cell_type}.txt"
    log: config["smr"]["create_besd"]["log"]
    shell:
        """
        {input.smr} --qfile {input.query} --add-n 133 --make-besd --out {params} > {log} 2>&1
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
            script = "scripts/smr_format_gwas.R"
    output: config["smr"]["format_gwas"]["ma"]
    params: config["smr"]["format_gwas"]["prefix"]
    resources: threads = 4, mem_mb = 20000
    singularity: config["containers"]["R"]
    message: "Prep GWAS into .ma format for SMR"
    benchmark: "reports/benchmarks/10smr.format_gwas_{gwas}.txt"
    log: config["smr"]["format_gwas"]["log"]
    script: "../scripts/smr_format_gwas.R"    

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
    message: "Run SMR"
    benchmark: "reports/benchmarks/10smr.smr_{cell_type}_{gwas}.txt"
    log:    config["smr"]["smr"]["log"]
    shell:  """
            {input.bin} --bfile {params.geno_prefix} \
                --gwas-summary {input.gwas} \
                --beqtl-summary {params.besd_prefix} \
                --out {params.out_prefix} \
                --peqtl-smr 5e-8 >> {log} 2>&1 
            """

# These rules can be used to generate code for smr plots
# Note that smr_plot doesn't produce output for snp / gene pairs
# If using these rules need to get code to handel this
#rule prep_gene_lists:
#    input: rules.create_query.output.gene_lst
#    output: ensembl = config['smr']['prep_gene_lists']['ensembl'],
#            symbol = config['smr']['prep_gene_lists']['symbol']
#    singularity: config["containers"]["R"]
#    message: "Creating the gene lists file for SMR plotting"
#    benchmark: "reports/benchmarks/10smr.prep_gene_lists_{cell_type}.txt"
#    log: config["smr"]["prep_gene_lists"]["log"]
#    script: "../scripts/prep_gene_lists_for_smr_plts.R"

#rule smr_plot:
#    input:  gwas = rules.format_gwas.output,
#            gene_list = rules.prep_gene_lists.output.ensembl,
#            bin = rules.get_smr_binary.output,
#    output: config["smr"]["smr_plot"]["output"]
#    params: geno_prefix = config["smr"]["smr"]["geno_prefix"],
#            besd_prefix = config["smr"]["create_besd"]["prefix"],
#            out_prefix = config["smr"]["smr_plot"]["out_prefix"],
#            smr_gene = config['smr_gene'],
#            window =  config['smr_window']
#    resources: threads = 4, mem_mb = 20000
#    envmodules: "compiler/gnu/5/5.0"
#    message: "Generate data for SMR plotting"
#    benchmark: "reports/benchmarks/10smr.smr_plot_{cell_type}_{gwas}_{smr_gene}.txt"
#    log:    config["smr"]["smr_plot"]["log"]
#    shell:  """
#            {input.bin} --bfile {params.geno_prefix} \
#              --gwas-summary {input.gwas} \
#              --beqtl-summary {params.besd_prefix} \
#              --out {params.out_prefix} \
#              --plot \
#              --probe {params.smr_gene} \
#              --probe-wind {params.window} \
#              --gene-list {input.gene_list} >> {log} 2>&1
#            """ 

rule smr_report:
    # Note diff paths for output and out_file; Rmarkdown needs outfile to be relative to Rmd file
    input:  smr_files = expand(rules.smr.output, cell_type = config["cell_types"], gwas = config["gwas"]),
#            plt_files = expand(rules.smr_plot.output, cell_type = config["cell_types"], gwas = config["gwas"], smr_gene = config['smr_gene']),
            rmd_script = "scripts/smr_report.Rmd"
    output: config["smr"]["smr_report"]["html"]
    params: cell_types = ','.join(['\'{}\''.format(x) for x in config["cell_types"]]),
            in_dir = config["smr"]["smr_report"]["in_dir"],
            output_file = "../reports/10SMR/10smr_report.html",
            p_smr = config["p_smr"],
            p_heidi =  config["p_heidi"]    
    resources: threads = 4, mem_mb = 20000
    singularity: config["containers"]["r_eqtl"]
    message: "Generate SMR report"
    benchmark: "reports/benchmarks/10smr.smr_report.txt"
    log: config["smr"]["smr_report"]["log"]
    shell:
        """
        Rscript -e "rmarkdown::render('{input.rmd_script}', \
            output_file = '{params.output_file}', \
            params = list(cell_types = c({params.cell_types}), in_dir = '{params.in_dir}', p_smr = '{params.p_smr}', p_heidi = '{params.p_heidi}', bmark_dir = '../reports/benchmarks/'))" > {log} 2>&1
        """

