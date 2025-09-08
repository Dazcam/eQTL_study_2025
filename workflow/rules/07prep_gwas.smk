localrules: get_gwas_sumstats

rule get_gwas_sumstats:
    # Download GWAS sumsatts files
    output:  config["slsdr"]["get_gwas_sumstats"]["output"]
    params:  lambda wildcards: config['gwas'][wildcards.gwas]
    message: "Download {wildcards.gwas} sumstats file"
    log:     config["slsdr"]["get_gwas_sumstats"]["log"]
    run:

             if wildcards.gwas in ("mdd", "neuroticism"):

                 shell("""

                 cp {params} temp; gunzip -c temp > {output}; rm temp 2> {log}

                 """)
             
             elif wildcards.gwas in ("scz", "bd"):

                 shell("""

                 wget -O - {params} | gunzip -c | sed '/##/d' > {output} 

                 """)

             else:
                 
                 shell("""

                 wget -O - {params} | gunzip -c > {output} 

                 """)


rule standardise_sumstats:
    # Standardises sumstats: SNP, CHR. BP, PVAL, A1, A2 + additional GWAS dependant cols
    # python convert available here: https://github.com/precimed/python_convert/tree/master
    input:   rules.get_gwas_sumstats.output
    output:  config["slsdr"]["standardise_sumstats"]["output"]
    message: "Standardising {input}"
    params: config["slsdr"]["standardise_sumstats"]["temp"]
    log:    config["slsdr"]["standardise_sumstats"]["log"] 
    run:

             if wildcards.gwas in ("bd", "scz"):

                 shell("""

                 cat {input} | sed 's/ID/SNP/g' | sed 's/#CHROM/CHR/g' > {params};
                 python ../resources/python_convert/sumstats.py csv \
                   --sumstats {params} \
                   --out {output} --force --auto --head 5 \
                   --log {log};

                 rm {params}

                  """)

             elif "height" in wildcards.gwas:

                 shell("""
                 cat {input} | sed 's/Tested_Allele/A1/g' | sed 's/Other_Allele/A2/g' > {params};
    
                 python ../resources/python_convert/sumstats.py csv \
                   --sumstats {params} \
       	           --out {output} --force --auto --head 5 \
                   --log {log};

                 rm {params}
             
                  """)

             elif "neuroticism" in wildcards.gwas:

                 shell("""
 
                 cat {input} | sed 's/POS/BP/g' | sed 's/RSID_UKB/SNP/g' | sed 's/REF/A1/g' | sed 's/ALT/A2/g' > {params};

                 python ../resources/python_convert/sumstats.py csv \
                   --sumstats {params} \
                   --out {output} --force --auto --head 5 \
                   --log {log};

                 rm {params}

                 """)

             else:

                 shell("""

                 python ../resources/python_convert/sumstats.py csv \
                   --sumstats {input} \
                   --out {output} --force --auto --head 5 \
                   --log {log}

                  """)

rule add_z_score:
    # Adds z-scores to GWAS sumstats lacking
    input:   rules.standardise_sumstats.output
    output:  config["slsdr"]["add_z_score"]["output"]
    message: "Adding Z score to {input} if required"
    log:     config["slsdr"]["add_z_score"]["log"]
    shell:
             """
             python ../resources/python_convert/sumstats.py zscore \
             --sumstats {input} \
             --out {output} --force \
             --log {log} \
             --a1-inc
             """

rule add_N:
    # N to GWAS sumstats lacking 
    input:   rules.add_z_score.output
    output:  config["slsdr"]["add_N"]["output"]
    message: "Adding N to {input} if required"
    log:     config["slsdr"]["add_N"]["log"]
    run:

             if "scz" in wildcards.gwas:

                 shell("""

                 awk -v OFS='\t' '{{{{s=(NR==1)?"N":"130644"; $0=$0 OFS s}}}}1' {input} > {output} 2> {log}

                 """)

             elif "adhd" in wildcards.gwas:

                 shell("""

                 awk -v OFS='\t' '{{s=(NR==1)?"N":"225534";$0=$0 OFS s}}1' {input} > {output} 2> {log}
 
                 """)

             elif "asd" in wildcards.gwas:

                 shell("""

                 awk -v OFS='\t' '{{s=(NR==1)?"N":"46350";$0=$0 OFS s}}1' {input} > {output} 2> {log}
 
                 """)
             
             elif "bd" in wildcards.gwas:

                 shell("""

                 awk -v OFS='\t' '{{s=(NR==1)?"N":"413466";$0=$0 OFS s}}1' {input} > {output} 2> {log}

                 """)

             elif "neuroticism" in wildcards.gwas:

                 shell("""

                 awk -v OFS='\t' '{{s=(NR==1)?"N":"313467";$0=$0 OFS s}}1' {input} > {output} 2> {log}

                 """)

             else:

                 shell("cp {input} {output}")

rule make_gwas_bed_hg19:
    # Generate gwas hg19 bed file for lift over
    input:   rules.add_N.output
    output:  config["slsdr"]["make_gwas_bed_hg19"]["output"]
    message: "Create bed input file for {input} hg19 to hg38 LiftOver"
    log: config["slsdr"]["make_gwas_bed_hg19"]["log"]
    run:

        shell("""
        
        awk 'BEGIN {{OFS="\t"}} NR > 1 {{print "chr"$2, $3-1, $3, $1}}' {input} > {output} 2> {log}

        """)

rule liftover_gwas_to_hg38:
    # Liftover sumstats from hg19 to hg38
    input: bed = rules.make_gwas_bed_hg19.output,
           chain = config["slsdr"]["liftover_gwas_to_hg38"]["chain"]
    output:  config["slsdr"]["liftover_gwas_to_hg38"]["output"]
    singularity: config["containers"]["ubuntu"]
    params: config["slsdr"]["liftover_gwas_to_hg38"]["unlifted"]
    message: "Lifting over {input} from hg19 to hg38"
    log: config["slsdr"]["liftover_gwas_to_hg38"]["log"]
    shell:
        """

        ../resources/liftover/liftOver {input.bed} {input.chain} {output} {params} 2> {log}
    
        """

rule add_hg38_coords_to_gwas:
    #Restore summary statistics file with hg38 coords
    input: lifted = rules.liftover_gwas_to_hg38.output,
           sumstats = rules.add_N.output
    output: config["slsdr"]["add_hg38_coords_to_gwas"]["output"]
    message: "Adding hg38 coords to {input.sumstats} sumstats"
    log: config["slsdr"]["add_hg38_coords_to_gwas"]["log"]
    script: "../scripts/add_hg38_coords_to_gwas.py"


rule munge_sumstats:
    # Format sumstats for LDSR input
    input:   snps = "../results/05SLDSR/sldsr/1000G.EUR.hg38.w_hm3_test.snplist",
             gwas = rules.add_hg38_coords_to_gwas.output
    output:  config["slsdr"]["munge_sumstats"]["output"]
#    conda:   "../envs/ldsr.yml" # Struggling to get this to work, can create env, but snakemake can't
    message: "Munging sumstats for LDSR compatibility: {input.gwas}"
    params:  config["slsdr"]["munge_sumstats"]["prefix"]
    log:     config["slsdr"]["munge_sumstats"]["log"]
    shell:
        """
        eval "$(/apps/languages/miniforge3/24.3.0-0/bin/conda shell.bash hook)"
        conda activate ldsr
        echo "Starting Snakemake script..." > {log}
        which python >> {log}
        python ../resources/ldsr/ldsc/munge_sumstats.py  --sumstats {input.gwas} \
          --merge-alleles {input.snps} \
          --out {params} \
          --a1-inc \
          --p PVAL >> {log} 2>&1
        echo "Finished at $(date)" >> {log}

        """


