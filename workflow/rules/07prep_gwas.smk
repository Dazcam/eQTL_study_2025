configfile: '../config/config.yaml' 

localrules: get_gwas_sumstats

#rule all:
#    input:
#        expand("../results/07PREP-GWAS/{gwas}_hg38_ldsr_ready.sumstats.gz", gwas = config['gwas'])

rule get_gwas_sumstats:
    # Download GWAS sumsatts files
    output:  config["prep_gwas"]["get_gwas_sumstats"]["output"]
    params:  lambda wildcards: config['gwas'][wildcards.gwas]
    message: "Download {wildcards.gwas} sumstats file"
    benchmark: "reports/benchmarks/07prep_gwas.get_gwas_sumstats_{gwas}.txt"
    log:     config["prep_gwas"]["get_gwas_sumstats"]["log"]
    run:
        
             user_agent = '"Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/91.0.4472.124 Safari/537.36"'

             if wildcards.gwas in ("scz", "bpd", "mdd"):

                 shell("""

                 wget -U {user_agent} -O - {params} | gunzip -c | sed '/##/d' > {output} 2> {log}

                 """)

             elif wildcards.gwas == "adhd":
            
                 shell("""
                 
                 wget -U {user_agent} -O temp_adhd.zip {params} 2> {log}
                 unzip -p temp_adhd.zip ADHD2022_iPSYCH_deCODE_PGC.meta.gz | gunzip -c | sed '/##/d' > {output} 2>> {log}
                 rm -f temp_adhd.zip
                 
                 """)

             elif wildcards.gwas == "ptsd":

                 shell("""

                 wget -U {user_agent} -O temp_ptsd.zip {params} 2> {log}
                 unzip -p temp_ptsd.zip eur_ptsd_pcs_v4_aug3_2021.vcf.gz | gunzip -c | sed '/##/d' > {output} 2>> {log}
                 rm -f temp_ptsd.zip
       	       	 
                 """)

             elif wildcards.gwas == "ocd":

                 shell("""

                 wget -U {user_agent} -O temp_ocd.zip {params} 2> {log}
                 unzip -p temp_ocd.zip daner_OCD_full_wo23andMe_190522.gz | gunzip -c | sed '/##/d' > {output} 2>> {log}
                 rm -f temp_ocd.zip

                 """)

             else:
                 shell("""
                 
                 wget -U {user_agent} -O - {params} | gunzip -c > {output} 2> {log}
                  
                 """)

rule standardise_sumstats:
    # Standardises sumstats: SNP, CHR. BP, PVAL, A1, A2 + additional GWAS dependant cols
    # python convert available here: https://github.com/precimed/python_convert/tree/master
    input:   rules.get_gwas_sumstats.output
    output:  config["prep_gwas"]["standardise_sumstats"]["output"]
    message: "Standardising {input}"
    benchmark: "reports/benchmarks/07prep_gwas.standardise_sumstats_{gwas}.txt"
    params: config["prep_gwas"]["standardise_sumstats"]["temp"]
    log:    config["prep_gwas"]["standardise_sumstats"]["log"] 
    run:

             if wildcards.gwas in ("mdd", "scz"):

                 shell("""

                 cat {input} | sed 's/ID/SNP/g' | sed 's/#CHROM/CHR/g' | sed 's/CHROM/CHR/g' > {params};
                 python ../resources/python_convert/sumstats.py csv \
                   --sumstats {params} \
                   --out {output} --force --auto --head 5 \
                   --keep-cols NCAS NCON \
                   --log {log};

                 rm {params}

                  """)

             elif "bpd" in wildcards.gwas:

                 shell("""
                 cat {input} | sed 's/EFFECT_ALLELE/A1/g' | sed 's/OTHER_ALLELE/A2/g' > {params};
    
                 python ../resources/python_convert/sumstats.py csv \
                   --sumstats {params} \
       	           --out {output} --force --auto --head 5 \
                   --keep-cols NCA NCO \
                   --log {log};

                 rm {params}
             
                  """)

             elif wildcards.gwas in ("adhd", "ocd"):

                 shell("""
              
                 cat {input} | sed 's/Nca/NCA/g' | sed 's/Nco/NCO/g' > {params};
                 python ../resources/python_convert/sumstats.py csv \
                   --sumstats {params} \
                   --out {output} --force --auto --head 5 \
                   --keep-cols NCA NCO \
                   --log {log};

                 rm {params}

                  """)

             elif "ptsd" in wildcards.gwas:

                 shell("""

                 cat {input} | sed 's/NEFF/N/g' | sed 's/ID/SNP/g' | sed 's/#CHROM/CHR/g' | sed 's/POS/BP/g'> {params};
                 python ../resources/python_convert/sumstats.py csv \
                   --sumstats {params} \
                   --out {output} --force --auto --head 5 \
                   --keep-cols N \
                   --log {log};

                 rm {params}

                  """)


             else:

                 shell("""

                 python ../resources/python_convert/sumstats.py csv \
                   --sumstats {input} \
                   --out {output} --force --auto --head 5 \
                   --keep-cols N \
                   --log {log}

                  """)

rule add_z_score:
    # Adds z-scores to GWAS sumstats lacking
    input:   rules.standardise_sumstats.output
    output:  config["prep_gwas"]["add_z_score"]["output"]
    message: "Adding Z score to {input} if required"
    benchmark: "reports/benchmarks/07prep_gwas.add_z_score_{gwas}.txt"
    log:     config["prep_gwas"]["add_z_score"]["log"]
    shell:
             """
             python scripts/add_z_score_to_sumstats.py \
             --sumstats {input} \
             --out {output} > {log} 2>&1
             """

rule add_N:
    # N to GWAS sumstats lacking 
    input:   rules.add_z_score.output
    output:  config["prep_gwas"]["add_N"]["output"]
    message: "Adding N to {input} if required"
    benchmark: "reports/benchmarks/07prep_gwas.add_N_{gwas}.txt"
    log:     config["prep_gwas"]["add_N"]["log"]
    run:
        if wildcards.gwas in ("mdd", "scz"):
            
            shell("""
            
            awk -v OFS='\t' 'BEGIN{{FS=OFS="\t"}} NR==1{{for(i=1;i<=NF;i++) {{if($i=="NCAS") c1=i; if($i=="NCON") c2=i}} print $0,"N"}} NR>1{{print $0,$c1+$c2}}' {input} > {output} 2> {log}
            
            """)
        
        elif wildcards.gwas in ("bpd", "adhd", "ocd"):
            
            shell("""
            
            awk -v OFS='\t' 'BEGIN{{FS=OFS="\t"}} NR==1{{for(i=1;i<=NF;i++) {{if($i=="NCA") c1=i; if($i=="NCO") c2=i}} print $0,"N"}} NR>1{{print $0,$c1+$c2}}' {input} > {output} 2> {log}
            
            """)
        
        else:
            
            shell("""
            
            cp {input} {output} 2> {log}
            
            """)

rule make_gwas_bed_hg19:
    # Generate gwas hg19 bed file for lift over
    input:   rules.add_N.output
    output:  config["prep_gwas"]["make_gwas_bed_hg19"]["output"]
    message: "Create bed input file for {input} hg19 to hg38 LiftOver"
    benchmark: "reports/benchmarks/07prep_gwas.make_gwas_bed_hg19_{gwas}.txt"
    log: config["prep_gwas"]["make_gwas_bed_hg19"]["log"]
    run:

        shell("""
        
        awk 'BEGIN {{OFS="\t"}} NR > 1 {{print "chr"$2, $3-1, $3, $1}}' {input} > {output} 2> {log}

        """)

rule liftover_gwas_to_hg38:
    # Liftover sumstats from hg19 to hg38
    input: bed = rules.make_gwas_bed_hg19.output,
           chain = config["prep_gwas"]["liftover_gwas_to_hg38"]["chain"]
    output:  config["prep_gwas"]["liftover_gwas_to_hg38"]["output"]
    singularity: config["containers"]["ubuntu"]
    params: config["prep_gwas"]["liftover_gwas_to_hg38"]["unlifted"]
    benchmark: "reports/benchmarks/07prep_gwas.liftover_gwas_to_hg38_{gwas}.txt"
    message: "Lifting over {input} from hg19 to hg38"
    log: config["prep_gwas"]["liftover_gwas_to_hg38"]["log"]
    shell:
        """

        ../resources/liftover/liftOver {input.bed} {input.chain} {output} {params} 2> {log}
    
        """

rule add_hg38_coords_to_gwas:
    #Restore summary statistics file with hg38 coords
    input: lifted = rules.liftover_gwas_to_hg38.output,
           sumstats = rules.add_N.output
    output: config["prep_gwas"]["add_hg38_coords_to_gwas"]["output"]
    resources: threads = 5, mem_mb = 20000, time="5:00:00"
    message: "Adding hg38 coords to {input.sumstats} sumstats"
    benchmark: "reports/benchmarks/07prep_gwas.add_hg38_coords_to_gwas_{gwas}.txt"
    log: config["prep_gwas"]["add_hg38_coords_to_gwas"]["log"]
    script: "../scripts/add_hg38_coords_to_gwas.py"


rule munge_sumstats:
    # Format sumstats for LDSR input
    input:   snps = "../resources/ldsr/ldsr_hg38_refs/1000G.EUR.hg38.w_hm3_test.snplist",
             gwas = rules.add_hg38_coords_to_gwas.output
    output:  config["prep_gwas"]["munge_sumstats"]["output"]
#    conda:   "../envs/ldsr.yml" # Struggling to get this to work, can create env, but snakemake can't
    message: "Munging sumstats for LDSR compatibility: {input.gwas}"
    benchmark: "reports/benchmarks/07prep_gwas.munge_sumstats_{gwas}.txt"
    params:  config["prep_gwas"]["munge_sumstats"]["prefix"]
    log:     config["prep_gwas"]["munge_sumstats"]["log"]
    shell:
        """
        eval "$(/apps/languages/miniforge3/24.3.0-0/bin/conda shell.bash hook)"
        conda activate ldsr
        echo "Starting Snakemake script..." > {log}
        which python >> {log}
        python ../resources/ldsr/ldsc/munge_sumstats.py  --sumstats {input.gwas} \
          --merge-alleles {input.snps} \
          --out {params} \
          --signed-sumstats Z,0 \
          --p PVAL >> {log} 2>&1
        echo "Finished at $(date)" >> {log}
        """


