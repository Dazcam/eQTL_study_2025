##  Settings  ------------------------------------------------------------------------------
PG = config["prep_gwas"]


##  Functions  -----------------------------------------------------------------------------
def raw_gwas_file(wildcards):
    gwas_dir = config["setup"]["get_gwas_sumstats"]["dir"].format(gwas = wildcards.gwas)
    return f"{gwas_dir}/{config['gwas'][wildcards.gwas]['file']}"


##  Rules  ---------------------------------------------------------------------------------

rule extract_gwas_sumstats:
    # Decompress the downloaded GWAS file and drop ## header lines
    input:   done = config["setup"]["get_gwas_sumstats"]["output"]
    output:  PG["extract_gwas_sumstats"]["output"]
    wildcard_constraints: gwas = "|".join(config["gwas"])
    params:  raw = raw_gwas_file
    message: "Extract {wildcards.gwas} sumstats file"
    benchmark: PG["extract_gwas_sumstats"]["benchmark"]
    log:     PG["extract_gwas_sumstats"]["log"]
    shell:
             """
             exec 2> {log}
             gunzip -c {params.raw} | sed "/##/d" > {output}
             """

rule standardise_sumstats:
    # Standardises sumstats: SNP, CHR. BP, PVAL, A1, A2 + additional GWAS dependant cols
    # python convert available here: https://github.com/precimed/python_convert/tree/master
    input:   rules.extract_gwas_sumstats.output
    output:  PG["standardise_sumstats"]["output"]
    message: "Standardising {input}"
    benchmark: PG["standardise_sumstats"]["benchmark"]
    params: PG["standardise_sumstats"]["temp"]
    log:    PG["standardise_sumstats"]["log"] 
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
    output:  PG["add_z_score"]["output"]
    params:  script = PG["add_z_score"]["script"]
    message: "Adding Z score to {input} if required"
    benchmark: PG["add_z_score"]["benchmark"]
    log:     PG["add_z_score"]["log"]
    shell:
             """
             python {params.script} \
             --sumstats {input} \
             --out {output} > {log} 2>&1
             """

rule add_N:
    # N to GWAS sumstats lacking 
    input:   rules.add_z_score.output
    output:  PG["add_N"]["output"]
    message: "Adding N to {input} if required"
    benchmark: PG["add_N"]["benchmark"]
    log:     PG["add_N"]["log"]
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
    output:  PG["make_gwas_bed_hg19"]["output"]
    message: "Create bed input file for {input} hg19 to hg38 LiftOver"
    benchmark: PG["make_gwas_bed_hg19"]["benchmark"]
    log: PG["make_gwas_bed_hg19"]["log"]
    run:

        shell("""
        
        awk 'BEGIN {{OFS="\t"}} NR > 1 {{print "chr"$2, $3-1, $3, $1}}' {input} > {output} 2> {log}

        """)

rule liftover_gwas_to_hg38:
    # Liftover sumstats from hg19 to hg38
    input: bed = rules.make_gwas_bed_hg19.output,
           chain = PG["liftover_gwas_to_hg38"]["chain"]
    output:  PG["liftover_gwas_to_hg38"]["output"]
    singularity: config["containers"]["ubuntu"]
    params: PG["liftover_gwas_to_hg38"]["unlifted"]
    benchmark: PG["liftover_gwas_to_hg38"]["benchmark"]
    message: "Lifting over {input} from hg19 to hg38"
    log: PG["liftover_gwas_to_hg38"]["log"]
    shell:
        """

        ../resources/liftover/liftOver {input.bed} {input.chain} {output} {params} 2> {log}
    
        """

rule add_hg38_coords_to_gwas:
    #Restore summary statistics file with hg38 coords
    input: lifted = rules.liftover_gwas_to_hg38.output,
           sumstats = rules.add_N.output
    output: PG["add_hg38_coords_to_gwas"]["output"]
    resources: threads = 5, mem_mb = 20000, time="5:00:00"
    message: "Adding hg38 coords to {input.sumstats} sumstats"
    benchmark: PG["add_hg38_coords_to_gwas"]["benchmark"]
    log: PG["add_hg38_coords_to_gwas"]["log"]
    script: PG["add_hg38_coords_to_gwas"]["script"]


rule munge_sumstats:
    # Format sumstats for LDSR input
    input:   snps = PG["munge_sumstats"]["snplist"],
             gwas = rules.add_hg38_coords_to_gwas.output
    output:  PG["munge_sumstats"]["output"]
#    conda:   "../envs/ldsr.yml" # Struggling to get this to work, can create env, but snakemake can't
    message: "Munging sumstats for LDSR compatibility: {input.gwas}"
    benchmark: PG["munge_sumstats"]["benchmark"]
    params:  PG["munge_sumstats"]["prefix"]
    log:     PG["munge_sumstats"]["log"]
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


##  Final outputs, collected by rule all in the Snakefile  ---------------------------------
TARGETS = globals().setdefault("TARGETS", {})
TARGETS["09prep_gwas"] = [expand(rules.munge_sumstats.output, gwas = config["gwas"])]
