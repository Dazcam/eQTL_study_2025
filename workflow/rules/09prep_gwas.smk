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
    # Standardise sumstats: SNP, CHR, BP, PVAL, A1, A2 + additional GWAS dependent cols
    # python_convert: https://github.com/precimed/python_convert/tree/master
    input:   rules.extract_gwas_sumstats.output
    output:  PG["standardise_sumstats"]["output"]
    params:  tmp = PG["standardise_sumstats"]["temp"],
             script = config["tools"]["python_convert"]
    message: "Standardising {input}"
    benchmark: PG["standardise_sumstats"]["benchmark"]
    log:     PG["standardise_sumstats"]["log"]
    run:
             if wildcards.gwas in ("mdd", "scz"):
                 shell("""
                 sed '1s/ID/SNP/g; 1s/#CHROM/CHR/g; 1s/CHROM/CHR/g' {input} > {params.tmp};
                 python {params.script} csv \
                   --sumstats {params.tmp} \
                   --out {output} --force --auto --head 5 \
                   --keep-cols NCAS NCON \
                   --log {log};
                 rm {params.tmp}
                 """)

             elif wildcards.gwas == "bpd":
                 shell("""
                 python {params.script} csv \
                   --sumstats {input} \
                   --out {output} --force --auto --head 5 \
                   --keep-cols NCA NCO \
                   --log {log}
                 """)

             elif wildcards.gwas in ("adhd", "ocd"):
                 shell("""
                 sed '1s/Nca/NCA/g; 1s/Nco/NCO/g' {input} > {params.tmp};
                 python {params.script} csv \
                   --sumstats {params.tmp} \
                   --out {output} --force --auto --head 5 \
                   --keep-cols NCA NCO \
                   --log {log};
                 rm {params.tmp}
                 """)

             else:
                 raise ValueError(f"No standardise_sumstats branch for GWAS {wildcards.gwas}")

rule add_z_score:
    # Add Z to sumstats that lack it
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
    # Add N = cases + controls
    input:   rules.add_z_score.output
    output:  PG["add_N"]["output"]
    message: "Adding N to {input}"
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
                 raise ValueError(f"No add_N branch for GWAS {wildcards.gwas}")

rule make_gwas_bed_hg19:
    # Generate GWAS hg19 bed file for liftOver
    input:   rules.add_N.output
    output:  PG["make_gwas_bed_hg19"]["output"]
    message: "Create bed input file for {input} hg19 to hg38 LiftOver"
    benchmark: PG["make_gwas_bed_hg19"]["benchmark"]
    log:     PG["make_gwas_bed_hg19"]["log"]
    shell:
             """
             awk 'BEGIN {{OFS="\t"}} NR > 1 {{print "chr"$2, $3-1, $3, $1}}' {input} > {output} 2> {log}
             """

rule liftover_gwas_to_hg38:
    # Liftover sumstats from hg19 to hg38
    input:   bed = rules.make_gwas_bed_hg19.output,
             chain = PG["liftover_gwas_to_hg38"]["chain"]
    output:  PG["liftover_gwas_to_hg38"]["output"]
    singularity: config["containers"]["ubuntu"]
    params:  unlifted = PG["liftover_gwas_to_hg38"]["unlifted"],
             liftover = config["tools"]["liftover"]
    message: "Lifting over {input} from hg19 to hg38"
    benchmark: PG["liftover_gwas_to_hg38"]["benchmark"]
    log:     PG["liftover_gwas_to_hg38"]["log"]
    shell:
             """
             {params.liftover} {input.bed} {input.chain} {output} {params.unlifted} 2> {log}
             """

rule add_hg38_coords_to_gwas:
    # Restore sumstats file with hg38 coords
    input:   lifted = rules.liftover_gwas_to_hg38.output,
             sumstats = rules.add_N.output
    output:  PG["add_hg38_coords_to_gwas"]["output"]
    resources: mem_mb = 8000, time="1:00:00"
    message: "Adding hg38 coords to {input.sumstats} sumstats"
    benchmark: PG["add_hg38_coords_to_gwas"]["benchmark"]
    log:     PG["add_hg38_coords_to_gwas"]["log"]
    script:  PG["add_hg38_coords_to_gwas"]["script"]

rule munge_sumstats:
    # Format sumstats for LDSR input
    input:   snps = config["setup"]["make_hm3_snplist"]["output"],
             gwas = rules.add_hg38_coords_to_gwas.output
    output:  PG["munge_sumstats"]["output"]
    envmodules: config["envmodules"]["miniforge"]
    params:  prefix = PG["munge_sumstats"]["prefix"],
             script = config["tools"]["ldsc_munge"],
             env = config["tools"]["ldsr_env"]
    message: "Munging sumstats for LDSR compatibility: {input.gwas}"
    benchmark: PG["munge_sumstats"]["benchmark"]
    log:     PG["munge_sumstats"]["log"]
    shell:
             """
             eval "$(conda shell.bash hook)"
             conda activate {params.env}
             echo "Starting Snakemake script..." > {log}
             which python >> {log}
             python {params.script} --sumstats {input.gwas} \
               --merge-alleles {input.snps} \
               --out {params.prefix} \
               --signed-sumstats Z,0 \
               --p PVAL >> {log} 2>&1
             echo "Finished at $(date)" >> {log}
             """


rule prep_gwas_report:
    input:   munged = expand(rules.munge_sumstats.output, gwas = config["gwas"]),
             rmd_script = PG["prep_gwas_report"]["script"]
    output:  PG["prep_gwas_report"]["html"]
    params:  gwas = ",".join(config["gwas"]),
             gwas_files = ",".join(config["gwas"][g]["file"] for g in config["gwas"]),
             gwas_md5 = ",".join(config["gwas"][g]["md5"] for g in config["gwas"]),
             unlifted = lambda w: PG["liftover_gwas_to_hg38"]["unlifted"],
             chain = PG["liftover_gwas_to_hg38"]["chain"],
             snplist = config["setup"]["make_hm3_snplist"]["output"],
             log_dir = PG["prep_gwas_report"]["log_dir"],
             bmark_dir = PG["prep_gwas_report"]["bmark_dir"]
    singularity: config["containers"]["r_eqtl"]
    message: "Generate GWAS prep report"
    benchmark: PG["prep_gwas_report"]["benchmark"]
    log:     PG["prep_gwas_report"]["log"]
    shell:
             """
             Rscript -e "wd <- getwd(); rmarkdown::render('{input.rmd_script}', \
                 output_file = file.path(wd, '{output}'), \
                 knit_root_dir = wd, \
                 params = list(gwas = '{params.gwas}', \
                               gwas_files = '{params.gwas_files}', \
                               gwas_md5 = '{params.gwas_md5}', \
                               unlifted = '{params.unlifted}', \
                               chain = '{params.chain}', \
                               snplist = '{params.snplist}', \
                               log_dir = '{params.log_dir}', \
                               bmark_dir = '{params.bmark_dir}'))" > {log} 2>&1
             """


##  Final outputs, collected by rule all in the Snakefile  ---------------------------------
TARGETS = globals().setdefault("TARGETS", {})
TARGETS["09prep_gwas"] = [expand(rules.munge_sumstats.output, gwas = config["gwas"]),
                          rules.prep_gwas_report.output]
