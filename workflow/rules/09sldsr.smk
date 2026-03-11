configfile: "../config/config.yaml"
localrules: ldsr_strat_summary

rule all:
    input:
        "reports/09SLDSR/09ldsr_report.html"

rule get_hg38_refs:
    output: baseline = expand(config["sldsr"]["get_hg38_refs"]["baseline"], chr=range(1, 23)),
            weights = expand(config["sldsr"]["get_hg38_refs"]["weights"], chr=range(1, 23)),
            frq = expand(config["sldsr"]["get_hg38_refs"]["frq"], chr=range(1, 23)),
            bim = expand(config["sldsr"]["get_hg38_refs"]["bim"], chr=range(1, 23))
#            snps = config["sldsr"]["get_hg38_refs"]["snps"]
    params: config["sldsr"]["get_hg38_refs"]["params"]
    benchmark: "reports/benchmarks/09sldsr.get_hg38_refs.txt"
    log:    config["sldsr"]["get_hg38_refs"]["log"]
    shell: "scripts/ldsr_get_hg38_refs.sh {params} > {log} 2>&1"

rule get_hg19_refs:
    output: baseline = expand(config["sldsr"]["get_hg19_refs"]["baseline"], chr=range(1, 23)),
            weights = expand(config["sldsr"]["get_hg19_refs"]["weights"], chr=range(1, 23)),
            frq = expand(config["sldsr"]["get_hg19_refs"]["frq"], chr=range(1, 23)),
            bim = expand(config["sldsr"]["get_hg19_refs"]["bim"], chr=range(1, 23))
#            snps = config["sldsr"]["get_hg19_refs"]["snps"]
    params: config["sldsr"]["get_hg19_refs"]["params"]
    benchmark: "reports/benchmarks/09sldsr.get_hg19_refs.txt"
    log:    config["sldsr"]["get_hg19_refs"]["log"]
    shell: "scripts/ldsr_get_hg19_refs.sh {params} > {log} 2>&1"

rule make_annot:
    input:  cred_set = config["susie"]["sort_susie"]["output"],
            bim = lambda wildcards: config["sldsr"]["get_hg38_refs"]["bim"].format(chr=wildcards.chr)
    output: maxcpp = config["sldsr"]["make_annot"]["maxcpp"],
            cs95 = config["sldsr"]["make_annot"]["cs95"]
    log:    config["sldsr"]["make_annot"]["log"]
    singularity: config["containers"]["r_eqtl"]
    benchmark: "reports/benchmarks/09sldsr.make_annot_{chr}_{cell_type}.txt"
    message: "Generating maxCPP and CS95 annotations for {wildcards.cell_type}, chr {wildcards.chr}"
    script: "../scripts/ldsr_make_annot.R"

rule lift_hapmap3_snps:
   # Lift hapmap3 SNPs to hg38
   input:  bim_hg38 = lambda wildcards: expand(rules.get_hg38_refs.output.bim, chr=range(1,23)),
           bim_hg19 = lambda wildcards: expand(rules.get_hg19_refs.output.bim, chr=range(1,23))
   output: config["sldsr"]["lift_hapmap3_snps"]["output"],
   singularity: config["containers"]["ubuntu"]   
   log: config["sldsr"]["lift_hapmap3_snps"]["log"]
   shell: "/scratch/c.c1477909/eQTL_study_2025/workflow/scripts/liftover_hapmap3_snps.sh {output} > {log} 2>&1"


rule ldsr_ld_scores_hg38:
    input:  annot = lambda wildcards: config["sldsr"]["make_annot"]["maxcpp"].format(cell_type=wildcards.cell_type, chr=wildcards.chr) if wildcards.annot_type == "maxCPP" else config["sldsr"]["make_annot"]["cs95"].format(cell_type=wildcards.cell_type, chr=wildcards.chr),
            bfile = rules.get_hg38_refs.output.bim,
            snps = rules.lift_hapmap3_snps.output
    output: config["sldsr"]["ldsr_ld_scores_hg38"]["output"]
    params: bfile = config["sldsr"]["ldsr_ld_scores_hg38"]["bfile"],
            ldscores = config["sldsr"]["ldsr_ld_scores_hg38"]["ldscores"]
    message: "Generating LD scores on hg38 for {wildcards.cell_type}, {wildcards.annot_type}, chr {wildcards.chr}" 
    log: config["sldsr"]["ldsr_ld_scores_hg38"]["log"]
    shell:
        """
        eval "$(/apps/languages/miniforge3/24.3.0-0/bin/conda shell.bash hook)"
        conda activate ldsr
        python ../resources/ldsr/ldsc/ldsc.py --l2 \
          --bfile {params.bfile} --ld-wind-cm 1 \
          --annot {input.annot} \
          --out {params.ldscores} \
          --print-snps {input.snps} 2> {log} 2>&1
        """

rule ldsr_strat_hg38_bl_v12:
    input:  gwas = config["prep_gwas"]["munge_sumstats"]["output"],
            ldsr = lambda wildcards: expand(rules.ldsr_ld_scores_hg38.output, 
                                            cell_type=config["cell_types"], 
                                            annot_type=config["annot_types"], 
                                            chr=range(1,23))
    output: config["sldsr"]["ldsr_strat_hg38_bl_v12"]["output"]
    params: weights = config["sldsr"]["ldsr_strat_hg38_bl_v12"]["weights"],
            baseline = config["sldsr"]["ldsr_strat_hg38_bl_v12"]["baseline"],
            frq = config["sldsr"]["ldsr_strat_hg38_bl_v12"]["frq"],
            ldscores = config["sldsr"]["ldsr_strat_hg38_bl_v12"]["ldscores"],
            out_prefix = config["sldsr"]["ldsr_strat_hg38_bl_v12"]["out_prefix"]
    message: "Running stratified LDSR with {wildcards.cell_type} using hg38 refs, baseline 1.2, {wildcards.annot_type} and {wildcards.gwas} GWAS"
    log: config["sldsr"]["ldsr_strat_hg38_bl_v12"]["log"]
    shell:
        """
        eval "$(/apps/languages/miniforge3/24.3.0-0/bin/conda shell.bash hook)"
        conda activate ldsr
        echo "Starting Snakemake script..." > {log}
        which python >> {log}

        python ../resources/ldsr/ldsc/ldsc.py --h2 {input.gwas} \
          --w-ld-chr {params.weights} \
          --ref-ld-chr {params.baseline},{params.ldscores} \
          --overlap-annot \
          --frqfile-chr {params.frq} \
          --out {params.out_prefix} \
          --print-coefficients >> {log} 2>&1
       """

rule ldsr_strat_summary:
    input: lambda wildcards: expand(rules.ldsr_strat_hg38_bl_v12.output, cell_type=config["cell_types"], annot_type=wildcards.annot_type, gwas=config["gwas"])
    output: config["sldsr"]["ldsr_strat_summary"]["output"]
    log: config["sldsr"]["ldsr_strat_summary"]["log"]
    shell: """
           head -1 {input[0]} > {output}
           grep L2_1 {input} | cut -f 6- -d/ >> {output}
           """

rule ldsr_report:
    input:  ldsr = lambda wildcards: expand(rules.ldsr_strat_summary.output, annot_type=config["annot_types"]),
            rmd_script = "scripts/ldsr_report.Rmd"
    output: config["sldsr"]["ldsr_report"]["html"]
    params: in_dir = config["sldsr"]["ldsr_report"]["in_dir"],
#            in_files = lambda wildcards: " ".join(expand("../../results/09SLDSR/strat_bl_v12/ldsr_strat_hg38_bl_v12.{annot_type}.summary.tsv", annot_type=config["annot_types"])),
            output_file = "../reports/09SLDSR/09ldsr_report.html"
    singularity: config["containers"]["r_eqtl"]
    message: "Generate SLDSR report"
    benchmark: "reports/benchmarks/09sldsr.sldsr_report.txt"
    log: config["sldsr"]["ldsr_report"]["log"]
    shell:
        """
        Rscript -e "rmarkdown::render('{input.rmd_script}', \
            output_file = '{params.output_file}', \
            params = list(in_dir = c('{params.in_dir}')))" > {log} 2>&1
        """
