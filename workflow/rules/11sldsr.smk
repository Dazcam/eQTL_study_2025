import os

localrules: ldsr_strat_summary


##  Settings  ------------------------------------------------------------------------------
LDS = config["sldsr"]

# The 00 setup get_hg38_refs / get_hg19_refs config patterns are the single source of truth
# for reference paths; the directories and ldsc prefixes below are derived from them.
HG38_REFS = config["setup"]["get_hg38_refs"]
HG19_REFS = config["setup"]["get_hg19_refs"]


##  Functions  -----------------------------------------------------------------------------
def strip_suffix(path, suffix):
    """Remove suffix from a config path pattern, failing loudly if it isn't there."""
    if not path.endswith(suffix):
        raise ValueError(f"Expected '{path}' to end with '{suffix}'")
    return path[: -len(suffix)]

def get_annot(wildcards):
    """Annotation file for the annot_type wildcard (maxCPP or CS95)."""
    key = {"maxCPP": "maxcpp", "CS95": "cs95"}[wildcards.annot_type]
    return LDS["make_annot"][key].format(cell_type=wildcards.cell_type,
                                         chr=wildcards.chr)


##  Rules  ---------------------------------------------------------------------------------

rule make_annot:
    input:  cred_set = config["susie"]["sort_susie"]["output"],
            bim = HG38_REFS["bim"]
    output: maxcpp = LDS["make_annot"]["maxcpp"],
            cs95 = LDS["make_annot"]["cs95"]
    log:    LDS["make_annot"]["log"]
    singularity: config["containers"]["r_eqtl"]
    benchmark: LDS["make_annot"]["benchmark"]
    message: "Generating maxCPP and CS95 annotations for {wildcards.cell_type}, chr {wildcards.chr}"
    script: LDS["make_annot"]["script"]

rule lift_hapmap3_snps:
    # Lift HapMap3 SNPs to hg38
    input:  snps = HG19_REFS["snps"],
            bim_hg19 = expand(HG19_REFS["bim"], chr=range(1, 23)),
            bim_hg38 = expand(HG38_REFS["bim"], chr=range(1, 23))
    output: LDS["lift_hapmap3_snps"]["output"]
    params: hg19_dir = os.path.dirname(HG19_REFS["bim"]),
            hg38_dir = os.path.dirname(HG38_REFS["bim"]),
            script = LDS["lift_hapmap3_snps"]["script"]
    singularity: config["containers"]["ubuntu"]
    benchmark: LDS["lift_hapmap3_snps"]["benchmark"]
    log:    LDS["lift_hapmap3_snps"]["log"]
    shell:
        """
        {params.script} {output} {input.snps} \
          {params.hg19_dir} {params.hg38_dir} > {log} 2>&1
        """

rule ldsr_ld_scores_hg38:
    input:  annot = get_annot,
            bfile = expand(HG38_REFS["bim"], chr=range(1, 23)),
            snps = rules.lift_hapmap3_snps.output
    output: LDS["ldsr_ld_scores_hg38"]["output"]
    params: bfile = strip_suffix(HG38_REFS["bim"], ".bim"),
            ldscores = strip_suffix(LDS["ldsr_ld_scores_hg38"]["output"],
                                    ".l2.ldscore.gz"),
            script = config["tools"]["ldsc"],
            env = config["tools"]["ldsr_env"]
    envmodules: config["envmodules"]["miniforge"]
    benchmark: LDS["ldsr_ld_scores_hg38"]["benchmark"]
    message: "Generating LD scores on hg38 for {wildcards.cell_type}, {wildcards.annot_type}, chr {wildcards.chr}"
    log:    LDS["ldsr_ld_scores_hg38"]["log"]
    shell:
        """
        eval "$(conda shell.bash hook)"
        conda activate {params.env}
        python {params.script} --l2 \
          --bfile {params.bfile} --ld-wind-cm 1 \
          --annot {input.annot} \
          --out {params.ldscores} \
          --print-snps {input.snps} > {log} 2>&1
        """

rule ldsr_strat_hg38_bl_v12:
    input:  gwas = config["prep_gwas"]["munge_sumstats"]["output"],
            ldsr = expand(rules.ldsr_ld_scores_hg38.output, chr=range(1, 23),
                          allow_missing=True)
    output: LDS["ldsr_strat_hg38_bl_v12"]["output"]
    params: weights = strip_suffix(HG38_REFS["weights"], "{chr}.l2.ldscore.gz"),
            baseline = strip_suffix(HG38_REFS["baseline"], "{chr}.l2.ldscore.gz"),
            frq = strip_suffix(HG38_REFS["frq"], "{chr}.frq"),
            ldscores = strip_suffix(LDS["ldsr_ld_scores_hg38"]["output"],
                                    "{chr}.l2.ldscore.gz"),
            out_prefix = strip_suffix(LDS["ldsr_strat_hg38_bl_v12"]["output"],
                                      ".results"),
            script = config["tools"]["ldsc"],
            env = config["tools"]["ldsr_env"]
    envmodules: config["envmodules"]["miniforge"]
    message: "Running stratified LDSR with {wildcards.cell_type} using hg38 refs, baseline 1.2, {wildcards.annot_type} and {wildcards.gwas} GWAS"
    benchmark: LDS["ldsr_strat_hg38_bl_v12"]["benchmark"]
    log:    LDS["ldsr_strat_hg38_bl_v12"]["log"]
    shell:
        """
        eval "$(conda shell.bash hook)"
        conda activate {params.env}
        echo "Starting Snakemake script..." > {log}
        which python >> {log}

        python {params.script} --h2 {input.gwas} \
          --w-ld-chr {params.weights} \
          --ref-ld-chr {params.baseline},{params.ldscores} \
          --overlap-annot \
          --frqfile-chr {params.frq} \
          --out {params.out_prefix} \
          --print-coefficients >> {log} 2>&1
        """

rule ldsr_strat_summary:
    input:  expand(rules.ldsr_strat_hg38_bl_v12.output, cell_type=config["cell_types"], gwas=config["gwas"], allow_missing=True)
    output: LDS["ldsr_strat_summary"]["output"]
    benchmark: LDS["ldsr_strat_summary"]["benchmark"]
    log:    LDS["ldsr_strat_summary"]["log"]
    shell:
        """
        exec > {log} 2>&1
        head -1 {input[0]} > {output}
        grep L2_1 {input} | cut -f 6- -d/ >> {output}
        """

rule ldsr_report:
    input:  ldsr = expand(rules.ldsr_strat_summary.output, annot_type=config["annot_types"]),
            rmd_script = LDS["ldsr_report"]["script"]
    output: LDS["ldsr_report"]["html"]
    params: in_dir = LDS["ldsr_report"]["in_dir"],
            output_file = LDS["ldsr_report"]["output_file"],
            bmark_dir = LDS["ldsr_report"]["bmark_dir"],
            cell_types = ",".join(config["cell_types"]),
            gwas = ",".join(config["gwas"]),
            annot_types = ",".join(config["annot_types"]),
            ld_baseline = os.path.dirname(HG38_REFS["baseline"]),
            ld_weights = os.path.dirname(HG38_REFS["weights"]),
            ld_frq = os.path.dirname(HG38_REFS["frq"]),
            ld_snps = LDS["lift_hapmap3_snps"]["output"]
    singularity: config["containers"]["r_eqtl"]
    message: "Generate SLDSR report"
    benchmark: LDS["ldsr_report"]["benchmark"]
    log:    LDS["ldsr_report"]["log"]
    shell:
        """
        Rscript -e "rmarkdown::render('{input.rmd_script}', \
            output_file = '{params.output_file}', \
            params = list(in_dir = '{params.in_dir}', \
                          bmark_dir = '{params.bmark_dir}', \
                          cell_types = '{params.cell_types}', \
                          gwas = '{params.gwas}', \
                          annot_types = '{params.annot_types}', \
                          ld_baseline = '{params.ld_baseline}', \
                          ld_weights = '{params.ld_weights}', \
                          ld_frq = '{params.ld_frq}', \
                          ld_snps = '{params.ld_snps}'))" > {log} 2>&1
        """


##  Final outputs, collected by rule all in the Snakefile  ---------------------------------
TARGETS = globals().setdefault("TARGETS", {})
TARGETS["11sldsr"] = [LDS["ldsr_report"]["html"]]
