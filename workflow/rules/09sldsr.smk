configfile: "../config/config.yaml"
localrules: ldsr_strat_summary

import os

##  Reference paths  ------------------------------------------------------------------
# The get_hg38_refs / get_hg19_refs config patterns are the single source of truth for
# reference paths; the directories and ldsc prefixes below are derived from them.

def strip_suffix(path, suffix):
    """Remove suffix from a config path pattern, failing loudly if it isn't there."""
    if not path.endswith(suffix):
        raise ValueError(f"Expected '{path}' to end with '{suffix}'")
    return path[: -len(suffix)]

def ref_dir(build):
    """Top-level reference directory for a build (parent of baseline_v1.2/)."""
    return os.path.dirname(os.path.dirname(config["sldsr"][f"get_{build}_refs"]["baseline"]))

def get_annot(wildcards):
    """Annotation file for the annot_type wildcard (maxCPP or CS95)."""
    key = {"maxCPP": "maxcpp", "CS95": "cs95"}[wildcards.annot_type]
    return config["sldsr"]["make_annot"][key].format(cell_type=wildcards.cell_type,
                                                     chr=wildcards.chr)

HG38_REFS = config["sldsr"]["get_hg38_refs"]
HG19_REFS = config["sldsr"]["get_hg19_refs"]

##  Rules  ----------------------------------------------------------------------------
rule get_hg38_refs:
    output: baseline = expand(HG38_REFS["baseline"], chr=range(1, 23)),
            weights = expand(HG38_REFS["weights"], chr=range(1, 23)),
            frq = expand(HG38_REFS["frq"], chr=range(1, 23)),
            bim = expand(HG38_REFS["bim"], chr=range(1, 23)),
            snps = HG38_REFS["snps"]
    params: out_dir = ref_dir("hg38")
    benchmark: "reports/benchmarks/09sldsr.get_hg38_refs.benchmark.txt"
    log:    HG38_REFS["log"]
    shell: "scripts/ldsr_get_hg38_refs.sh {params.out_dir} > {log} 2>&1"

rule get_hg19_refs:
    output: baseline = expand(HG19_REFS["baseline"], chr=range(1, 23)),
            weights = expand(HG19_REFS["weights"], chr=range(1, 23)),
            frq = expand(HG19_REFS["frq"], chr=range(1, 23)),
            bim = expand(HG19_REFS["bim"], chr=range(1, 23)),
            snps = HG19_REFS["snps"]
    params: out_dir = ref_dir("hg19")
    benchmark: "reports/benchmarks/09sldsr.get_hg19_refs.benchmark.txt"
    log:    HG19_REFS["log"]
    shell: "scripts/ldsr_get_hg19_refs.sh {params.out_dir} > {log} 2>&1"

rule make_annot:
    input:  cred_set = config["susie"]["sort_susie"]["output"],
            bim = HG38_REFS["bim"]
    output: maxcpp = config["sldsr"]["make_annot"]["maxcpp"],
            cs95 = config["sldsr"]["make_annot"]["cs95"]
    log:    config["sldsr"]["make_annot"]["log"]
    singularity: config["containers"]["r_eqtl"]
    benchmark: "reports/benchmarks/09sldsr.make_annot_{cell_type}_{chr}.benchmark.txt"
    message: "Generating maxCPP and CS95 annotations for {wildcards.cell_type}, chr {wildcards.chr}"
    script: "../scripts/ldsr_make_annot.R"

rule lift_hapmap3_snps:
    # Lift HapMap3 SNPs to hg38
    input:  snps = rules.get_hg19_refs.output.snps,
            bim_hg19 = rules.get_hg19_refs.output.bim,
            bim_hg38 = rules.get_hg38_refs.output.bim
    output: config["sldsr"]["lift_hapmap3_snps"]["output"]
    params: hg19_dir = os.path.dirname(HG19_REFS["bim"]),
            hg38_dir = os.path.dirname(HG38_REFS["bim"])
    singularity: config["containers"]["ubuntu"]
    benchmark: "reports/benchmarks/09sldsr.lift_hapmap3_snps.benchmark.txt"
    log:    config["sldsr"]["lift_hapmap3_snps"]["log"]
    shell:
        """
        scripts/ldsr_liftover_hapmap3_snps.sh {output} {input.snps} \
          {params.hg19_dir} {params.hg38_dir} > {log} 2>&1
        """

rule ldsr_ld_scores_hg38:
    input:  annot = get_annot,
            bfile = rules.get_hg38_refs.output.bim,
            snps = rules.lift_hapmap3_snps.output
    output: config["sldsr"]["ldsr_ld_scores_hg38"]["output"]
    params: bfile = strip_suffix(HG38_REFS["bim"], ".bim"),
            ldscores = strip_suffix(config["sldsr"]["ldsr_ld_scores_hg38"]["output"], ".l2.ldscore.gz")
    benchmark: "reports/benchmarks/09sldsr.ldsr_ld_scores_hg38_{cell_type}_{annot_type}_{chr}.benchmark.txt"
    message: "Generating LD scores on hg38 for {wildcards.cell_type}, {wildcards.annot_type}, chr {wildcards.chr}"
    log:    config["sldsr"]["ldsr_ld_scores_hg38"]["log"]
    shell:
        """
        eval "$(/apps/languages/miniforge3/24.3.0-0/bin/conda shell.bash hook)"
        conda activate ldsr
        python ../resources/ldsr/ldsc/ldsc.py --l2 \
          --bfile {params.bfile} --ld-wind-cm 1 \
          --annot {input.annot} \
          --out {params.ldscores} \
          --print-snps {input.snps} > {log} 2>&1
        """

rule ldsr_strat_hg38_bl_v12:
    input:  gwas = config["prep_gwas"]["munge_sumstats"]["output"],
            ldsr = expand(rules.ldsr_ld_scores_hg38.output, chr=range(1, 23), allow_missing=True)
    output: config["sldsr"]["ldsr_strat_hg38_bl_v12"]["output"]
    params: weights = strip_suffix(HG38_REFS["weights"], "{chr}.l2.ldscore.gz"),
            baseline = strip_suffix(HG38_REFS["baseline"], "{chr}.l2.ldscore.gz"),
            frq = strip_suffix(HG38_REFS["frq"], "{chr}.frq"),
            ldscores = strip_suffix(config["sldsr"]["ldsr_ld_scores_hg38"]["output"], "{chr}.l2.ldscore.gz"),
            out_prefix = strip_suffix(config["sldsr"]["ldsr_strat_hg38_bl_v12"]["output"], ".results")
    message: "Running stratified LDSR with {wildcards.cell_type} using hg38 refs, baseline 1.2, {wildcards.annot_type} and {wildcards.gwas} GWAS"
    benchmark: "reports/benchmarks/09sldsr.ldsr_strat_hg38_bl_v12_{cell_type}_{annot_type}_{gwas}.benchmark.txt"
    log:    config["sldsr"]["ldsr_strat_hg38_bl_v12"]["log"]
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
    input:  expand(rules.ldsr_strat_hg38_bl_v12.output, cell_type=config["cell_types"], gwas=config["gwas"], allow_missing=True)
    output: config["sldsr"]["ldsr_strat_summary"]["output"]
    benchmark: "reports/benchmarks/09sldsr.ldsr_strat_summary_{annot_type}.benchmark.txt"
    log:    config["sldsr"]["ldsr_strat_summary"]["log"]
    shell:
        """
        head -1 {input[0]} > {output}
        grep L2_1 {input} | cut -f 6- -d/ >> {output}
        """

rule ldsr_report:
    input:  ldsr = expand(rules.ldsr_strat_summary.output, annot_type=config["annot_types"]),
            rmd_script = "scripts/ldsr_report.Rmd"
    output: config["sldsr"]["ldsr_report"]["html"]
    params: in_dir = config["sldsr"]["ldsr_report"]["in_dir"],
            output_file = config["sldsr"]["ldsr_report"]["output_file"],
            bmark_dir = config["sldsr"]["ldsr_report"]["bmark_dir"]
    singularity: config["containers"]["r_eqtl"]
    message: "Generate SLDSR report"
    benchmark: "reports/benchmarks/09sldsr.ldsr_report.benchmark.txt"
    log:    config["sldsr"]["ldsr_report"]["log"]
    shell:
        """
        Rscript -e "rmarkdown::render('{input.rmd_script}', \
            output_file = '{params.output_file}', \
            params = list(in_dir = c('{params.in_dir}'),
            bmark_dir = '{params.bmark_dir}'))" > {log} 2>&1
        """
