configfile: "../config/config.yaml"

RB = config["replication_bulk"]


def perm_files(cell_types):
    """tensorQTL permutation output per cell type, at the final expPC (exp_pc_map)."""
    return [config["tensorQTL"]["tensorqtl_perm"]["output"].format(
                cell_type=ct, norm_method=config["tensorQTL"]["norm_methods"][0],
                geno_pc=config["tensorQTL"]["geno_pcs"], exp_pc=config["exp_pc_map"][ct])
            for ct in cell_types]


# Benchmark file per rule, shown in the report
REPL_BULK_BENCHMARKS = [(rule, RB[key]["benchmark"]) for rule, key in [
    ("extract_unique_egenes",   "extract_unique_egenes"),
    ("overlap_obrien",          "overlap_obrien"),
    ("isoform_eqtl_comparison", "isoform_eqtl"),
    ("repl_bulk_report",        "report")]]


rule extract_unique_egenes:
    input:  perm_files = perm_files(config["cell_types"])
    output: RB["extract_unique_egenes"]["out_file"]
    params: cell_types  = config["cell_types"],
            l1_types    = config["fetal_l1_types"],
            eqtl_fdr    = config["eqtl_fdr"],
            geno_pc     = config["tensorQTL"]["geno_pcs"],
            norm_method = config["tensorQTL"]["norm_methods"][0]
    singularity: config["containers"]["r_eqtl"]
    message:   "Extract unique FDR-sig eGenes across all L1+L2+pseudotime cell types"
    benchmark: RB["extract_unique_egenes"]["benchmark"]
    log:       RB["extract_unique_egenes"]["log"]
    script:    RB["extract_unique_egenes"]["script"]

rule overlap_obrien:
    input:  fetal_egenes = rules.extract_unique_egenes.output[0],
            obrien_supp  = RB["overlap_obrien"]["obrien_supp"]
    output: RB["overlap_obrien"]["out_file"]
    singularity: config["containers"]["r_eqtl"]
    message:   "Overlap unique fetal eGenes against O'Brien 2018 bulk gene- and transcript-level significant eQTL"
    benchmark: RB["overlap_obrien"]["benchmark"]
    log:       RB["overlap_obrien"]["log"]
    script:    RB["overlap_obrien"]["script"]

rule isoform_eqtl_comparison:
    input:  fetal_egenes = rules.extract_unique_egenes.output[0],
            gtf          = config["parse"]["get_refs"]["gtf_out"],
            bulk_eqtl    = RB["isoform_eqtl"]["bulk_eqtl"],
            gene_lookup  = config["gene_lookup"]
    output: RB["isoform_eqtl"]["summary"]
    params: nominal_threshold = RB["isoform_eqtl"]["nominal_threshold"]
    conda:  config["scanpy"]["env"]
    resources: threads = 8, mem_mb = 48000, time = "3:00:00"
    threads: 8
    message:   "Comparing snRNA-seq eGenes against bulk transcript-level eQTL (O'Brien et al.)"
    benchmark: RB["isoform_eqtl"]["benchmark"]
    log:       RB["isoform_eqtl"]["log"]
    script:    RB["isoform_eqtl"]["script"]

rule repl_bulk_report:
    input:  egenes_per_celltype = rules.extract_unique_egenes.output[0],
            overlap             = rules.overlap_obrien.output[0],
            isoform_summary     = rules.isoform_eqtl_comparison.output[0],
            rmd_script          = RB["report"]["script"]
    output: RB["report"]["html"]
    params: egenes_per_celltype = lambda wc, input: "../" + input.egenes_per_celltype,
            overlap             = lambda wc, input: "../" + input.overlap,
            isoform_summary     = lambda wc, input: "../" + input.isoform_summary,
            eqtl_fdr            = config["eqtl_fdr"],
            nominal_threshold   = RB["isoform_eqtl"]["nominal_threshold"],
            benchmarks          = ",".join(f"{rule}=../{path}" for rule, path in REPL_BULK_BENCHMARKS),
            output_file         = RB["report"]["output_file"]
    singularity: config["containers"]["r_eqtl"]
    message:   "Generate replication_bulk report"
    benchmark: RB["report"]["benchmark"]
    log:       RB["report"]["log"]
    shell:
        """
        Rscript -e "rmarkdown::render('{input.rmd_script}', \
            output_file = '{params.output_file}', \
            params = list(egenes_per_celltype = '{params.egenes_per_celltype}', \
                          overlap = '{params.overlap}', \
                          isoform_summary = '{params.isoform_summary}', \
                          eqtl_fdr = '{params.eqtl_fdr}', \
                          nominal_threshold = '{params.nominal_threshold}', \
                          benchmarks = '{params.benchmarks}'))" > {log} 2>&1
        """
