import re


##  Settings  ------------------------------------------------------------------------------
TQ = config["tensorQTL"]

# Pseudotime bins (NPC-to-*) excluded from the tensorQTL pipeline: they are mapped in
# 06PSEUDOTIME and exported into the 05TENSORQTL folders for downstream pipelines.
TQTL_CELL_TYPES = [ct for ct in config["cell_types"] if "-to-" not in ct]

# Per-rule constraint so these rules never claim the pseudotime cell types; global
# constraints would also apply to downstream pipelines using {cell_type}
TQTL_CELL_TYPE_REGEX = "|".join(map(re.escape, TQTL_CELL_TYPES))

# Rules whose benchmarks the report summarises; templates are expanded in the report
TQTL_BENCHMARK_RULES = ["prep_tensorQTL_input", "zip_pblk_cnts", "convert_genotypes", "split_covariates",
                        "tensorqtl_nom", "tensorqtl_perm", "tensorqtl_independent", "tensorqtl_report"]

TQTL_FINAL_NORM = TQ["norm_methods"][0]
TQTL_L1_TYPES   = [ct for ct in TQTL_CELL_TYPES if ct in config["fetal_l1_types"]]
TQTL_L2_TYPES   = [ct for ct in TQTL_CELL_TYPES if ct not in config["fetal_l1_types"]]


##  Targets  -------------------------------------------------------------------------------
# Runs the report needs. Sensitivity (perm): L1 cell types with every norm method,
# L2 with the final norm method only; both over every expPC. Independent: final
# norm method at each cell type's exp_pc_map value.
TQTL_PERM_RUNS = (
    expand(TQ["tensorqtl_perm"]["output"], cell_type = TQTL_L1_TYPES,
           norm_method = TQ["norm_methods"],
           geno_pc = TQ["geno_pcs"], exp_pc = TQ["exp_pcs"]) +
    expand(TQ["tensorqtl_perm"]["output"], cell_type = TQTL_L2_TYPES,
           norm_method = TQTL_FINAL_NORM,
           geno_pc = TQ["geno_pcs"], exp_pc = TQ["exp_pcs"]))
TQTL_INDEP_RUNS = [TQ["tensorqtl_independent"]["output"].format(
                       cell_type = ct, norm_method = TQTL_FINAL_NORM,
                       geno_pc = TQ["geno_pcs"], exp_pc = config["exp_pc_map"][ct])
                   for ct in TQTL_CELL_TYPES]

# Nominal results used downstream (SMR, replication, SuSiE, TWAS): final norm method at
# each cell type's exp_pc_map value
TQTL_NOM_RUNS = [TQ["tensorqtl_nom"]["output"].format(
                     cell_type = ct, norm_method = TQTL_FINAL_NORM,
                     geno_pc = TQ["geno_pcs"], exp_pc = config["exp_pc_map"][ct])
                 for ct in TQTL_CELL_TYPES]


##  Rules  ---------------------------------------------------------------------------------

rule prep_tensorQTL_input:
    input:  cov_file    = config["geno_post_impute"]["calc_genotype_pcs"]["output"],
            sex_file    = config["scanpy"]["metadata_per_sample"],
            gene_lookup = config["gene_lookup"],
            pseudobulk  = config["scanpy"]["pseudobulk_csv"]
    output: cov_out     = TQ["prep_tensorQTL_input"]["cov_out"],
            exp_out     = TQ["prep_tensorQTL_input"]["exp_out"],
            report_rds  = TQ["prep_tensorQTL_input"]["report_rds"]
    params: cell_type   = "{cell_type}",
            norm_method = "{norm_method}"
    singularity: config["containers"]["r_eqtl"]
    resources: threads = 1, mem_mb = 6000, time="5:00:00"
    message: "Prep pseudoblk GeX and covariate matricies for tensorQTL with norm: {wildcards.norm_method}"
    benchmark: TQ["prep_tensorQTL_input"]["benchmark"]
    log:    TQ["prep_tensorQTL_input"]["log"]
    wildcard_constraints: cell_type = TQTL_CELL_TYPE_REGEX
    script: TQ["prep_tensorQTL_input"]["script"]

rule zip_pblk_cnts:
    input:   rules.prep_tensorQTL_input.output.exp_out
    output:  TQ["zip_pblk_cnts"]["output"]
    envmodules: "HTSlib/1.19.1-GCC-13.2.0"   # bgzip and tabix
    message: "bgzip and index pseudoblk counts for norm: {wildcards.norm_method}"
    benchmark: TQ["zip_pblk_cnts"]["benchmark"]
    log:     TQ["zip_pblk_cnts"]["log"]
    wildcard_constraints: cell_type = TQTL_CELL_TYPE_REGEX
    shell:
             """
             bgzip -c {input} > {output} 2> {log}
             tabix -p bed {output} >> {log} 2>&1
             """

rule convert_genotypes:
    input:  config["geno_post_impute"]["exclude_SNPs"]["output"]
    output: TQ["convert_genotypes"]["output"]
    params: TQ["convert_genotypes"]["prefix_out"]
    envmodules: "plink/2.0"
    message: "Convert genotypes to plink2 format for tensorQTL"
    benchmark: TQ["convert_genotypes"]["benchmark"]
    log:    TQ["convert_genotypes"]["log"]
    shell: "plink2 --vcf {input} --make-pgen --out {params} > {log} 2>&1"

rule split_covariates:
    input:  rules.prep_tensorQTL_input.output.cov_out
    output:  TQ["split_covariates"]["output"],
    params:  cell_type = "{cell_type}"
    singularity: config["containers"]["R"]
    resources: threads = 1, mem_mb = 6000, time="5:00:00"
    message: "Divide covariate file to test different PC thresholds in tensorQTL for norm: {wildcards.norm_method}"
    benchmark: TQ["split_covariates"]["benchmark"]
    log:    TQ["split_covariates"]["log"]
    wildcard_constraints: cell_type = TQTL_CELL_TYPE_REGEX
    script: TQ["split_covariates"]["script"]

rule tensorqtl_nom:
    input:  genotypes = rules.convert_genotypes.output,
            counts = rules.zip_pblk_cnts.output,
            covariates = rules.split_covariates.output
    output: TQ["tensorqtl_nom"]["output"]
    params: prefix_in = TQ["tensorqtl_nom"]["prefix_in"],
            prefix_out = TQ["tensorqtl_nom"]["prefix_out"],
            window = TQ["window"]
    singularity: config["containers"]["tensorqtl"]
    resources: threads = 10, mem_mb = 100000, time="5:00:00"
    message: "Run tensorQTL nominal for norm: {wildcards.norm_method}"
    benchmark: TQ["tensorqtl_nom"]["benchmark"]
    log:    TQ["tensorqtl_nom"]["log"]
    wildcard_constraints: cell_type = TQTL_CELL_TYPE_REGEX
    shell:
            """
            python3 -m tensorqtl {params.prefix_in} {input.counts} {params.prefix_out} \
               --covariates {input.covariates} \
               --window {params.window} \
               --mode cis_nominal > {log} 2>&1
            """

rule tensorqtl_perm:
    input:  genotypes = rules.convert_genotypes.output,
            counts = rules.zip_pblk_cnts.output,
            covariates = rules.split_covariates.output
    output: TQ["tensorqtl_perm"]["output"]
    params: prefix_in = TQ["tensorqtl_perm"]["prefix_in"],
            prefix_out = TQ["tensorqtl_perm"]["prefix_out"],
            window = TQ["window"]
    singularity: config["containers"]["tensorqtl"]
    resources: threads = 10, mem_mb = 20000, time="5:00:00"
    message: "Run tensorQTL permutation for norm: {wildcards.norm_method}"
    benchmark: TQ["tensorqtl_perm"]["benchmark"]
    log:    TQ["tensorqtl_perm"]["log"]
    wildcard_constraints: cell_type = TQTL_CELL_TYPE_REGEX
    shell:
            """
            python3 -m tensorqtl {params.prefix_in} {input.counts} {params.prefix_out} \
               --covariates {input.covariates} \
               --window {params.window} \
               --mode cis > {log} 2>&1
            """

rule tensorqtl_independent:
    input:  genotypes = rules.convert_genotypes.output,
            counts = rules.zip_pblk_cnts.output,
            covariates = rules.split_covariates.output,
            cis_perm = rules.tensorqtl_perm.output
    output: TQ["tensorqtl_independent"]["output"]
    params: prefix_in = TQ["tensorqtl_independent"]["prefix_in"],
            prefix_out = TQ["tensorqtl_independent"]["prefix_out"],
            window = TQ["window"]
    singularity: config["containers"]["tensorqtl"]
    resources: threads = 10, mem_mb = 100000, time="5:00:00"
    message: "Run tensorQTL independent analysis for norm: {wildcards.norm_method}"
    benchmark: TQ["tensorqtl_independent"]["benchmark"]
    log:    TQ["tensorqtl_independent"]["log"]
    wildcard_constraints: cell_type = TQTL_CELL_TYPE_REGEX
    shell:
            """
            python3 -m tensorqtl {params.prefix_in} {input.counts} {params.prefix_out} \
              --covariates {input.covariates} \
              --window {params.window} \
              --mode cis_independent \
              --cis_output {input.cis_perm} > {log} 2>&1
            """

rule tensorqtl_report:
    # Note diff paths for output and out_file; Rmarkdown needs outfile to be relative to Rmd file
    input:  perm  = TQTL_PERM_RUNS,
            indep = TQTL_INDEP_RUNS,
            cell_counts = config["scanpy"]["l1_cell_counts"],
            rmd_script = TQ["tensorqtl_report"]["script"]
    output: TQ["tensorqtl_report"]["output"]
    params: in_dir       = TQ["tensorqtl_report"]["in_dir"],
            ind_dir      = TQ["tensorqtl_report"]["ind_dir"],
            cell_counts  = "../" + config["scanpy"]["l1_cell_counts"],   # the Rmd knits in scripts/
            output_file  = TQ["tensorqtl_report"]["output_file"],
            l1_types     = ",".join(config["fetal_l1_types"]),
            l2_types     = ",".join(TQTL_L2_TYPES),
            exp_pcs      = ",".join(map(str, TQ["exp_pcs"])),
            geno_pc      = TQ["geno_pcs"],
            norm_methods = ",".join(TQ["norm_methods"]),
            window       = TQ["window"],
            exp_pc_map   = ",".join(f"{ct}:{config['exp_pc_map'][ct]}" for ct in TQTL_CELL_TYPES),
            eqtl_fdr     = config["eqtl_fdr"],
            benchmarks   = lambda w: ",".join(f"{r}=../{TQ[r]['benchmark']}" for r in TQTL_BENCHMARK_RULES)
    singularity: config["containers"]["r_eqtl"]
    message: "Generate tensorQTL report"
    benchmark: TQ["tensorqtl_report"]["benchmark"]
    log:     TQ["tensorqtl_report"]["log"]
    shell:
        """
        Rscript -e "rmarkdown::render('{input.rmd_script}', \
            output_file = '{params.output_file}', \
            params = list(in_dir = '{params.in_dir}', ind_dir = '{params.ind_dir}', \
                          cell_counts = '{params.cell_counts}', \
                          l1_types = '{params.l1_types}', \
                          l2_types = '{params.l2_types}', exp_pcs = '{params.exp_pcs}', \
                          geno_pc = '{params.geno_pc}', norm_methods = '{params.norm_methods}', \
                          window = '{params.window}', exp_pc_map = '{params.exp_pc_map}', \
                          eqtl_fdr = '{params.eqtl_fdr}', benchmarks = '{params.benchmarks}'))" > {log} 2>&1
        """


##  Final outputs, collected by rule all in the Snakefile  ---------------------------------
TARGETS = globals().setdefault("TARGETS", {})
TARGETS["05tensorQTL"] = TQTL_NOM_RUNS + [TQ["tensorqtl_report"]["output"]]
