import re

configfile: "../config/config.yaml"

localrules: dwnld_ziffra, snp_lookup

QR          = config["replication_sc"]
NORM_METHOD = config["tensorQTL"]["norm_methods"][0]
GENO_PC     = config["tensorQTL"]["geno_pcs"]

# Pseudotime bins (NPC-to-*) excluded from the replication pipeline
REP_CELL_TYPES = [ct for ct in config["cell_types"] if "-to-" not in ct]
REP_FETAL_L1   = config["fetal_l1_types"]

# Jang et al. 2026 association files (downloaded by setup get_jang_singlebrain)
JANG_FULL = config["setup"]["get_jang_singlebrain"]["full_assoc"]
JANG_TOP  = config["setup"]["get_jang_singlebrain"]["top_assoc"]

# beta_correlation: fetal-adult pairs, and the unique cell types on each side
# (order matters: the script names its input files in this order)
JANG_COMPARISONS = QR["jang_comparisons"]
BETA_FETAL_TYPES = list(dict.fromkeys(c["fetal"] for c in JANG_COMPARISONS))
BETA_JANG_TYPES  = list(dict.fromkeys(c["jang"] for c in JANG_COMPARISONS))


def perm_file(cell_type, norm_method=NORM_METHOD, geno_pc=GENO_PC, exp_pc=None):
    """tensorQTL permutation output; expPC defaults to the final value in exp_pc_map."""
    exp_pc = config["exp_pc_map"][cell_type] if exp_pc is None else exp_pc
    return config["tensorQTL"]["tensorqtl_perm"]["output"].format(
        cell_type=cell_type, norm_method=norm_method, geno_pc=geno_pc, exp_pc=exp_pc)


def nom_file(cell_type, norm_method=NORM_METHOD, geno_pc=GENO_PC, exp_pc=None):
    """Concatenated tensorQTL nominal output (cat_nom_qtl); expPC as perm_file()."""
    exp_pc = config["exp_pc_map"][cell_type] if exp_pc is None else exp_pc
    return QR["cat_nom_qtl"]["output"].format(
        cell_type=cell_type, norm_method=norm_method, geno_pc=geno_pc, exp_pc=exp_pc)


def nom_parquets(wc):
    """Per-chromosome tensorQTL nominal parquets (chr 1-22) for one cat_nom_qtl job."""
    return expand(QR["cat_nom_qtl"]["parquet"], cell_type=wc.cell_type,
                  norm_method=wc.norm_method, geno_pc=wc.geno_pc,
                  exp_pc=wc.exp_pc, chr=range(1, 23))


def l1_parent(cell_type):
    """L1 parent of a subcluster (GABA-1 -> GABA); L1 types return themselves."""
    return re.sub(r"-[0-9]+$", "", cell_type)


def ziffra_file(ct):
    return QR["atac_enrich"]["out_file"].format(cell_type=ct, norm_method=NORM_METHOD,
                                                geno_pc=GENO_PC, exp_pc=config["exp_pc_map"][ct])


def internal_file(kind, ct, ref_ct):
    """kind: 'enrich' or 'pi1'."""
    return QR["pi1_enrich_internal"][kind].format(cell_type=ct, ref_cell_type=ref_ct,
                                                  norm_method=NORM_METHOD, geno_pc=GENO_PC,
                                                  exp_pc=config["exp_pc_map"][ct],
                                                  ref_exp_pc=config["exp_pc_map"][ref_ct])


def jang_file(kind, ct, jct):
    """kind: 'enrich' or 'pi1'."""
    return QR["pi1_enrich_jang"][kind].format(cell_type=ct, jang_cell_type=jct)


# Targets  ---------------------------------------------------------------------------

# Ziffra OCR enrichment: all L1 + L2 cell types, each L1 followed by its subclusters
REP_L1_CONFIG_ORDER = [ct for ct in REP_CELL_TYPES if l1_parent(ct) == ct]
ZIFFRA_CELL_TYPES   = [ct for l1 in REP_L1_CONFIG_ORDER
                       for ct in REP_CELL_TYPES if l1_parent(ct) == l1]
ziffra_inputs = [ziffra_file(ct) for ct in ZIFFRA_CELL_TYPES]

# Internal pi1: fetal L1 x fetal L1
internal_inputs = [internal_file(kind, ct, ref_ct) for kind in ["enrich", "pi1"]
                   for ct in REP_FETAL_L1 for ref_ct in REP_FETAL_L1]

# Jang pi1: fetal L1 x Jang L1
JANG_L1_ORDER = config["jang_l1_order"]
assert sorted(JANG_L1_ORDER) == sorted(config["jang_l1_types"]), \
    "jang_l1_order must list the same types as jang_l1_types"
jang_pi1_inputs = [jang_file(kind, ct, jct) for kind in ["enrich", "pi1"]
                   for ct in REP_FETAL_L1 for jct in config["jang_l1_types"]]

# Benchmark files of the current runs, by rule (read by the report if present)
REP_BENCHMARKS = (
    [("dwnld_ziffra", QR["dwnld_ziffra"]["benchmark"])]
    + [(rule, QR[rule]["benchmark"].format(cell_type=ct, norm_method=NORM_METHOD,
                                           geno_pc=GENO_PC, exp_pc=config["exp_pc_map"][ct]))
       for rule in ["snp_lookup", "atac_enrich"] for ct in ZIFFRA_CELL_TYPES]
    + [("ziffra_summary", QR["ziffra_summary"]["benchmark"])]
    + [("cat_nom_qtl", QR["cat_nom_qtl"]["benchmark"].format(cell_type=ct, norm_method=NORM_METHOD,
                                                            geno_pc=GENO_PC, exp_pc=config["exp_pc_map"][ct]))
       for ct in REP_FETAL_L1]
    + [("pi1_enrich_internal", QR["pi1_enrich_internal"]["benchmark"].format(
            cell_type=ct, ref_cell_type=ref_ct, norm_method=NORM_METHOD, geno_pc=GENO_PC,
            exp_pc=config["exp_pc_map"][ct], ref_exp_pc=config["exp_pc_map"][ref_ct]))
       for ct in REP_FETAL_L1 for ref_ct in REP_FETAL_L1]
    + [("pi1_enrich_jang", QR["pi1_enrich_jang"]["benchmark"].format(cell_type=ct, jang_cell_type=jct))
       for ct in REP_FETAL_L1 for jct in config["jang_l1_types"]]
    + [("beta_correlation", QR["beta_correlation"]["benchmark"]),
       ("repl_sc_report", QR["report"]["benchmark"])]
)


def rmd_path(path):
    """Paths passed into the Rmd are relative to workflow/scripts/, where it knits."""
    return "../" + path


# Rules  -----------------------------------------------------------------------------

rule dwnld_ziffra:
    output: QR["dwnld_ziffra"]["output"]
    params: web_link = QR["dwnld_ziffra"]["web_link"]
    message: "Download snATAC-seq peaks file from Ziffra 2021, PMID:34616060"
    benchmark: QR["dwnld_ziffra"]["benchmark"]
    log:    QR["dwnld_ziffra"]["log"]
    shell:  """
            wget -O {output} {params.web_link} &>> {log}
            """

rule snp_lookup:
    input:  qtl_perm = config["tensorQTL"]["tensorqtl_perm"]["output"]
    output: QR["snp_lookup"]["out_file"]
    message: "Generate snp lookup for Ziffra ATAC-seq peaks; Net access req. - run local"
    singularity: config["containers"]["r_eqtl"]
    benchmark: QR["snp_lookup"]["benchmark"]
    log:    QR["snp_lookup"]["log"]
    script: QR["snp_lookup"]["script"]

rule atac_enrich:
    input:  qtl_perm = config["tensorQTL"]["tensorqtl_perm"]["output"],
            peaks    = rules.dwnld_ziffra.output,
            snp_file = rules.snp_lookup.output
    output: QR["atac_enrich"]["out_file"]
    params: peak_sets = QR["ziffra_peak_sets"],
            seed      = QR["atac_enrich"]["seed"]
    singularity: config["containers"]["r_eqtl"]
    message: "Test for sig. eQTL enrichments in Ziffra ATAC-seq peaks"
    benchmark: QR["atac_enrich"]["benchmark"]
    log:    QR["atac_enrich"]["log"]
    script: QR["atac_enrich"]["script"]

rule ziffra_summary:
    input:  enrich = ziffra_inputs
    output: QR["ziffra_summary"]["out_file"]
    params: cell_types = ZIFFRA_CELL_TYPES
    singularity: config["containers"]["r_eqtl"]
    message: "Combine Ziffra ATAC-seq peak enrichments across cell types (Fig 2E)"
    benchmark: QR["ziffra_summary"]["benchmark"]
    log:    QR["ziffra_summary"]["log"]
    script: QR["ziffra_summary"]["script"]

rule cat_nom_qtl:
    input:  parquets = nom_parquets
    output: QR["cat_nom_qtl"]["output"]
    params: script = QR["cat_nom_qtl"]["script"]
    resources: threads = 1, mem_mb = 10000, time="1:00:00"
    message: "Combine per-chr parquet files into gzip file for pi1 enrichment"
    benchmark: QR["cat_nom_qtl"]["benchmark"]
    log:    QR["cat_nom_qtl"]["log"]
    shell:  "{params.script} {output} {input.parquets} > {log} 2>&1"

rule pi1_enrich_internal:
    input:  public_all = lambda wc: nom_file(wc.ref_cell_type, wc.norm_method, wc.geno_pc, wc.ref_exp_pc),
            public_top = lambda wc: perm_file(wc.ref_cell_type, wc.norm_method, wc.geno_pc, wc.ref_exp_pc),
            qtl_all    = lambda wc: nom_file(wc.cell_type, wc.norm_method, wc.geno_pc, wc.exp_pc),
            qtl_top    = lambda wc: perm_file(wc.cell_type, wc.norm_method, wc.geno_pc, wc.exp_pc)
    output: enrich = QR["pi1_enrich_internal"]["enrich"],
            pi1    = QR["pi1_enrich_internal"]["pi1"]
    resources: threads = 4, mem_mb = 20000, time="1:00:00"
    message: "Calc internal pi1 enrichments between {wildcards.cell_type} and sn-eQTL {wildcards.ref_cell_type}"
    singularity: config["containers"]["r_eqtl"]
    benchmark: QR["pi1_enrich_internal"]["benchmark"]
    log:    QR["pi1_enrich_internal"]["log"]
    script: QR["pi1_enrich_internal"]["script"]

rule pi1_enrich_jang:
    input:  public_all = JANG_FULL,
            public_top = JANG_TOP,
            qtl_all    = lambda wc: nom_file(wc.cell_type),
            qtl_top    = lambda wc: perm_file(wc.cell_type)
    output: enrich = QR["pi1_enrich_jang"]["enrich"],
            pi1    = QR["pi1_enrich_jang"]["pi1"]
    resources: threads = 4, mem_mb = 60000, time="1:00:00"
    message: "Calc pi1 enrichments between fetal L1 sn-eQTL and Jang et al. 2026 L1 eQTL"
    singularity: config["containers"]["r_eqtl"]
    benchmark: QR["pi1_enrich_jang"]["benchmark"]
    log:    QR["pi1_enrich_jang"]["log"]
    script: QR["pi1_enrich_jang"]["script"]

rule beta_correlation:
    input:  fetal_perm        = [perm_file(ct) for ct in BETA_FETAL_TYPES],
            jang_l1_files     = expand(JANG_FULL, jang_cell_type=BETA_JANG_TYPES),
            jang_l1_top_files = expand(JANG_TOP, jang_cell_type=BETA_JANG_TYPES),
            pvar_file         = config["tensorQTL"]["convert_genotypes"]["pvar"]
    output: QR["beta_correlation"]["out_file"]
    params: fetal_cell = [c["fetal"] for c in JANG_COMPARISONS],
            jang_cell  = [c["jang"] for c in JANG_COMPARISONS]
    singularity: config["containers"]["r_eqtl"]
    resources: threads = 4, mem_mb = 64000, time = "04:00:00"
    threads: 4
    message: "Compute fetal-adult eQTL beta correlations and flag opposite-effect pairs with allele-flip annotations"
    benchmark: QR["beta_correlation"]["benchmark"]
    log:    QR["beta_correlation"]["log"]
    script: QR["beta_correlation"]["script"]

rule repl_sc_report:
    input:  ziffra_tbl       = rules.ziffra_summary.output[0],
            internal         = internal_inputs,
            jang_pi1         = jang_pi1_inputs,
            beta_correlation = rules.beta_correlation.output[0],
            gene_lookup      = config["gene_lookup"],
            rmd_script       = QR["report"]["script"]
    output: QR["report"]["html"]
    params: ziffra_tbl       = rmd_path(QR["ziffra_summary"]["out_file"]),
            internal_enrich  = ",".join(rmd_path(internal_file("enrich", ct, r)) for ct in REP_FETAL_L1 for r in REP_FETAL_L1),
            internal_pi1     = ",".join(rmd_path(internal_file("pi1", ct, r)) for ct in REP_FETAL_L1 for r in REP_FETAL_L1),
            jang_enrich      = ",".join(rmd_path(jang_file("enrich", ct, j)) for ct in REP_FETAL_L1 for j in JANG_L1_ORDER),
            jang_pi1         = ",".join(rmd_path(jang_file("pi1", ct, j)) for ct in REP_FETAL_L1 for j in JANG_L1_ORDER),
            beta_correlation = rmd_path(QR["beta_correlation"]["out_file"]),
            gene_lookup      = rmd_path(config["gene_lookup"]),
            cell_types       = ",".join(REP_CELL_TYPES),
            fetal_l1         = ",".join(REP_FETAL_L1),
            jang_l1          = ",".join(JANG_L1_ORDER),
            exp_pcs          = ",".join(f"{ct}={config['exp_pc_map'][ct]}" for ct in ZIFFRA_CELL_TYPES),
            geno_pc          = GENO_PC,
            norm_method      = NORM_METHOD,
            seed             = QR["atac_enrich"]["seed"],
            benchmarks       = ",".join(f"{rule}={rmd_path(path)}" for rule, path in REP_BENCHMARKS),
            output_file      = QR["report"]["output_file"],
            render           = QR["report"]["render"]
    singularity: config["containers"]["r_eqtl"]
    message: "Generate eQTL replication report"
    resources: time="1:00:00"
    benchmark: QR["report"]["benchmark"]
    log:    QR["report"]["log"]
    shell:
        """
        Rscript {params.render} {input.rmd_script} {params.output_file} \
            ziffra_tbl={params.ziffra_tbl} \
            internal_enrich={params.internal_enrich} \
            internal_pi1={params.internal_pi1} \
            jang_enrich={params.jang_enrich} \
            jang_pi1={params.jang_pi1} \
            beta_correlation={params.beta_correlation} \
            gene_lookup={params.gene_lookup} \
            cell_types={params.cell_types} \
            fetal_l1={params.fetal_l1} \
            jang_l1={params.jang_l1} \
            exp_pcs={params.exp_pcs} \
            geno_pc={params.geno_pc} \
            norm_method={params.norm_method} \
            seed={params.seed} \
            benchmarks={params.benchmarks} > {log} 2>&1
        """
