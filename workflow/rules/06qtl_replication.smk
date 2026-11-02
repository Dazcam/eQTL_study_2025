configfile: "../config/config.yaml"

localrules: dwnld_ziffra, snp_lookup

QR          = config["qtl_rep"]
NORM_METHOD = config["tensorQTL"]["norm_methods"][0]
GENO_PC     = config["tensorQTL"]["geno_pcs"]

# Pseudotime bins (NPC-to-*) excluded from the replication pipeline
REP_CELL_TYPES = [ct for ct in config["cell_types"] if "-to-" not in ct]
REP_FETAL_L1   = QR["fetal_l1_types"]

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


# Targets  ---------------------------------------------------------------------------

# Ziffra OCR enrichment: all L1 + L2 cell types
ziffra_inputs = [
    QR["atac_enrich"]["out_file"].format(cell_type=ct, norm_method=NORM_METHOD,
                                         geno_pc=GENO_PC, exp_pc=config["exp_pc_map"][ct])
    for ct in REP_CELL_TYPES
]

# Internal pi1: fetal L1 x fetal L1
internal_inputs = [
    QR["pi1_enrich_internal"]["pi1"].format(cell_type=ct, ref_cell_type=ref_ct,
                                            norm_method=NORM_METHOD, geno_pc=GENO_PC,
                                            exp_pc=config["exp_pc_map"][ct],
                                            ref_exp_pc=config["exp_pc_map"][ref_ct])
    for ct in REP_FETAL_L1 for ref_ct in REP_FETAL_L1
]

# Jang pi1: fetal L1 x Jang L1
jang_pi1_inputs = [
    QR["pi1_enrich_jang"]["pi1"].format(cell_type=ct, jang_cell_type=jct)
    for ct in REP_FETAL_L1 for jct in config["jang_l1_types"]
]


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

# Freeze 2 outputs of this rule were made from the 10SMR nominal export, which
# used expPC 40 for every cell type; qtl_all now points at the final expPC.
rule pi1_enrich_jang:
    input:  public_all = QR["pi1_enrich_jang"]["public_all"],
            public_top = QR["pi1_enrich_jang"]["public_top"],
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
            jang_l1_files     = expand(QR["pi1_enrich_jang"]["public_all"], jang_cell_type=BETA_JANG_TYPES),
            jang_l1_top_files = expand(QR["pi1_enrich_jang"]["public_top"], jang_cell_type=BETA_JANG_TYPES),
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

# Report params are rebuilt with the report (stage 6); output_file is relative to the Rmd
rule pi1_enrichment_report:
    input:  ziffra = ziffra_inputs,
            internal = internal_inputs,
            jang_pi1 = jang_pi1_inputs,
            beta_correlation = rules.beta_correlation.output[0],
            rmd_script = config["qtl_rep"]["report"]["script"]
    output: config["qtl_rep"]["report"]["html"]
    params: ziffra_dir = "../../results/06QTL-REPLICATION/atac_enrich/",
            ziffra_indep_dir = "../../results/06QTL-REPLICATION/atac_enrich_indep/",
            obrien_dir = "../../results/06QTL-REPLICATION/obrien_bulk/",
            wen_dir = "../../results/06QTL-REPLICATION/wen_bulk/",
            bryois_dir = "../../results/06QTL-REPLICATION/bryois/",
            fugita_dir = "../../results/06QTL-REPLICATION/fugita/",
            internal_dir = "../../results/06QTL-REPLICATION/internal/",
            beta_dir = "../../results/06QTL-REPLICATION/beta_cor/",
            tbl_dir = "../../results/13MANUSCRIPT_PLOTS_TABLES/tables/",
            bmark_dir = "../reports/benchmarks/",
            output_file = config["qtl_rep"]["report"]["output_file"]
    singularity: config["containers"]["r_eqtl"]
    message: "Generate pi1 enrichments report report"
    resources: time="1:00:00"
    benchmark: config["qtl_rep"]["report"]["benchmark"]
    log:     config["qtl_rep"]["report"]["log"]
    shell:
        """
        Rscript -e "rmarkdown::render('{input.rmd_script}', \
            output_file = '{params.output_file}', \
            params = list(ziffra_dir = '{params.ziffra_dir}',
            ziffra_indep_dir = '{params.ziffra_indep_dir}',
            obrien_dir = '{params.obrien_dir}',
            wen_dir = '{params.wen_dir}',
            bryois_dir = '{params.bryois_dir}',
            internal_dir = '{params.internal_dir}',
            fugita_dir = '{params.fugita_dir}',
            beta_dir = '{params.beta_dir}',
            tbl_dir = '{params.tbl_dir}',
            bmark_dir = '{params.bmark_dir}'))" > {log} 2>&1
        """

