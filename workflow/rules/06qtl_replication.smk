configfile: "../config/config.yaml"

localrules: dwnld_ziffra, snp_lookup

# Pseudotime bins (NPC-to-*) excluded from the replication pipeline
REP_CELL_TYPES = [ct for ct in config["cell_types"] if "-to-" not in ct]
REP_FETAL_L1   = ["Glu-UL", "Glu-DL", "GABA", "NPC", "MG", "OPC", "Endo-Peri"]

CELL_EXP_PC_TUPLES = [(ct, config["exp_pc_map"][ct]) for ct in REP_CELL_TYPES]
L1_EXP_PC_TUPLES   = [(ct, config["exp_pc_map"][ct]) for ct in REP_FETAL_L1]

# Ziffra OCR enrichment: all L1 + L2 cell types
ziffra_inputs = []
for ct, pc in CELL_EXP_PC_TUPLES:
    ziffra_inputs.extend(
        expand(
            config["qtl_rep"]["atac_enrich"]["out_file"],
            cell_type=[ct],
            exp_pc=[pc],
            geno_pc=config["tensorQTL"]["geno_pcs"],
            norm_method=config["tensorQTL"]["norm_methods"][0]
        )
    )

# Internal pi1: fetal L1 x fetal L1
internal_inputs = []
for ct, ct_pc in L1_EXP_PC_TUPLES:
    for ref_ct, ref_pc in L1_EXP_PC_TUPLES:
        internal_inputs.extend(
            expand(
                config["qtl_rep"]["pi1_enrich_internal"]["pi1"],
                cell_type=[ct],
                exp_pc=[ct_pc],
                ref_cell_type=[ref_ct],
                ref_exp_pc=[ref_pc],
                geno_pc=config["tensorQTL"]["geno_pcs"],
                norm_method=config["tensorQTL"]["norm_methods"][0],
            )
        )

# Jang pi1: fetal L1 x Jang L1
jang_pi1_inputs = []
for ct in REP_FETAL_L1:
    jang_pi1_inputs.extend(
        expand(
            config["dev_specificity"]["pi1_enrich_jang"]["pi1"],
            cell_type=[ct],
            jang_cell_type=config["jang_l1_types"],
            exp_pc=[config["exp_pc_map"][ct]],
            geno_pc=config["tensorQTL"]["geno_pcs"],
            norm_method=config["tensorQTL"]["norm_methods"][0]
        )
    )

def get_jang_beta_corr_full_inputs():
    jang_dir = config["dev_specificity"]["classify_specificity"]["jang_dir"]
    return [f"{jang_dir}{ct}_eqtl_full_assoc.tsv.gz" for ct in ["Ext", "IN", "MG", "OPC", "End"]]

def get_jang_beta_corr_top_inputs():
    jang_dir = config["dev_specificity"]["classify_specificity"]["jang_dir"]
    return [f"{jang_dir}{ct}_eqtl_top_assoc.tsv.gz" for ct in ["Ext", "IN", "MG", "OPC", "End"]]


rule dwnld_ziffra:
    output: config["qtl_rep"]["dwnld_ziffra"]["output"]
    params: web_link = config["qtl_rep"]["dwnld_ziffra"]["web_link"],
    message: "Download snATAC-seq peaks file from Ziffra 2021, PMID:34616060"
    benchmark: "reports/benchmarks/06qtl_replication.dwnld_ziffra.txt"
    log:    config["qtl_rep"]["dwnld_ziffra"]["log"]
    shell:  """
            wget -O {output} {params.web_link} &>> {log}
            """

rule snp_lookup:
    input:  qtl_perm = config["qtl_rep"]["snp_lookup"]["qtl_perm"]
    output: config["qtl_rep"]["snp_lookup"]["out_file"]
    message: "Generate snp lookup for Ziffra ATAC-seq peaks; Net access req. - run local"
    singularity: config["containers"]["r_eqtl"]
    benchmark: "reports/benchmarks/qtl_replication.snp_lookup_{cell_type}_{norm_method}_genPC_{geno_pc}_expPC_{exp_pc}.txt"
    log:    config["qtl_rep"]["snp_lookup"]["log"]
    script: "../scripts/repl_sc_snp_lookup.R"

rule atac_enrich:
    input:  qtl_perm = config["qtl_rep"]["snp_lookup"]["qtl_perm"],
            peaks = rules.dwnld_ziffra.output,
            snp_file = rules.snp_lookup.output
    output: config["qtl_rep"]["atac_enrich"]["out_file"]
    params: peak_dir = config["qtl_rep"]["atac_enrich"]["peak_dir"]
    singularity: config["containers"]["r_eqtl"]
    message: "Test for sig. eQTL enrichments in Ziffra ATAC-seq peaks"
    benchmark: "reports/benchmarks/qtl_replication.atac_enrich_{cell_type}_{norm_method}_genPC_{geno_pc}_expPC_{exp_pc}.txt"
    log:    config["qtl_rep"]["atac_enrich"]["log"]
    script: "../scripts/repl_sc_atac_enrichment.R"

rule cat_nom_qtl:
    input:  parquets = lambda wc: expand(config["qtl_rep"]["cat_nom_qtl"]["parquet"], cell_type=wc.cell_type, norm_method=wc.norm_method, geno_pc=wc.geno_pc, exp_pc=wc.exp_pc, chr=range(1, 23), allow_missing=True)
    output:  config["qtl_rep"]["cat_nom_qtl"]["output"]
    resources: threads = 1, mem_mb = 10000, time="1:00:00"
    message: "Combine per-chr parquet files into gzip file for pi1 enrichment"
    benchmark: "reports/benchmarks/06qtl_replication.cat_nom_qtl_{cell_type}_{norm_method}_genPC_{geno_pc}_expPC_{exp_pc}.txt"
    log:     config["qtl_rep"]["cat_nom_qtl"]["log"]
    shell:   "scripts/repl_sc_cat_nom_qtl.py {output} {input.parquets} > {log} 2>&1"

rule pi1_enrich_internal:
    input:
        public_all = "../results/06QTL-REPLICATION/cat_nom_qtl/{ref_cell_type}_{norm_method}_genPC_{geno_pc}_expPC_{ref_exp_pc}_nom.cis_qtl_pairs.tsv.gz",
        public_top = "../results/05TENSORQTL/tensorqtl_perm/{ref_cell_type}_{norm_method}_genPC_{geno_pc}_expPC_{ref_exp_pc}/{ref_cell_type}_{norm_method}_perm.cis_qtl.txt.gz",
        qtl_all    = "../results/06QTL-REPLICATION/cat_nom_qtl/{cell_type}_{norm_method}_genPC_{geno_pc}_expPC_{exp_pc}_nom.cis_qtl_pairs.tsv.gz",
        qtl_top    = "../results/05TENSORQTL/tensorqtl_perm/{cell_type}_{norm_method}_genPC_{geno_pc}_expPC_{exp_pc}/{cell_type}_{norm_method}_perm.cis_qtl.txt.gz"
    output: enrich = config["qtl_rep"]["pi1_enrich_internal"]["enrich"],
            pi1    = config["qtl_rep"]["pi1_enrich_internal"]["pi1"]
    resources: threads = 4, mem_mb = 20000, time="1:00:00"
    message: "Calc internal pi1 enrichments between  {wildcards.cell_type} and sn-eQTL {wildcards.ref_cell_type}"
    singularity: config["containers"]["r_eqtl"]
    benchmark: "reports/benchmarks/06qtl_replication.pi1_enrich_internal_{cell_type}_vs_{ref_cell_type}_{norm_method}_genPC_{geno_pc}_expPC_{exp_pc}_expPCref_{ref_exp_pc}.txt"
    log:    config["qtl_rep"]["pi1_enrich_internal"]["log"]
    script: "../scripts/repl_sc_pi1_enrichment.R"

# Freeze 2 outputs of this rule were made from the 10SMR nominal export, which
# used expPC 40 for every cell type; qtl_all now points at the final expPC.
rule pi1_enrich_jang:
    input: public_all = config["dev_specificity"]["pi1_enrich_jang"]["public_all"],
           public_top = config["dev_specificity"]["pi1_enrich_jang"]["public_top"],
           qtl_all    = lambda wildcards: config["qtl_rep"]["cat_nom_qtl"]["output"].format(
                            cell_type   = wildcards.cell_type,
                            norm_method = config["tensorQTL"]["norm_methods"][0],
                            geno_pc     = config["tensorQTL"]["geno_pcs"],
                            exp_pc      = config["exp_pc_map"][wildcards.cell_type]
                        ),
           qtl_top    = lambda wildcards: config["tensorQTL"]["tensorqtl_perm"]["output"].format(
                            cell_type   = wildcards.cell_type,
                            norm_method = config["tensorQTL"]["norm_methods"][0],
                            geno_pc     = config["tensorQTL"]["geno_pcs"],
                            exp_pc      = config["exp_pc_map"][wildcards.cell_type]
                        )
    output: enrich = config["dev_specificity"]["pi1_enrich_jang"]["enrich"],
            pi1    = config["dev_specificity"]["pi1_enrich_jang"]["pi1"]
    resources: threads = 4, mem_mb = 60000, time="1:00:00"
    message: "Calc pi1 enrichments between fetal L1 sn-eQTL and Jang et al. 2026 L1 eQTL"
    singularity: config["containers"]["r_eqtl"]
    benchmark: "reports/benchmarks/dev_specificity.pi1_enrich_jang_{cell_type}_vs_{jang_cell_type}.txt"
    log:    config["dev_specificity"]["pi1_enrich_jang"]["log"]
    script: "../scripts/repl_sc_pi1_enrichment.R"

rule beta_correlation:
    input:  jang_l1_files     = get_jang_beta_corr_full_inputs(),
            jang_l1_top_files = get_jang_beta_corr_top_inputs(),
            pvar_file         = config["dev_specificity"]["beta_correlation"]["pvar_file"]
    output: config["dev_specificity"]["beta_correlation"]["out_file"]
    params: nom_dir = config["dev_specificity"]["beta_correlation"]["nom_dir"]
    singularity: config["containers"]["r_eqtl"]
    resources: threads = 4, mem_mb = 64000, time = "04:00:00"
    threads: 4
    message: "Compute fetal-adult eQTL beta correlations and flag opposite-effect pairs with allele-flip annotations"
    benchmark: "reports/benchmarks/dev_specificity.beta_correlation.txt"
    log:    config["dev_specificity"]["beta_correlation"]["log"]
    script: "../scripts/repl_sc_beta_correlation.R"

rule pi1_enrichment_report:
    # Note diff paths for output and out_file; Rmarkdown needs outfile to be relative to Rmd file
    input:  ziffra = ziffra_inputs,
            internal = internal_inputs,
            jang_pi1 = jang_pi1_inputs,
            beta_correlation = rules.beta_correlation.output[0],
            rmd_script = "scripts/repl_sc_report.Rmd"
    output: "reports/06QTL-REPLICATION/06replication_pi1_enrichment_report.html"
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
            output_file = "../reports/06QTL-REPLICATION/06replication_pi1_enrichment_report.html"
    singularity: config["containers"]["r_eqtl"]
    message: "Generate pi1 enrichments report report"
    resources: time="1:00:00"
    benchmark: "reports/benchmarks/06qtl_replication.pi1_enrichment_report.benchmark.txt"
    log:     "../results/00LOG/06QTL-REPLICATION/pi1_enrichment_report.log"
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

