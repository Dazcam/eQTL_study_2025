configfile: '../config/config.yaml'

configfile: "../config/config.yaml"

localrules: dwnld_obrien, dwnld_bryois, dwnld_ziffra, dwnld_wen, snp_lookup

CELL_EXP_PC = {ct: config["exp_pc_map"][ct] for ct in config["cell_types"]}
CELL_EXP_PC_TUPLES = [(ct, pc) for ct, pc in CELL_EXP_PC.items()]

ziffra_all_inputs = []

for ct, pc in CELL_EXP_PC_TUPLES:
    ziffra_all_inputs.extend(
        expand(
            config["qtl_rep"]["atac_enrich_all"]["out_file"],
            cell_type=[ct],
            exp_pc=[pc],
            geno_pc=config["tensorQTL"]["geno_pcs"],
            norm_method=config["tensorQTL"]["norm_methods"]
        )
    )

rule all:
    input:
        ziffra_all_inputs
#        "../results/20MVARS/mVars_enrichment_bars.pdf"
#        expand("../results/06QTL-REPLICATION/atac_enrich_indep/{cell_type}_{norm_method}_genPC_{geno_pc}_expPC_{exp_pc}.ziffra_peak_enrich_indep.rds", cell_type = config["cell_types"], norm_method = config["tensorQTL"]["norm_methods"][0], geno_pc = config["tensorQTL"]["geno_pcs"], exp_pc = config["tensorQTL"]["exp_pcs"][0])    
#        "../results/20BETA_COR/replication_beta_correlation.tsv"

#rule mVars_enrich:
#    input: mvars = "../resources/public_datasets/Won_2025/emVar_list.xlsx"
#    output: "../results/20MVARS/mVars_enrichment_bars.pdf"
#    log:    "../results/00LOG/20MVARS/mVars_enrichment.log"
#    resources: threads = 6, mem_mb = 96000, time="1:00:00"
#    singularity: config["containers"]["r_eqtl"]
#    message: "Testing for enrichment of mVars in FDR sig. eQTL"
#    script:  "../scripts/mVar_enrichment.R"    

#rule snp_lookup_indep:
#    input:  qtl_indep = config["qtl_rep"]["snp_lookup_indep"]["qtl_indep"]
#    output: config["qtl_rep"]["snp_lookup_indep"]["out_file"]
#    message: "Generate snp lookup (indep eQTL) for Ziffra ATAC-seq peaks; Net access req. - run local"
#    singularity: config["containers"]["r_eqtl"]
#    benchmark: "reports/benchmarks/qtl_replication.snp_lookup_indep_{cell_type}_{norm_method}_genPC_{geno_pc}_expPC_{exp_pc}.txt"
#    log:    config["qtl_rep"]["snp_lookup_indep"]["log"]
#    script: "../scripts/replication_create_snp_lookup_indep.R"     

#rule atac_enrich_indep:
#    input:  qtl_indep = config["qtl_rep"]["snp_lookup_indep"]["qtl_indep"], 
#            peaks = config["qtl_rep"]["dwnld_ziffra"]["output"],
#            snp_file = rules.snp_lookup_indep.output 
#    output: config["qtl_rep"]["atac_enrich_indep"]["out_file"]
#    params: peak_dir = config["qtl_rep"]["atac_enrich"]["peak_dir"]
#    singularity: config["containers"]["r_eqtl"]
#    message: "Test for sig. eQTL indep enrichments in Ziffra ATAC-seq peaks"
#    benchmark: "reports/benchmarks/qtl_replication.atac_enrich_indep_{cell_type}_{norm_method}_genPC_{geno_pc}_expPC_{exp_pc}.txt"
#    log:    config["qtl_rep"]["atac_enrich_indep"]["log"]
#    script: "../scripts/replication_atac_enrichments_indep.R"

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
    script: "../scripts/replication_create_snp_lookup.R"

rule atac_enrich_all:
    input:  qtl_perm = config["qtl_rep"]["snp_lookup"]["qtl_perm"],
            peaks = rules.dwnld_ziffra.output,
            snp_file = rules.snp_lookup.output
    output: config["qtl_rep"]["atac_enrich_all"]["out_file"]
    params: peak_dir = config["qtl_rep"]["atac_enrich_all"]["peak_dir"]
    singularity: config["containers"]["r_eqtl"]
    message: "Test for sig. eQTL enrichments in Ziffra ATAC-seq peaks"
    benchmark: "reports/benchmarks/qtl_replication.atac_enrich_all_{cell_type}_{norm_method}_genPC_{geno_pc}_expPC_{exp_pc}.txt"
    log:    config["qtl_rep"]["atac_enrich_all"]["log"]
    script: "../scripts/replication_atac_enrichments_all_L1.R"
