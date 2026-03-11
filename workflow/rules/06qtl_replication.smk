configfile: "../config/config.yaml"

localrules: dwnld_obrien, dwnld_bryois, dwnld_ziffra, dwnld_wen, snp_lookup

CELL_EXP_PC = {ct: config["exp_pc_map"][ct] for ct in config["cell_types"]}
CELL_EXP_PC_TUPLES = [(ct, pc) for ct, pc in CELL_EXP_PC.items()]

ziffra_inputs = []
ziffra_indep_inputs = []
obrien_inputs = []
wen_inputs = []
bryois_inputs = []
fugita_inputs = []
internal_inputs = []
cat_nom_inputs = []

for ct, pc in CELL_EXP_PC_TUPLES:
    ziffra_inputs.extend(
        expand(
            config["qtl_rep"]["atac_enrich"]["out_file"],
            cell_type=[ct],
            exp_pc=[pc],
            geno_pc=config["tensorQTL"]["geno_pcs"],
            norm_method=config["tensorQTL"]["norm_methods"]
        )
    )

for ct, pc in CELL_EXP_PC_TUPLES:
    ziffra_indep_inputs.extend(
        expand(
            config["qtl_rep"]["atac_enrich_indep"]["out_file"],
            cell_type=[ct],
            exp_pc=[pc],
            geno_pc=config["tensorQTL"]["geno_pcs"],
            norm_method=config["tensorQTL"]["norm_methods"]
        )
    )

for ct, pc in CELL_EXP_PC_TUPLES:
    obrien_inputs.extend(
        expand(
            config["qtl_rep"]["pi1_enrich_obrien"]["pi1"],
            cell_type=[ct],
            exp_pc=[pc],
            geno_pc=config["tensorQTL"]["geno_pcs"],
            norm_method=config["tensorQTL"]["norm_methods"]
        )
    )

for ct, pc in CELL_EXP_PC_TUPLES:
    wen_inputs.extend(
        expand(
            config["qtl_rep"]["pi1_enrich_wen"]["pi1"],
            cell_type=[ct],
            exp_pc=[pc],
            geno_pc=config["tensorQTL"]["geno_pcs"],
            norm_method=config["tensorQTL"]["norm_methods"]
        )
    )

for ct, pc in CELL_EXP_PC_TUPLES:
    bryois_inputs.extend(
        expand(
            config["qtl_rep"]["pi1_enrich_bryois"]["pi1"],
            cell_type=[ct],
            exp_pc=[pc],
            geno_pc=config["tensorQTL"]["geno_pcs"],
            norm_method=config["tensorQTL"]["norm_methods"],
            ref_cell_type=config["cell_types_bryois"]
        )
    )

for ct, pc in CELL_EXP_PC_TUPLES:
    fugita_inputs.extend(
        expand(
            config["qtl_rep"]["pi1_enrich_fugita"]["pi1"],
            cell_type=[ct],
            exp_pc=[pc],
            geno_pc=config["tensorQTL"]["geno_pcs"],
            norm_method=config["tensorQTL"]["norm_methods"],
            fugita_cell_type=config["fugita_cell_types"]
        )
    )

internal_inputs = []

for ct, ct_pc in CELL_EXP_PC_TUPLES:
    for ref_ct, ref_pc in CELL_EXP_PC_TUPLES:
        internal_inputs.extend(
            expand(
                config["qtl_rep"]["pi1_enrich_internal"]["pi1"],
                cell_type=[ct],
                exp_pc=[ct_pc],
                ref_cell_type=[ref_ct],
                ref_exp_pc=[ref_pc],
                geno_pc=config["tensorQTL"]["geno_pcs"],
                norm_method=config["tensorQTL"]["norm_methods"],
            )
        )

cat_nom_inputs = []

for ct, pc in CELL_EXP_PC_TUPLES:
    cat_nom_inputs.extend(
        expand(
            config["qtl_rep"]["cat_nom_qtl"]["output"],
            cell_type=[ct],
            exp_pc=[pc],
            geno_pc=config["tensorQTL"]["geno_pcs"],
            norm_method=config["tensorQTL"]["norm_methods"],
        )
    )

rule all:
    input:
#        "reports/06QTL-REPLICATION/06replication_pi1_enrichment_report.html"
#        expand("../results/06QTL-REPLICATION/atac_enrich/{cell_type}_{norm_method}_genPC_{geno_pc}_expPC_{exp_pc}.ziffra_peak_enrich.rds", cell_type = config["cell_types"], norm_method = config["tensorQTL"]["norm_methods"][0], geno_pc = config["tensorQTL"]["geno_pcs"], exp_pc = config["tensorQTL"]["exp_pcs"][0])
#        expand(config["qtl_rep"]["atac_enrich"]["out_file"], zip,  cell_type=[ct for ct, pc in CELL_EXP_PC_TUPLES], exp_pc=[pc for ct, pc in CELL_EXP_PC_TUPLES], geno_pc=config["tensorQTL"]["geno_pcs"],norm_method=config["tensorQTL"]["norm_methods"])
        '../results/06QTL-REPLICATION/beta_cor/replication_beta_correlation_single.done',
#        ziffra_indep_inputs,
#        fugita_inputs,
#        internal_inputs 

rule dwnld_obrien:
    output: all_qtl = config["qtl_rep"]["dwnld_obrien"]["all_qtl"],
            top_qtl = config["qtl_rep"]["dwnld_obrien"]["top_qtl"]
    params: out_dir = config["qtl_rep"]["dwnld_obrien"]["out_dir"],
            web_link = config["qtl_rep"]["dwnld_obrien"]["web_link"],
            zip_out = config["qtl_rep"]["dwnld_obrien"]["zip_out"]
    message: "Download bulk brain gene eQTL file from O'Brien 2018, PMID:30419947"
    benchmark: "reports/benchmarks/06qtl_replication.dwnld_obrien.txt"
    log:    config["qtl_rep"]["dwnld_obrien"]["log"]
    shell:  """
            curl -L -b "cookies.txt" -A "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/91.0.4472.124 Safari/537.36" -o {params.zip_out} {params.web_link} &&
            unzip -j -o {params.zip_out} all_eqtls_gene.txt.gz -d {params.out_dir} &&
            unzip -j -o {params.zip_out} top_eqtls_gene.txt.gz -d {params.out_dir} &&
            rm -f {params.zip_out} 2>> {log}
            """
        
rule dwnld_bryois:
    output: touch(config["qtl_rep"]["dwnld_bryois"]["output"])
    params: json  = config["qtl_rep"]["dwnld_bryois"]["json"],
            out_dir = config["qtl_rep"]["dwnld_bryois"]["out_dir"],
            perm_file = config["qtl_rep"]["dwnld_bryois"]["perm_file"],
            web_link = config["qtl_rep"]["dwnld_bryois"]["web_link"]
    envmodules: "compiler/gnu/7/3.0", "jq"
    message: "Download and cat single cell adult brain eQTL file from Bryois 2022, PMID:35915177"
    benchmark: "reports/benchmarks/06qtl_replication.dwnld_bryois.txt"    
    log:    config["qtl_rep"]["dwnld_bryois"]["log"]
    shell:  """
            mkdir -p {params.out_dir}
            jq -r '.files.entries | to_entries[] | [.key, .value.links.content] | join(" ")' {params.json} > {out_dir}temp_download_list.txt
            while read -r file url; do
              [ -f {params.out_dir}/$file ] || curl -L -o {params.out_dir}/$file $url 2>> {log}
              echo "Downloaded: $file" >> {log}
            done < {out_dir}temp_download_list.txt
            rm -f {out_dir}temp_download_list.txt
        
            python scripts/concat_bryois.py {params.out_dir} 2>> {log}
            wget -O {params.perm_file} {params.web_link} &>> {log} 
            touch {output}
            """        

rule dwnld_ziffra:
    output: config["qtl_rep"]["dwnld_ziffra"]["output"]
    params: web_link = config["qtl_rep"]["dwnld_ziffra"]["web_link"],
    message: "Download snATAC-seq peaks file from Ziffra 2021, PMID:34616060"
    benchmark: "reports/benchmarks/06qtl_replication.dwnld_ziffra.txt"
    log:    config["qtl_rep"]["dwnld_ziffra"]["log"]
    shell:  """
            wget -O {output} {params.web_link} &>> {log}
            """

rule dwnld_wen:
    output: touch(config["qtl_rep"]["dwnld_wen"]["output"])
    params: json = config["qtl_rep"]["dwnld_wen"]["json"],
            out_dir = config["qtl_rep"]["dwnld_wen"]["out_dir"]
    envmodules: "compiler/gnu/7/3.0", "jq"
    message: "Download bulk brain eQTL file from Wen 2018, PMID:38781368"
    benchmark: "reports/benchmarks/06qtl_replication.dwnld_wen.txt"
    log:    config["qtl_rep"]["dwnld_wen"]["log"]
    shell:  """
	    mkdir -p {params.out_dir}
            jq -r '.files.entries | to_entries[] | [.key, .value.links.content] | join(" ")' {params.json} > {out_dir}temp_download_list.txt
            while read -r file url; do
              [ -f {params.out_dir}/$file ] || curl -L -o {params.out_dir}/$file $url 2>> {log}
              echo "Downloaded: $file" >> {log}
            done < {out_dir}temp_download_list.txt
            rm -f {out_dir}temp_download_list.txt
	
            touch {output}
            """


## These are optional rules to run the ziffra enrichemnts with proxies -----
## Need to add distinct labels in config if using
#rule get_perm_ld_proxies:
#    input:
#            perm = lambda wc: config["tensorQTL"]["tensorqtl_perm"]["output"].format(
#                cell_type = wc.cell_type,
#                norm_method = config["tensorQTL"]["norm_methods"][0],
#                geno_pc     = config["tensorQTL"]["geno_pcs"],
#                exp_pc      = config["tensorQTL"]["exp_pcs"][0]
#            )
#    output: config["qtl_rep"]["get_perm_ld_proxies"]["output"]
#    params: proxy_dir = config["qtl_rep"]["get_perm_ld_proxies"]["proxy_dir"],
#            ldlink_token = config["qtl_rep"]["get_perm_ld_proxies"]["ldlink_token"]
#    singularity: config["containers"]["r_eqtl"]
#    message: "Get LD proxies for tensorQTL perm eQTL for atac enrich rule"
#    benchmark: "reports/benchmarks/06qtl_replication.get_perm_ld_proxies_{cell_type}.txt"
#    log:    config["qtl_rep"]["get_perm_ld_proxies"]["log"]
#    script:  "../scripts/replication_get_perm_ld_proxies.R"

#rule snp_lookup_proxies:
#    input:  qtl_perm = config["qtl_rep"]["snp_lookup"]["qtl_perm"],
#            proxies = rules.get_perm_ld_proxies.output
#    output: config["qtl_rep"]["snp_lookup"]["out_file"]
#    params: config["qtl_rep"]["get_perm_ld_proxies"]["proxy_dir"]
#    message: "Generate snp lookup for cell-specific sig. eQTL and LD proxies; Net access req."
#    singularity: config["containers"]["r_eqtl"]
#    benchmark: "reports/benchmarks/06qtl_replication.snp_lookup_{cell_type}_{norm_method}_genPC_{geno_pc}_expPC_{exp_pc}.txt"
#    log:    config["qtl_rep"]["snp_lookup"]["log"]
#    script: "../scripts/replication_create_snp_lookup_proxies.R"     

#rule atac_enrich_proxies:
#    input:  peaks = rules.dwnld_ziffra.output,
#            snp_file = rules.snp_lookup.output 
#    output: config["qtl_rep"]["atac_enrich"]["out_file"]
#    params: peak_dir = config["qtl_rep"]["atac_enrich"]["peak_dir"]
#    singularity: config["containers"]["r_eqtl"]    
#    message: "Test for sig. eQTL enrichments in Ziffra ATAC-seq peaks"
#    benchmark: "reports/benchmarks/06qtl_replication.atac_enrich_{cell_type}_{norm_method}_genPC_{geno_pc}_expPC_{exp_pc}.txt"
#    log:    config["qtl_rep"]["atac_enrich"]["log"]
#    script: "../scripts/replication_atac_enrichments_proxies.R"
  
##  --------------------------------------------------------------------------

rule snp_lookup:
    input:  qtl_perm = config["qtl_rep"]["snp_lookup"]["qtl_perm"]
    output: config["qtl_rep"]["snp_lookup"]["out_file"]
    message: "Generate snp lookup for Ziffra ATAC-seq peaks; Net access req. - run local"
    singularity: config["containers"]["r_eqtl"]
    benchmark: "reports/benchmarks/qtl_replication.snp_lookup_{cell_type}_{norm_method}_genPC_{geno_pc}_expPC_{exp_pc}.txt"
    log:    config["qtl_rep"]["snp_lookup"]["log"]
    script: "../scripts/replication_create_snp_lookup.R"     

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
    script: "../scripts/replication_atac_enrichments.R"

rule snp_lookup_indep:
    input:  qtl_indep = config["qtl_rep"]["snp_lookup_indep"]["qtl_indep"]
    output: config["qtl_rep"]["snp_lookup_indep"]["out_file"]
    message: "Generate snp lookup (indep eQTL) for Ziffra ATAC-seq peaks; Net access req. - run local"
    singularity: config["containers"]["r_eqtl"]
    benchmark: "reports/benchmarks/qtl_replication.snp_lookup_indep_{cell_type}_{norm_method}_genPC_{geno_pc}_expPC_{exp_pc}.txt"
    log:    config["qtl_rep"]["snp_lookup_indep"]["log"]
    script: "../scripts/replication_create_snp_lookup_indep.R"     

rule atac_enrich_indep:
    input:  qtl_indep = config["qtl_rep"]["snp_lookup_indep"]["qtl_indep"], 
            peaks = config["qtl_rep"]["dwnld_ziffra"]["output"],
            snp_file = rules.snp_lookup_indep.output 
    output: config["qtl_rep"]["atac_enrich_indep"]["out_file"]
    params: peak_dir = config["qtl_rep"]["atac_enrich"]["peak_dir"]
    singularity: config["containers"]["r_eqtl"]
    resources: time="1:00:00"

    message: "Test for sig. eQTL indep enrichments in Ziffra ATAC-seq peaks"
    benchmark: "reports/benchmarks/qtl_replication.atac_enrich_indep_{cell_type}_{norm_method}_genPC_{geno_pc}_expPC_{exp_pc}.txt"
    log:    config["qtl_rep"]["atac_enrich_indep"]["log"]
    script: "../scripts/replication_atac_enrichments_indep.R"


rule cat_nom_qtl:
    input:  parquets = lambda wc: expand(config["qtl_rep"]["cat_nom_qtl"]["parquet"], cell_type=wc.cell_type, norm_method=wc.norm_method, geno_pc=wc.geno_pc, exp_pc=wc.exp_pc, chr=range(1, 23), allow_missing=True)
    output:  config["qtl_rep"]["cat_nom_qtl"]["output"]
    resources: threads = 1, mem_mb = 10000, time="1:00:00"
    message: "Combine per-chr parquet files into gzip file for pi1 enrichment"
    benchmark: "reports/benchmarks/06qtl_replication.cat_nom_qtl_{cell_type}_{norm_method}_genPC_{geno_pc}_expPC_{exp_pc}.txt"
    log:     config["qtl_rep"]["cat_nom_qtl"]["log"]
    shell:   "scripts/replication_cat_nom_qtl.py {output} {input.parquets} > {log} 2>&1"

rule pi1_enrich_obrien:
    # Currenly running all cell types together at a specific exp_pc, gen_pc, norm_method (norm method stil to be added)
    input:  public_all = rules.dwnld_obrien.output.all_qtl,
            public_top = rules.dwnld_obrien.output.top_qtl,
            qtl_all = rules.cat_nom_qtl.output,
            qtl_top = config["tensorQTL"]["tensorqtl_perm"]["output"]
    output: enrich = config["qtl_rep"]["pi1_enrich_obrien"]["enrich"],
            pi1 = config["qtl_rep"]["pi1_enrich_obrien"]["pi1"] 
    resources: threads = 4, mem_mb = 20000, time="1:00:00"
    message: "Calc pi1 enrichments between sn-eQTL and O'Brien 2018 bulk brain eQTL"
    singularity: config["containers"]["r_eqtl"]
    benchmark: "reports/benchmarks/06qtl_replication.pi1_enrich_obrien_{cell_type}_{norm_method}_genPC_{geno_pc}_expPC_{exp_pc}.txt"
    log:    config["qtl_rep"]["pi1_enrich_obrien"]["log"]    
    script: "../scripts/replication_pi1_enrichment.R"    

rule pi1_enrich_wen:
    # Currenly running all cell types together at a specific exp_pc, gen_pc, norm_method (norm method stil to be added)
    input:  public_all = config["qtl_rep"]["pi1_enrich_wen"]["public_all"],             
            public_top = config["qtl_rep"]["pi1_enrich_wen"]["public_top"],
            qtl_all = rules.cat_nom_qtl.output,
            qtl_top = config["tensorQTL"]["tensorqtl_perm"]["output"]
    output: enrich = config["qtl_rep"]["pi1_enrich_wen"]["enrich"],
            pi1 = config["qtl_rep"]["pi1_enrich_wen"]["pi1"]
    resources: threads = 4, mem_mb = 20000, time="1:00:00"
    message: "Calc pi1 enrichments between sn-eQTL and Wen 2024 bulk brain eQTL"
    singularity: config["containers"]["r_eqtl"]
    benchmark: "reports/benchmarks/06qtl_replication.pi1_enrich_wen_{cell_type}_{norm_method}_genPC_{geno_pc}_expPC_{exp_pc}.txt"
    log:    config["qtl_rep"]["pi1_enrich_wen"]["log"]
    script: "../scripts/replication_pi1_enrichment.R"

rule pi1_enrich_bryois:
    input:  public_all = config["qtl_rep"]["pi1_enrich_bryois"]["public_all"],
            public_top = config["qtl_rep"]["pi1_enrich_bryois"]["public_top"],
            qtl_all = rules.cat_nom_qtl.output,
            qtl_top = config["tensorQTL"]["tensorqtl_perm"]["output"]
    output: enrich = config["qtl_rep"]["pi1_enrich_bryois"]["enrich"],
            pi1    = config["qtl_rep"]["pi1_enrich_bryois"]["pi1"]
    resources: threads = 4, mem_mb = 20000, time="1:00:00"
    message: "Calc pi1 enrichments between sn-eQTL and Bryois 2022 single-cell eQTL"
    singularity: config["containers"]["r_eqtl"]
    benchmark: "reports/benchmarks/06qtl_replication.pi1_enrich_bryois_{cell_type}_vs_{ref_cell_type}_{norm_method}_genPC_{geno_pc}_expPC_{exp_pc}.txt"
    log:    config["qtl_rep"]["pi1_enrich_bryois"]["log"]
    script: "../scripts/replication_pi1_enrichment.R"    


rule pi1_enrich_fugita:
    input: public_all = config["qtl_rep"]["pi1_enrich_fugita"]["public_all"],
           public_top = config["qtl_rep"]["pi1_enrich_fugita"]["public_top"],
           qtl_all = rules.cat_nom_qtl.output,
           qtl_top = config["tensorQTL"]["tensorqtl_perm"]["output"]
    output: enrich = config["qtl_rep"]["pi1_enrich_fugita"]["enrich"],
            pi1 = config["qtl_rep"]["pi1_enrich_fugita"]["pi1"]
    resources: threads = 4, mem_mb = 20000, time="1:00:00"
    message: "Calc pi1 enrichments between sn-eQTL and Fugita 2024 single-cell eQTL"
    singularity: config["containers"]["r_eqtl"]
    benchmark: "reports/benchmarks/06qtl_replication.pi1_enrich_fugita_{cell_type}_vs_{fugita_cell_type}_{norm_method}_genPC_{geno_pc}_expPC_{exp_pc}.txt"
    log:    config["qtl_rep"]["pi1_enrich_fugita"]["log"]
    script: "../scripts/replication_pi1_enrichment.R"
    
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
    script: "../scripts/replication_pi1_enrichment.R"    

rule beta_correlation_all:
    output: '../results/06QTL-REPLICATION/beta_cor/replication_beta_correlation.tsv'
    params: in_dir = '../results/05TENSORQTL/tensorqtl_perm/',
            fugita_dir = '../resources/public_datasets/fugita_2024/'
    singularity: config["containers"]["r_eqtl"]
    resources: threads = 6, mem_mb = 96000
    message: "Generating data for eQTL beta correlation analysis for replication stage"
    benchmark: "reports/benchmarks/06qtl_replication.pi1_beta_correlation_all.txt"
    log:    "../results/00LOG/06QTL-REPLICATION/beta_correlation_all.log"
    script: "../scripts/replication_beta_correlation_all.R"

rule beta_correlation_single:
    output: '../results/06QTL-REPLICATION/beta_cor/replication_beta_correlation_single.done'
    params: in_dir = '../results/05TENSORQTL/tensorqtl_perm/',
            fugita_dir = '../resources/public_datasets/fugita_2024/'
    singularity: config["containers"]["r_eqtl"]
    resources: threads = 6, mem_mb = 96000, time="1-0:00:00"
    message: "Generating data for eQTL beta correlation analysis for replication (single N pairs)"
    benchmark: "reports/benchmarks/06qtl_replication.pi1_beta_correlation_single.txt"
    log:    "../results/00LOG/06QTL-REPLICATION/beta_correlation_single.log"
    script: "../scripts/replication_beta_correlation_single.R"

rule pi1_enrichment_report:
    # Note diff paths for output and out_file; Rmarkdown needs outfile to be relative to Rmd file
    input:  ziffra = ziffra_inputs,
            ziffra_indep = ziffra_indep_inputs,
#            obrien = obrien_inputs,
#            wen = wen_inputs,
#            bryois = bryois_inputs,
            fugita = fugita_inputs,
            internal = internal_inputs,
            beta_cor = '../results/06QTL-REPLICATION/beta_cor/replication_beta_correlation.tsv',
            beta_cor_all = '../results/06QTL-REPLICATION/beta_cor/replication_beta_correlation_single.done',
            rmd_script = "scripts/replication_pi1_enrichment_report.Rmd"
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
