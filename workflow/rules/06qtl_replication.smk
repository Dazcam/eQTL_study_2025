localrules: dwnld_obrien, dwnld_bryois, dwnld_ziffra, dwnld_wen, snp_lookup

rule dwnld_obrien:
    output: all_qtl = config["qtl_rep"]["dwnld_obrien"]["all_qtl"],
            top_qtl = config["qtl_rep"]["dwnld_obrien"]["top_qtl"]
    params: out_dir = config["qtl_rep"]["dwnld_obrien"]["out_dir"],
            web_link = config["qtl_rep"]["dwnld_obrien"]["web_link"],
            zip_out = config["qtl_rep"]["dwnld_obrien"]["zip_out"]
    message: "Download bulk brain gene eQTL file from O'Brien 2018, PMID:30419947"
    benchmark: "reports/benchmarks/qtl_replication.dwnld_obrien.txt"
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
            out_dir = config["qtl_rep"]["dwnld_bryois"]["out_dir"]
    envmodules: "compiler/gnu/7/3.0", "jq"
    message: "Download and cat single cell adult brain eQTL file from Bryois 2022, PMID:35915177"
    benchmark: "reports/benchmarks/qtl_replication.dwnld_bryois.txt"    
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
            touch {output}
            """        

rule dwnld_ziffra:
    output: config["qtl_rep"]["dwnld_ziffra"]["output"]
    params: web_link = config["qtl_rep"]["dwnld_ziffra"]["web_link"],
    message: "Download snATAC-seq peaks file from Ziffra 2021, PMID:34616060"
    benchmark: "reports/benchmarks/qtl_replication.dwnld_ziffra.txt"
    log:    config["qtl_rep"]["dwnld_ziffra"]["log"]
    shell:  """
            wget -O {output} {params.web_link} &>> {log}
            """

rule dwnld_wen:
    output: touch(config["qtl_rep"]["dwnld_wen"]["output"])
    params: json = config["qtl_rep"]["dwnld_wen"]["json"],
            out_dir = config["qtl_rep"]["dwnld_wen"]["out_dir"],
    envmodules: "compiler/gnu/7/3.0", "jq"
    message: "Download bulk brain eQTL file from Wen 2018, PMID:38781368"
    benchmark: "reports/benchmarks/qtl_replication.dwnld_wen.txt"
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

rule snp_lookup:
    input:  qtl_perm = config["qtl_rep"]["snp_lookup"]["qtl_perm"]
    output: config["qtl_rep"]["snp_lookup"]["out_file"]
    message: "Generate snp lookup for Ziffra ATAC-seq peaks; Net access req. - run local"
    singularity: config["containers"]["r_eqtl"]
    benchmark: "reports/benchmarks/qtl_replication.snp_lookup_{cell_type}_{norm_method}_genPC_{geno_pc}_expPC_{exp_pc}.txt"
    log:    config["qtl_rep"]["snp_lookup"]["log"]
    script: "../scripts/create_snp_lookup.R"     

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
    script: "../scripts/atac_enrichments_ziffra.R"
  
rule cat_nom_qtl:
    input:   parquets = expand(config["qtl_rep"]["cat_nom_qtl"]["parquet"], chr=range(1,23),allow_missing=True)
    output:  config["qtl_rep"]["cat_nom_qtl"]["output"]
    resources: threads = 1, mem_mb = 10000, time="1:00:00"
    message: "Combine per-chr parquet files into gzip file for pi1 enrichment"
    benchmark: "reports/benchmarks/qtl_replication.cat_nom_qtl_{cell_type}_{norm_method}_genPC_{geno_pc}_expPC_{exp_pc}.txt"
    log:     config["qtl_rep"]["cat_nom_qtl"]["log"]
    shell:   "scripts/cat_nom_qtl.py {output} {input.parquets} > {log} 2>&1"

rule pi1_enrich_obrien:
    # Currenly running all cell types together at a specific exp_pc, gen_pc, norm_method (norm method stil to be added)
    input:  public_all = rules.dwnld_obrien.output.all_qtl,
            public_top = rules.dwnld_obrien.output.top_qtl,
            qtl_all = rules.cat_nom_qtl.output,
            qtl_top = config["tensorQTL"]["tensorqtl_perm"]["output"]
    output: enrich = config["qtl_rep"]["pi1_enrich_obrien"]["enrich"],
            pi1 = config["qtl_rep"]["pi1_enrich_obrien"]["pi1"] 
    resources: threads = 4, mem_mb = 20000, time="1:00:00"
    message: "Calc pi1 enrichemnts between sn-eQTL and O'Brien 2018 bulk brain eQTL"
    singularity: config["containers"]["r_eqtl"]
    benchmark: "reports/benchmarks/qtl_replication.pi1_enrich_obrien_{cell_type}_{norm_method}_genPC_{geno_pc}_expPC_{exp_pc}.txt"
    log:    config["qtl_rep"]["pi1_enrich_obrien"]["log"]    
    script: "../scripts/pi1_enrichments_bulk.R"    

rule pi1_enrich_wen:
    # Currenly running all cell types together at a specific exp_pc, gen_pc, norm_method (norm method stil to be added)
    input:  public_all = config["qtl_rep"]["pi1_enrich_wen"]["public_all"],             
            public_top = config["qtl_rep"]["pi1_enrich_wen"]["public_top"],
            qtl_all = rules.cat_nom_qtl.output,
            qtl_top = config["tensorQTL"]["tensorqtl_perm"]["output"]
    output: enrich = config["qtl_rep"]["pi1_enrich_wen"]["enrich"],
            pi1 = config["qtl_rep"]["pi1_enrich_wen"]["pi1"]
    resources: threads = 4, mem_mb = 20000, time="1:00:00"
    message: "Calc pi1 enrichemnts between sn-eQTL and Wen 2024 bulk brain eQTL"
    singularity: config["containers"]["r_eqtl"]
    benchmark: "reports/benchmarks/qtl_replication.pi1_enrich_wen_{cell_type}_{norm_method}_genPC_{geno_pc}_expPC_{exp_pc}.txt"
    log:    config["qtl_rep"]["pi1_enrich_wen"]["log"]
    script: "../scripts/pi1_enrichments_bulk.R"
