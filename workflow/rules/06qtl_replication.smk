localrules: dwnld_obrien, dwnld_bryois

rule dwnld_obrien:
    output: config["qtl_rep"]["dwnld_obrien"]["output"]
    params: out_dir = config["qtl_rep"]["dwnld_obrien"]["out_dir"],
            web_link = config["qtl_rep"]["dwnld_obrien"]["web_link"],
            zip_out = config["qtl_rep"]["dwnld_obrien"]["zip_out"]
    message: "Download bulk brain gene eQTL file from O'Brien 2018, PMID: 30419947"
    benchmark: "reports/benchmarks/qtl_replication.dwnld_obrien.txt"
    log:    config["qtl_rep"]["dwnld_obrien"]["log"]
    shell:  """
            curl -L -b "cookies.txt" -A "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/91.0.4472.124 Safari/537.36" -o {params.zip_out} {params.web_link} &&
            unzip -j -o {params.zip_out} all_eqtls_gene.txt.gz -d {params.out_dir} &&
            rm -f {params.zip_out} 2>> {log}
            """
        
rule dwnld_bryois:
    output: touch(config["qtl_rep"]["dwnld_bryois"]["output"])
    params: json  = config["qtl_rep"]["dwnld_bryois"]["json"],
            out_dir = config["qtl_rep"]["dwnld_bryois"]["out_dir"]
    envmodules: "compiler/gnu/7/3.0", "jq"
    message: "Download single nuclei gene eQTL file from Bryois 2022, PMID: 35915177"
    benchmark: "reports/benchmarks/qtl_replication.dwnld_bryois.txt"    
    log:    config["qtl_rep"]["dwnld_bryois"]["log"]
    shell:
        """
        mkdir -p {params.out_dir}
        jq -r '.files.entries | to_entries[] | [.key, .value.links.content] | join(" ")' {params.json} > temp_download_list.txt
        while read -r file url; do
            [ -f {params.out_dir}/$file ] || curl -L -o {params.out_dir}/$file $url 2>> {log}
            echo "Downloaded: $file" >> {log}
        done < temp_download_list.txt
        rm -f temp_download_list.txt
        touch {output}
        """

rule dwnld_ziffra:
    output: config["qtl_rep"]["dwnld_ziffra"]["output"]
    params: out_dir = config["qtl_rep"]["dwnld_obrien"]["out_dir"],
            web_link = config["qtl_rep"]["dwnld_obrien"]["web_link"],
            zip_out = config["qtl_rep"]["dwnld_obrien"]["zip_out"]
    message: "Download bulk brain gene eQTL file from O'Brien 2018, PMID: 30419947"
    benchmark: "reports/benchmarks/qtl_replication.dwnld_obrien.txt"
    log:    config["qtl_rep"]["dwnld_obrien"]["log"]
    shell:  """
            curl -L -b "cookies.txt" -A "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/91.0.4472.124 Safari/537.36" -o {params.zip_out} {params.web_link} &&
            unzip -j -o {params.zip_out} all_eqtls_gene.txt.gz -d {params.out_dir} &&
            rm -f {params.zip_out} 2>> {log}
            """
     
