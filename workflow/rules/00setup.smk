import os

localrules: create_parse_json, setup_liftover, get_gwas_sumstats, get_ldsc


##  Settings  ------------------------------------------------------------------------------
SET    = config["setup"]
PLATES = list(config["fastq"].keys())


##  Functions  -----------------------------------------------------------------------------
def ref_dir(build):
    """Top-level LDSC reference directory for a build (parent of baseline_v1.2/)."""
    return os.path.dirname(os.path.dirname(SET[f"get_{build}_refs"]["baseline"]))


##  Rules  ---------------------------------------------------------------------------------

rule get_containers:
    output:  SET["get_containers"]["output"],
    params:  script = SET["get_containers"]["script"]
    message: "Downloading singulairty / apptainer containers"
    benchmark: SET["get_containers"]["benchmark"]
    log:     SET["get_containers"]["log"]
    shell:
             r"""
             {params.script} > {log} 2>&1
             """

rule create_parse_json:
    output:  SET["create_parse_json"]["output"],
    params:  dirs=lambda wildcards: config["fastq"][wildcards.plate],
             script = SET["create_parse_json"]["script"]
    message:  "Create json to map fastq files to sample IDs for Parse alignment"
    benchmark: SET["create_parse_json"]["benchmark"]
    log:     SET["create_parse_json"]["log"]
    shell:
             """
             python {params.script} \
             --plate {wildcards.plate} \
             --fastq_dirs {params.dirs} > {log} 2>&1
             """

rule setup_liftover:
    output:  SET["liftover"]["output"]
    params:  outdir = SET["liftover"]["outdir"]
    benchmark: SET["liftover"]["benchmark"]
    log:     SET["liftover"]["log"]
    shell:
             """
             exec > {log} 2>&1
             wget -nc -P {params.outdir} http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/liftOver
             wget -nc -P {params.outdir} http://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz
             wget -nc -P {params.outdir} http://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz
        
             chmod +x {params.outdir}/liftOver
             touch {output}
             """

# Requires scanpy env so need to set up envs first
rule get_fugita_data:
    output:  SET["get_fugita_data"]["output"]
    params:  outdir=SET["get_fugita_data"]["outdir"],
             token=SET["get_fugita_data"]["token"],
             script = SET["get_fugita_data"]["script"]
    conda:   config["scanpy"]["env"]
    message: "Downloading Fugita et al. cell cell type and subtype data from Synapse"
    benchmark: SET["get_fugita_data"]["benchmark"]
    log:     SET["get_fugita_data"]["log"]
    shell:
             """
             python {params.script} \
             --token {params.token} \
             --outdir {params.outdir} > {log} 2>&1
             touch {output}
             """

rule get_jang_singlebrain:
    output:  SET["get_jang_singlebrain"]["output"]
    params:  outdir=SET["get_jang_singlebrain"]["outdir"]
    conda:   config["scanpy"]["env"]
    message: "Downloading Jang et al. SingleBrain eQTL summary statistics from Zenodo"
    benchmark: SET["get_jang_singlebrain"]["benchmark"]
    log:     SET["get_jang_singlebrain"]["log"]
    shell:
             """
             mkdir -p {params.outdir}
             zenodo_get 10.5281/zenodo.14908182 -o {params.outdir} > {log} 2>&1
             touch {output}
             """

rule parse_nott_peaks:
    input:  xlsx = SET["parse_nott_peaks"]["xlsx_file"]
    output: SET["parse_nott_peaks"]["out_file"]
    singularity: config["containers"]["r_eqtl"]
    message: "Parse Nott et al. 2019 Table S5 (8 sheets) into a combined hg19 BED"
    benchmark: SET["parse_nott_peaks"]["benchmark"]
    log:    SET["parse_nott_peaks"]["log"]
    script: SET["parse_nott_peaks"]["script"]

rule liftover_nott_peaks:
    # Liftover Nott et al. 2019 peaks from hg19 to hg38
    input: sentinel = SET["liftover"]["output"],
           bed   = rules.parse_nott_peaks.output,
           chain = SET["liftover_nott_hg38"]["chain"]
    output:  SET["liftover_nott_hg38"]["output"]
    singularity: config["containers"]["ubuntu"]
    params: unlifted = SET["liftover_nott_hg38"]["unlifted"],
            liftover = config["tools"]["liftover"]
    benchmark: SET["liftover_nott_hg38"]["benchmark"]
    message: "Lifting over Nott et al. 2019 peaks from hg19 to hg38"
    log: SET["liftover_nott_hg38"]["log"]
    shell: "{params.liftover} {input.bed} {input.chain} {output} {params.unlifted} 2> {log}"


#rule get_obrien_data:
#    output: all_qtl = SET["get_obrien_data"]["all_qtl"],
#            top_qtl = SET["get_obrien_data"]["top_qtl"]
#    params: out_dir = SET["get_obrien_data"]["out_dir"],
#            web_link = SET["get_obrien_data"]["web_link"],
#            zip_out = SET["get_obrien_data"]["zip_out"]
#    message: "Download bulk brain gene eQTL files from O'Brien 2018, PMID:30419947"
#    log:    SET["get_obrien_data"]["log"]
#    shell:  """
#            curl -L -b "cookies.txt" -A "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/91.0.4472.124 Safari/537.36" -o {params.zip_out} {params.web_link} &&
#            unzip -j -o {params.zip_out} all_eqtls_gene.txt.gz -d {params.out_dir} &&
#            unzip -j -o {params.zip_out} top_eqtls_gene.txt.gz -d {params.out_dir} &&
#            rm -f {params.zip_out} 2>> {log}
#            """

rule get_obrien_supp:
    output: SET["get_obrien_supp"]["output"],
    params: web_link = SET["get_obrien_supp"]["web_link"],
    message: "Download bulk eQTL supp tables for from O'Brien 2018, PMID:30419947"
    benchmark: SET["get_obrien_supp"]["benchmark"]
    log:    SET["get_obrien_supp"]["log"]
    shell:   "wget -O {output} {params.web_link} > {log} 2>&1"

rule get_gwas_sumstats:
    output:  SET["get_gwas_sumstats"]["output"]
    wildcard_constraints: gwas = "|".join(config["gwas"])
    params:  url        = lambda w: config["gwas"][w.gwas]["url"],
             file       = lambda w: config["gwas"][w.gwas]["file"],
             md5        = lambda w: config["gwas"][w.gwas]["md5"],
             dir        = SET["get_gwas_sumstats"]["dir"],
             token_file = SET["get_gwas_sumstats"]["figshare_token"]
    message: "Download {wildcards.gwas} GWAS summary statistics from figshare"
    benchmark: SET["get_gwas_sumstats"]["benchmark"]
    log:     SET["get_gwas_sumstats"]["log"]
    shell:
             """
             exec > {log} 2>&1
             if [ ! -s {params.token_file} ]; then
                 echo "ERROR: figshare token file {params.token_file} missing or empty"
                 exit 1
             fi
             mkdir -p {params.dir}
             curl -fsSL -H "Authorization: token $(cat {params.token_file})" \
                  -o {params.dir}/{params.file} "{params.url}"
             ls -lh {params.dir}/{params.file}
             md5sum {params.dir}/{params.file}
             echo "{params.md5}  {params.dir}/{params.file}" | md5sum -c -
             touch {output}
             """

rule get_hg38_refs:
    output: baseline = expand(SET["get_hg38_refs"]["baseline"], chr=range(1, 23)),
            weights = expand(SET["get_hg38_refs"]["weights"], chr=range(1, 23)),
            frq = expand(SET["get_hg38_refs"]["frq"], chr=range(1, 23)),
            bim = expand(SET["get_hg38_refs"]["bim"], chr=range(1, 23)),
            bed = expand(SET["get_hg38_refs"]["bed"], chr=range(1, 23)),
            fam = expand(SET["get_hg38_refs"]["fam"], chr=range(1, 23)),
            snps = SET["get_hg38_refs"]["snps"]
    params: out_dir = ref_dir("hg38"),
            script  = SET["get_hg38_refs"]["script"]
    benchmark: SET["get_hg38_refs"]["benchmark"]
    log:    SET["get_hg38_refs"]["log"]
    shell: "{params.script} {params.out_dir} > {log} 2>&1"

rule get_hg19_refs:
    output: baseline = expand(SET["get_hg19_refs"]["baseline"], chr=range(1, 23)),
            weights = expand(SET["get_hg19_refs"]["weights"], chr=range(1, 23)),
            frq = expand(SET["get_hg19_refs"]["frq"], chr=range(1, 23)),
            bim = expand(SET["get_hg19_refs"]["bim"], chr=range(1, 23)),
            snps = SET["get_hg19_refs"]["snps"]
    params: out_dir = ref_dir("hg19"),
            script  = SET["get_hg19_refs"]["script"]
    benchmark: SET["get_hg19_refs"]["benchmark"]
    log:    SET["get_hg19_refs"]["log"]
    shell: "{params.script} {params.out_dir} > {log} 2>&1"

rule make_hm3_snplist:
    # HapMap3 no-MHC SNPs in the hg38 weights, with A1 and A2 from the hg38 bim files
    input:   bim = expand(SET["get_hg38_refs"]["bim"], chr=range(1, 23)),
             weights = expand(SET["get_hg38_refs"]["weights"], chr=range(1, 23))
    output:  SET["make_hm3_snplist"]["output"]
    benchmark: SET["make_hm3_snplist"]["benchmark"]
    log:     SET["make_hm3_snplist"]["log"]
    shell:
             """
             exec 2> {log}
             cat {input.bim} | awk 'BEGIN {{OFS="\t"; print "SNP", "A1", "A2"}} {{print $2, $5, $6}}' > {output}.bim_tmp
             first=1
             for f in {input.weights}; do
                 if [ $first = 1 ]; then zcat $f; first=0; else zcat $f | tail -n +2; fi
             done > {output}.weights_tmp
             awk 'NR==FNR {{snp[$2]; next}} $1 in snp' {output}.weights_tmp {output}.bim_tmp > {output}
             rm {output}.bim_tmp {output}.weights_tmp
             """

rule get_ldsc:
    # Clone ldsc at the pinned commit; an existing clone is checked out, not cloned again
    output:  SET["get_ldsc"]["output"]
    params:  repo = SET["get_ldsc"]["repo"],
             commit = SET["get_ldsc"]["commit"],
             dir = SET["get_ldsc"]["dir"]
    message: "Clone ldsc"
    benchmark: SET["get_ldsc"]["benchmark"]
    log:     SET["get_ldsc"]["log"]
    shell:
             """
             exec > {log} 2>&1
             if [ ! -d {params.dir}/.git ]; then
                 git clone {params.repo} {params.dir}
             fi
             git -C {params.dir} checkout {params.commit}
             git -C {params.dir} log -1 --format="%H %ad %s"
             touch {output}
             """


##  Final outputs, collected by rule all in the Snakefile  ---------------------------------
TARGETS = globals().setdefault("TARGETS", {})
TARGETS["00setup"] = [rules.get_containers.output,
                      expand(rules.create_parse_json.output, plate = PLATES),
                      rules.setup_liftover.output,
                      rules.get_fugita_data.output,
                      rules.get_jang_singlebrain.output,
                      rules.get_obrien_supp.output,
                      expand(rules.get_gwas_sumstats.output, gwas = config["gwas"]),
                      rules.liftover_nott_peaks.output,
                      rules.get_hg38_refs.output,
                      rules.get_hg19_refs.output,
                      rules.make_hm3_snplist.output,
                      rules.get_ldsc.output]
