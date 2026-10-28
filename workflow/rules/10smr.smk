configfile: "../config/config.yaml"


rule get_smr_binary:
    output: config["smr"]["get_smr_binary"]["smr"]
    message: "Download the smr binary file"
    benchmark: "reports/benchmarks/10smr.get_smr_binary.benchmark.txt"
    log: config["smr"]["get_smr_binary"]["log"]
    shell:
        """
        (
        wget https://yanglab.westlake.edu.cn/software/smr/download/smr-1.4.0-linux-x86_64.zip
        unzip smr-1.4.0-linux-x86_64.zip 'smr-1.4.0-linux-x86_64/smr'
        mv smr-1.4.0-linux-x86_64/smr {output}
        rm -rf smr-1.4.0-linux-x86_64.zip smr-1.4.0-linux-x86_64
        ) > {log} 2>&1
        """

rule cat_refs:
    input:  frq = expand(config["smr"]["cat_refs"]["frq_in"], chr = range(1, 23)),
            bed = expand(config["smr"]["cat_refs"]["bed_in"], chr = range(1, 23)),
            bim = expand(config["smr"]["cat_refs"]["bim_in"], chr = range(1, 23)),
            fam = config["smr"]["cat_refs"]["fam_in"]
    output: frq = config["smr"]["cat_refs"]["cat_frq"],
            bed = config["smr"]["cat_refs"]["cat_bed"],
            bim = config["smr"]["cat_refs"]["cat_bim"],
            fam = config["smr"]["cat_refs"]["cat_fam"]
    message: "Cat hg38 refs into a single file"
    benchmark: "reports/benchmarks/10smr.cat_refs.benchmark.txt"
    log: config["smr"]["cat_refs"]["log"]
    shell:
        """
        head -n 1 {input.frq[0]} > {output.frq}
        for file in {input.frq}; do
            tail -n +2 $file >> {output.frq}
        done

        cat {input.bed} > {output.bed}
        cat {input.bim} > {output.bim}
        cp {input.fam} {output.fam}

        echo "Concatenated frq, bed, bim, fam files" > {log}
        """

def tensorqtl_nom_dir(w):
    """tensorQTL nominal output directory for one cell type (first norm method)."""
    norm_method = config['tensorQTL']['norm_methods'][0]
    geno_pc = config['tensorQTL']['geno_pcs']
    exp_pc = config['exp_pc_map'][w.cell_type]
    return (f"../results/05TENSORQTL/tensorqtl_nom/"
            f"{w.cell_type}_{norm_method}_genPC_{geno_pc}_expPC_{exp_pc}/")

rule cat_tensorqtl_nom_snps:
    input:  lambda w: tensorqtl_nom_dir(w) +
                f"{w.cell_type}_{config['tensorQTL']['norm_methods'][0]}_nom.cis_qtl_pairs.4.parquet"
    output: cat_snps = config["smr"]["cat_tensorqtl_nom_snps"]["cat_snps"],
            summary = config["smr"]["cat_tensorqtl_nom_snps"]["summary"]
    params: dir = tensorqtl_nom_dir,
            norm_method = config['tensorQTL']['norm_methods'][0]
    message: "Cat TensorQTL parquet files into single file for SMR"
    benchmark: "reports/benchmarks/10smr.cat_tensorqtl_nom_snps_{cell_type}.benchmark.txt"
    log: config["smr"]["cat_tensorqtl_nom_snps"]["log"]
    shell:
        """
        python scripts/smr_cat_tensorqtl_nom_snps.py \
          --nom_dir {params.dir} \
          --cell_type {wildcards.cell_type} \
          --norm_method {params.norm_method} \
          --concat_out {output.cat_snps} \
          --summary_out {output.summary}  >> {log} 2>&1
        """

rule create_query:
    input:  eqtl = rules.cat_tensorqtl_nom_snps.output.cat_snps,
            snps = rules.cat_refs.output.bim,
            genes = config["susie"]["prep_susie_gene_meta"]["output"],
            frq = rules.cat_refs.output.frq,
            gene_lookup = config["tensorQTL"]["prep_tensorQTL_input"]["gene_lookup"]
    output: query = config["smr"]["create_query"]["query"],
            gene_lst = config["smr"]["create_query"]["gene_lst"]
    params: remove_mhc = config["smr"]["create_query"]["remove_mhc"],
            mhc_chr = config["mhc_region"]["chr"],
            mhc_start = config["mhc_region"]["start"],
            mhc_end = config["mhc_region"]["end"]
    singularity: config["containers"]["R"]
    resources: mem_mb = 20000
    message: "Create query file for SMR"
    benchmark: "reports/benchmarks/10smr.create_query_{cell_type}.benchmark.txt"
    log: config["smr"]["create_query"]["log"]
    script: "../scripts/smr_create_query.R"

rule create_besd:
    input:  query = rules.create_query.output.query,
            smr = rules.get_smr_binary.output
    output: config["smr"]["create_besd"]["besd"]
    params: prefix = config["smr"]["create_besd"]["prefix"],
            n_samples = config["smr_n_samples"]
    resources: mem_mb = 20000
    envmodules: "compiler/gnu/5/5.0"
    message: "Create besd file for SMR"
    benchmark: "reports/benchmarks/10smr.create_besd_{cell_type}.benchmark.txt"
    log: config["smr"]["create_besd"]["log"]
    shell:
        """
        {input.smr} --qfile {input.query} --add-n {params.n_samples} --make-besd --out {params.prefix} > {log} 2>&1
        """

rule format_gwas:
    input:  gwas = config["smr"]["format_gwas"]["gwas"],
            frq  = rules.cat_refs.output.frq
    output: config["smr"]["format_gwas"]["ma"]
    resources: mem_mb = 20000
    singularity: config["containers"]["R"]
    message: "Prep GWAS into .ma format for SMR"
    benchmark: "reports/benchmarks/10smr.format_gwas_{gwas}.benchmark.txt"
    log: config["smr"]["format_gwas"]["log"]
    script: "../scripts/smr_format_gwas.R"

rule smr:
    input:  bin = rules.get_smr_binary.output,
            gwas = rules.format_gwas.output,
            besd = rules.create_besd.output,
            ref_bed = rules.cat_refs.output.bed,
            ref_bim = rules.cat_refs.output.bim,
            ref_fam = rules.cat_refs.output.fam
    output: config["smr"]["smr"]["smr"]
    params: geno_prefix = config["smr"]["smr"]["geno_prefix"],
            besd_prefix = config["smr"]["create_besd"]["prefix"],
            out_prefix = config["smr"]["smr"]["smr_prefix"],
            peqtl_smr = config["peqtl_smr"]
    resources: mem_mb = 20000
    envmodules: "compiler/gnu/5/5.0"
    message: "Run SMR"
    benchmark: "reports/benchmarks/10smr.smr_{cell_type}_{gwas}.benchmark.txt"
    log:    config["smr"]["smr"]["log"]
    shell:  """
            {input.bin} --bfile {params.geno_prefix} \
                --gwas-summary {input.gwas} \
                --beqtl-summary {params.besd_prefix} \
                --out {params.out_prefix} \
                --peqtl-smr {params.peqtl_smr} > {log} 2>&1
            """

rule smr_report:
    input:  smr_files = expand(rules.smr.output, cell_type = config["cell_types"], gwas = config["gwas"]),
            gene_lookup = config["tensorQTL"]["prep_tensorQTL_input"]["gene_lookup"],
            rmd_script = "scripts/smr_report.Rmd"
    output: config["smr"]["smr_report"]["html"]
    params: cell_types = ",".join(config["cell_types"]),
            gwas = ",".join(config["gwas"]),
            in_dir = config["smr"]["smr_report"]["in_dir"],
            log_dir = config["smr"]["smr_report"]["log_dir"],
            bmark_dir = config["smr"]["smr_report"]["bmark_dir"],
            p_smr = config["p_smr"],
            p_heidi = config["p_heidi"],
            peqtl_smr = config["peqtl_smr"],
            smr_n_samples = config["smr_n_samples"]
    resources: mem_mb = 20000
    singularity: config["containers"]["r_eqtl"]
    message: "Generate SMR report"
    benchmark: "reports/benchmarks/10smr.smr_report.benchmark.txt"
    log: config["smr"]["smr_report"]["log"]
    shell:
        """
        Rscript -e "wd <- getwd(); rmarkdown::render('{input.rmd_script}', \
            output_file = file.path(wd, '{output}'), \
            knit_root_dir = wd, \
            params = list(cell_types = '{params.cell_types}', \
                          gwas = '{params.gwas}', \
                          in_dir = '{params.in_dir}', \
                          log_dir = '{params.log_dir}', \
                          gene_lookup = '{input.gene_lookup}', \
                          bmark_dir = '{params.bmark_dir}', \
                          p_smr = '{params.p_smr}', \
                          p_heidi = '{params.p_heidi}', \
                          peqtl_smr = '{params.peqtl_smr}', \
                          smr_n_samples = '{params.smr_n_samples}'))" > {log} 2>&1
        """
