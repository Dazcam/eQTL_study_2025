localrules: prep_exp_data

def weights_input_function(wildcards):
    """
    Return the list of all .wgt.RDat files needed for a given cell_type
    by reading the cell-type-specific permanent coordinate file.
    """
    import pandas as pd
    
    coord_file = config["twas"]["prep_exp_data"]["coord"] \
                 .format(cell_type=wildcards.cell_type)
    # → becomes ../results/11TWAS/fusion_input/MyCellType_gene_coord.txt

    df = pd.read_csv(coord_file, sep=r'\s+', header=0)
    # Your file has header: chr start end gene_id
    genes = df["gene_id"].tolist()

    return expand(
        "../results/11TWAS/weights/{cell_type}/{gene}.wgt.RDat",
        cell_type=wildcards.cell_type,
        gene=genes
    )

def pos_input_function(wildcards):
    """
    Return list of successfully computed .wgt.RDat files for a given cell_type.
    Only includes genes that actually produced a weight file (some may be empty).
    """
    import glob, os
    weight_dir = f"../results/11TWAS/weights/{wildcards.cell_type}"
    rdats = glob.glob(f"{weight_dir}/*.wgt.RDat")
    # Filter out zero-byte files (genes that were skipped)
    rdats = [f for f in rdats if os.path.getsize(f) > 0]
    return rdats

rule get_gemma:
   output:  config["twas"]["get_gemma"]["output"] 
   message: "Download gemma binary for FUSION"
   benchmark: "reports/benchmarks/11twas.get_gemma.benchmark.txt"   
   shell:  
           """
           wget https://github.com/genetics-statistics/GEMMA/releases/download/v0.98.5/gemma-0.98.5-linux-static-AMD64.gz
           gunzip gemma-0.98.5-linux-static-AMD64.gz
           mv gemma-0.98.5-linux-static-AMD64 {output}
           chmod +x {output}
           """

rule prep_exp_data:
    input:  config["twas"]["prep_exp_data"]["input"]
    output: exp = config["twas"]["prep_exp_data"]["exp"],
            coord = config["twas"]["prep_exp_data"]["coord"]
    params: run_test = config['twas_run_test']
    singularity: config["containers"]["R"]
    message: "Creating plink ready expression data and gene coordinate file for FUSION weight calculation"
    benchmark: "reports/benchmarks/11twas.prep_exp_data_{cell_type}.benchmark.txt"    
    log:    config["twas"]["prep_exp_data"]["log"]
    script: "../scripts/prep_expr_for_twas.R"

# Rule to convert VCF to PLINK
rule convert_vcf:
    input: config["twas"]["convert_vcf"]["input"]
    output: bed = config["twas"]["convert_vcf"]["bed"],
            bim = config["twas"]["convert_vcf"]["bim"],
            fam = config["twas"]["convert_vcf"]["fam"]
    params: config["twas"]["convert_vcf"]["prefix"]
    envmodules: "plink/2.0"
    message: "Converting genotypes from VCF to PLINK2 format for FUSION weight calculation"
    benchmark: "reports/benchmarks/11twas.convert_vcf.benchmark.txt"
    log:    config["twas"]["convert_vcf"]["log"]
    shell:  """plink2 --vcf {input} --make-bed --out {params} > {log} 2>&1"""   

rule restrict_geno_to_ldref:
    input:  bed = config["twas"]["convert_vcf"]["bed"],
            bim = config["twas"]["convert_vcf"]["bim"],
            fam = config["twas"]["convert_vcf"]["fam"],
            ldref_bim = expand("../resources/ldsr/ldsr_hg38_refs/plink_files/1000G.EUR.hg38.{chr}.bim", chr = range(1, 23))
    output: bed = config["twas"]["restrict_geno_to_ldref"]["bed"],
            bim = config["twas"]["restrict_geno_to_ldref"]["bim"],
            fam = config["twas"]["restrict_geno_to_ldref"]["fam"],
            snp_list = config["twas"]["restrict_geno_to_ldref"]["snp_lst"]
    params: input = config["twas"]["restrict_geno_to_ldref"]["prefix_in"],
            output = config["twas"]["restrict_geno_to_ldref"]["prefix_out"]
    envmodules: "plink/2.0"
    message: "Restrict genotyped SNPs to LD reference SNPs"
    benchmark: "reports/benchmarks/11twas.restrict_geno_to_ldref.benchmark.txt"
    log:    "../results/00LOG/11TWAS/restrict_genotypes_to_ldref.log"
    shell:  """
            # Extract SNP IDs from LD reference .bim files
            cat {input.ldref_bim} | cut -f 2 > {output.snp_list}
            # Restrict genotypes to LD reference SNPs
            plink2 --bfile {params.input} \
              --extract {output.snp_list} \
              --make-bed --out {params.output} > {log} 2>&1
            """

rule extract_cis_snps:
    input:  geno_bed = config["twas"]["restrict_geno_to_ldref"]["bed"],
            geno_bim = config["twas"]["restrict_geno_to_ldref"]["bim"],
            geno_fam = config["twas"]["restrict_geno_to_ldref"]["fam"],
            coord_file = config["twas"]["prep_exp_data"]["coord"]   # ← new input
    output: cis_bed = temp("../results/11TWAS/fusion_input/{cell_type}/{gene_id}_cis.bed"),
            cis_bim = temp("../results/11TWAS/fusion_input/{cell_type}/{gene_id}_cis.bim"),
            cis_fam = temp("../results/11TWAS/fusion_input/{cell_type}/{gene_id}_cis.fam")
    params: prefix_in = config["twas"]["restrict_geno_to_ldref"]["prefix_out"],
            prefix_out = "../results/11TWAS/fusion_input/{cell_type}/{gene_id}_cis"
    envmodules: "plink/2.0"
    message:  "Extract cis-SNPs within cis-window for {wildcards.gene_id}"
    benchmark: "reports/benchmarks/11twas.extract_cis_snps_{cell_type}_{gene_id}.benchmark.txt"
    shell:
        """
        # Look up chr, start, end for this gene_id directly from the coord file
        LINE=$(grep -w {wildcards.gene_id} {input.coord_file})
        GENE_CHR=$(echo "$LINE" | awk '{{print $1}}')
        GENE_START=$(echo "$LINE" | awk '{{print $2}}')
        GENE_END=$(echo "$LINE" | awk '{{print $3}}')

        CIS_START=$((GENE_START - 500000))
        CIS_END=$((GENE_END + 500000))
        [[ $CIS_START -lt 0 ]] && CIS_START=0

        plink2 --bfile {params.prefix_in} \
               --chr $GENE_CHR \
               --from-bp $CIS_START \
               --to-bp $CIS_END \
               --make-bed \
               --out {params.prefix_out} || true

        # Create empty stubs if no SNPs
        for ext in bed bim fam; do
            [[ -s "{params.prefix_out}.$ext" ]] || touch "{params.prefix_out}.$ext"
        done
        """

rule prepare_gene_pheno:
    input:  config["twas"]["prep_exp_data"]["exp"],
    output: temp("../results/11TWAS/fusion_input/{cell_type}/{gene_id}_pheno.txt")
    params: "{gene_id}"
    message:  "Prepare individual gene expression data for FUSION"
    benchmark: "reports/benchmarks/11twas.prepare_gene_pheno.benchmark_{cell_type}_{gene_id}.txt"
    shell:
        """
        # Get the column index of the gene_id
        COL_IDX=$(head -n 1 {input} | tr '\t' '\n' | grep -n -w {params} | cut -f1 -d:)
        # Extract FID, IID, and the gene_id column
        cut -f 1,2,${{COL_IDX}} {input} > {output}
        """

rule prepare_covar:
    input:  config["tensorQTL"]["split_covariates"]["output"].format(cell_type="{cell_type}", norm_method="quantile", geno_pc=4, exp_pc=40)
    output: "../results/11TWAS/fusion_input/{cell_type}_covariates.txt"
    message: "Prepare covariate gene expression data for FUSION"
    singularity: config["containers"]["r_eqtl"]
    benchmark: "reports/benchmarks/11twas.prepare_covariates.benchmark_{cell_type}.txt"
    log:      "../results/00LOG/11TWAS/prepare_covar_{cell_type}.log"
    script:   "../scripts/prep_covar_for_FUSION.R"    

rule compute_weights:
    # Need to fix blup and bslmm in GEMMA. Running on 3 models only for now
    input:  geno_bed = "../results/11TWAS/fusion_input/{cell_type}/{gene_id}_cis.bed",
            geno_bim = "../results/11TWAS/fusion_input/{cell_type}/{gene_id}_cis.bim",
            geno_fam = "../results/11TWAS/fusion_input/{cell_type}/{gene_id}_cis.fam",
            gene_pheno = "../results/11TWAS/fusion_input/{cell_type}/{gene_id}_pheno.txt",
            covar = "../results/11TWAS/fusion_input/{cell_type}_covariates.txt",
            gemma = config["twas"]["get_gemma"]["output"]
    output: "../results/11TWAS/weights/{cell_type}/{gene_id}.wgt.RDat" # Docs say output is .RDat, which is not the case
    params: indir = "../results/11TWAS/fusion_input/{cell_type}/{gene_id}_cis",
            outdir = "../results/11TWAS/weights/{cell_type}",
            out_rdat = "../results/11TWAS/weights/{cell_type}/{gene_id}"
    singularity: config["containers"]["twas"]
    message:  "Compute expression weights with FUSION"
    benchmark: "reports/benchmarks/11twas.compute_weights_{cell_type}_{gene_id}.benchmark.txt"
    log:    "../results/00LOG/11TWAS/compute_weights_{cell_type}_{gene_id}.log"
    shell:
            r"""
            if [ -s {input.geno_bed} ]; then
              echo "Running compute_weights for {wildcards.gene_id}" >> {log}
              Rscript ../resources/fusion/FUSION.compute_weights.R \
                --bfile {params.indir} \
                --pheno {input.gene_pheno} \
                --hsq_p 0.01 \
                --crossval 5 \
                --tmp {params.outdir}/{wildcards.gene_id}.tmp \
                --PATH_plink /apps/genomics/plink/1.9/el7/AVX512/intel-2018/serial/plink-1.9/usr/local/bin/plink \
                --PATH_gcta ../resources/fusion/gcta_nr_robust \
                --PATH_gemma {input.gemma} \
                --models top1,lasso,enet \
                --verbose 2 \
                --out {params.out_rdat} >> {log} 2>&1 || true
              if [ ! -f {output} ]; then
                echo "{wildcards.gene_id} was skipped (likely due to low heritability)." >> {log}
                touch {output}
              fi
           else
             echo "No SNPs for {wildcards.gene_id}, skipping..." >> {log}
             touch {output}
           fi
           """

rule aggregate_per_cell_type:
    input: weights_input_function
    output: touch("../results/11TWAS/{cell_type}/all_genes_done.txt")

rule aggregate_all:
    input:  expand("../results/11TWAS/{cell_type}/all_genes_done.txt", cell_type = config['cell_types']),
    output: "../results/11TWAS/all_genes_in_all_cell_types_done.txt",
    message:  "Aggregate all (cell- and gene-specific) data after checkpoint"
    benchmark: "reports/benchmarks/11twas.aggregate_all.benchmark.txt"
    shell:  "cat {input} > {output}"

rule make_pos_file:
    input:  pos_input_function
    output: "../results/11TWAS/weights/{cell_type}/{cell_type}.pos"
    params: coord_file = config["twas"]["prep_exp_data"]["coord"],
            cis_window = config["tensorQTL"]["window"]
    message:  "Generate FUSION .pos file for {wildcards.cell_type}"
    benchmark: "reports/benchmarks/11twas.make_pos_file_{cell_type}.benchmark.txt"
    log:      "../results/00LOG/11TWAS/make_pos_file_{cell_type}.log"
    script:   "../scripts/make_fusion_pos_file.py"

# Not got base TWAS running yet, using CTWAS anyway

#rule clean_sumstats:
#    # FUSION wants LDSR munge_sumstats.py format, but can't handle blank rows so need to omit these first
#    input:    "../results/07PREP-GWAS/{gwas}_hg38_ldsr_ready.sumstats.gz"
#    output:   "../results/11TWAS/gwas/{gwas}_hg38_ldsr_ready.sumstats.tsv"
#    message: "Removing blank rows from sumstats file"
#    benchmark: "reports/benchmarks/11twas.clean_sumstats_{gwas}.benchmark.txt"
#    log:    "../results/00LOG/11TWAS/clean_sumstats_{gwas}.log"
#    shell:    """
#              zcat {input} | awk 'NR==1 || (NF >= 5 && $1 ~ /^rs/)' > {output} 2>&1 | tee {log}
#              """

#rule run_twas:
#    input:  pos = "../results/11TWAS/weights/{cell_type}/{cell_type}.pos",
#            sumstats = "../results/11TWAS/gwas/{gwas}_hg38_ldsr_ready.sumstats.tsv"
#    output: "../results/11TWAS/associations/{cell_type}/{cell_type}.{gwas}.chr{chr}.twas"
#    params: ref_ld = "../resources/ldsr/ldsr_hg38_refs/plink_files/1000G.EUR.hg38.",
#            weights_dir = "../results/11TWAS/weights/{cell_type}/"
#    resources: threads = 5, mem_mb = 20000, time="5:00:00"
#    singularity: config["containers"]["twas"]
#    log:    "../results/00LOG/11TWAS/run_twas_{cell_type}_{gwas}_chr{chr}.log"
#    benchmark: "reports/benchmarks/11twas.run_twas_{cell_type}_{gwas}_chr{chr}.benchmark.txt"
#    shell:  """
#            Rscript ../resources/fusion/FUSION.assoc_test.R \
#              --sumstats {input.sumstats} \
#              --weights {input.pos} \
#              --weights_dir {params.weights_dir} \
#              --ref_ld_chr {params.ref_ld} \
#              --chr {wildcards.chr} \
#              --out {output} > {log} 2>&1
#            """

#rule merge_twas:
#    input:  expand("../results/11TWAS/associations/{cell_type}/{cell_type}.{gwas}.chr{chr}.twas", cell_type = config['cell_types'], gwas = config['gwas'], chr=range(1, 23))
#    output: "../results/11TWAS/associations/{cell_type}/{cell_type}.{gwas}.twas"
#    shell:  """
#            head -n 1 {input[0]} > {output}
#            for f in {input}; do tail -n +2 "$f"; done >> {output}
#            """

rule twas_weights_report:
    # Note diff paths for output and out_file; Rmarkdown needs outfile to be relative to Rmd file
    input:  weights_done = "../results/11TWAS/all_genes_in_all_cell_types_done.txt",
#            twas_results = expand("../results/11TWAS/associations/{cell_type}/{cell_type}.{gwas}.twas", cell_type = config['cell_types'], gwas = config['gwas']),
            rmd_script = "scripts/twas_weights_summary.Rmd"
    output: "reports/11TWAS/11twas_weights_report.html"
    params: cell_types = ','.join(['\'{}\''.format(x) for x in config["cell_types"]]),
            log_dir = "../results/00LOG/11TWAS/", 
            output_file = "../reports/11TWAS/11twas_weights_report.html",
    singularity: config["containers"]["R"]
    message:  "Generate TWAS weights report data for report"
    benchmark: "reports/benchmarks/11twas.twas_weights_summary.benchmark.txt"
    log: "../results/00LOG/11TWAS/twas_weights_summary.log" 
    shell:
        """
        Rscript -e "rmarkdown::render('{input.rmd_script}', \
            output_file = '{params.output_file}', \
            params = list(cell_types = c({params.cell_types}), log_dir = '{params.log_dir}'))" > {log} 2>&1
        """

