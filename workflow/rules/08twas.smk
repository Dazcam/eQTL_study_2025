localrules: prep_exp_data

# Rule to download gemma binary for FUSION
rule get_gemma:
   output:  config["twas"]["get_gemma"]["output"] 
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
    singularity: config["containers"]["R"]
    log:    config["twas"]["prep_exp_data"]["log"]
    script: "../scripts/prep_expr_for_FUSION.R"

# Rule to convert VCF to PLINK
rule convert_vcf:
    input: config["twas"]["convert_vcf"]["input"]
    output: bed = config["twas"]["convert_vcf"]["bed"],
            bim = config["twas"]["convert_vcf"]["bim"],
            fam = config["twas"]["convert_vcf"]["fam"]
    params: config["twas"]["convert_vcf"]["prefix"]
    envmodules: "plink/2.0"
    log:    config["twas"]["convert_vcf"]["log"]
    shell:  """plink2 --vcf {input} --make-bed --out {params} > {log} 2>&1"""   

# Rule to restrict genotypes to LD reference SNPs
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
    log:
        "../results/00LOG/06TWAS/restrict_genotypes_to_ldref.log"
    shell:
        """
        # Extract SNP IDs from LD reference .bim files
        cat {input.ldref_bim} | cut -f 2 > {output.snp_list}
        # Restrict genotypes to LD reference SNPs
        plink2 --bfile {params.input} \
              --extract {output.snp_list} \
              --make-bed --out {params.output} > {log} 2>&1
        """

checkpoint get_gene_info:
    input:  config["twas"]["prep_exp_data"]["coord"]
    output: directory("../results/06TWAS/fusion_input/{cell_type}/gene_info")
    shell:
        """
        mkdir {output}
        tail -n +2 {input} | while IFS=$'\t' read -r line; do
        gene=$(echo "$line" | cut -f4)
        echo "$line" > "{output}/${{gene}}_gene_info.txt"
        done
        """

rule extract_cis_snps:
    input:
        geno_bed = config["twas"]["restrict_geno_to_ldref"]["bed"],
        geno_bim = config["twas"]["restrict_geno_to_ldref"]["bim"],
        geno_fam = config["twas"]["restrict_geno_to_ldref"]["fam"],
        gene_info = "../results/06TWAS/fusion_input/{cell_type}/gene_info/{gene_id}_gene_info.txt"
    output:
        cis_bed = "../results/06TWAS/fusion_input/{cell_type}/{gene_id}_cis.bed",
        cis_bim = "../results/06TWAS/fusion_input/{cell_type}/{gene_id}_cis.bim",
        cis_fam = "../results/06TWAS/fusion_input/{cell_type}/{gene_id}_cis.fam"
    params:
        prefix_in = config["twas"]["restrict_geno_to_ldref"]["prefix_out"],
        prefix_out = "../results/06TWAS/fusion_input/{cell_type}/{gene_id}_cis"
    envmodules: "plink/2.0"
    shell:
        """
        GENE_CHR=$(awk '{{print $1}}' {input.gene_info})
        GENE_START=$(awk '{{print $2}}' {input.gene_info})
        GENE_END=$(awk '{{print $3}}' {input.gene_info})
        CIS_START=$((GENE_START - 500000))
        CIS_END=$((GENE_END + 500000))

        if (( CIS_START < 0 )); then
          CIS_START=0
        fi

        plink2 --bfile {params.prefix_in} \
               --chr $GENE_CHR \
               --from-bp $CIS_START \
               --to-bp $CIS_END \
               --make-bed \
               --out {params.prefix_out} || true

        # Create empty stub files if PLINK didn't output anything
        for ext in bed bim fam; do
            f="{params.prefix_out}.$ext"
            if [ ! -s "$f" ]; then
                echo -n "" > "$f"
            fi
        done
        """

rule prepare_gene_pheno:
    input:  config["twas"]["prep_exp_data"]["exp"],
    output: "../results/06TWAS/fusion_input/{cell_type}/{gene_id}_pheno.txt"
    params: "{gene_id}"
    shell:
        """
        # Get the column index of the gene_id
        COL_IDX=$(head -n 1 {input} | tr '\t' '\n' | grep -n -w {params} | cut -f1 -d:)
        # Extract FID, IID, and the gene_id column
        cut -f 1,2,${{COL_IDX}} {input} > {output}
        """


rule compute_weights:
    input:
        geno_bed = "../results/06TWAS/fusion_input/{cell_type}/{gene_id}_cis.bed",
        geno_bim = "../results/06TWAS/fusion_input/{cell_type}/{gene_id}_cis.bim",
        geno_fam = "../results/06TWAS/fusion_input/{cell_type}/{gene_id}_cis.fam",
        gene_pheno = "../results/06TWAS/fusion_input/{cell_type}/{gene_id}_pheno.txt",
#        covar = "../results/03SCANPY/pseudobulk/{cell_type}_covariates.txt", # Need to get this into correct format and consider whether to include them or not
        gemma = config["twas"]["get_gemma"]["output"]
    output:
        "../results/06TWAS/weights/{cell_type}/{gene_id}.RDat"
    params:
        indir = "../results/06TWAS/fusion_input/{cell_type}/{gene_id}_cis",
        outdir = "../results/06TWAS/weights/{cell_type}"
    singularity: config["containers"]["twas"]
    log: "../results/00LOG/06TWAS/compute_weights_{cell_type}_{gene_id}.log"
    shell:
        r"""
        if [ -s {input.geno_bed} ]; then
          echo "Running compute_weights for {wildcards.gene_id}" >> {log}
          Rscript ../resources/fusion/FUSION.compute_weights.R \
            --bfile {params.indir} \
            --pheno {input.gene_pheno} \
            --hsq_p 0.01 \
            --crossval 5 \
            --PATH_plink /apps/genomics/plink/1.9/el7/AVX512/intel-2018/serial/plink-1.9/usr/local/bin/plink \
            --PATH_gcta ../resources/fusion/gcta_nr_robust \
            --PATH_gemma {input.gemma} \
            --out {output} >> {log} 2>&1 || true
          if [ ! -f {output} ]; then
            echo "{wildcards.gene_id} was skipped (likely due to low heritability)." >> {log}
            echo "{wildcards.gene_id} was skipped due to low heritability or other issues." > ../results/06TWAS/weights/{wildcards.cell_type}/{wildcards.gene_id}_skipped.txt
            touch {output}
          fi
        else
          echo "No SNPs for {wildcards.gene_id}, skipping..." >> {log}
          echo "{wildcards.gene_id} has no SNPs in the cis region." > ../results/06TWAS/weights/{wildcards.cell_type}/{wildcards.gene_id}_no_snps.txt
          touch {output}
        fi
        """

def aggregate_input(wildcards):
    checkpoint_output = checkpoints.get_gene_info.get(**wildcards).output[0]
    return expand(
        "../results/06TWAS/weights/{cell_type}/{gene_id}.RDat",
        cell_type = wildcards.cell_type,
        gene_id = glob_wildcards(os.path.join(checkpoint_output, "{gene_id}_gene_info.txt")).gene_id,
    )

rule aggregate_per_cell_type:
    input: aggregate_input
    output: "../results/06TWAS/{cell_type}/all_genes_done.txt"
    shell:
        """
        echo "Processed genes for {wildcards.cell_type}:" > {output}
        for file in {input}; do
            basename $file .RDat >> {output}
        done
        """

rule aggregate_all:
    input:  expand("../results/06TWAS/{cell_type}/all_genes_done.txt", cell_type = config['cell_types']),
    output: "../results/06TWAS/all_genes_in_all_cell_types_done.txt",
    shell:  "cat {input} > {output}"


rule twas_weights_summary:
    # Note diff paths for output and out_file; Rmarkdown needs outfile to be relative to Rmd file
    input:  agg_file = "../results/06TWAS/all_genes_in_all_cell_types_done.txt",
            rmd_file = "scripts/twas_weights_summary.Rmd"
    output: "reports/06TWAS/twas_weights_summary.html"
    params: cell_types = config["cell_types"], 
            log_dir = "../results/00LOG/06TWAS", 
            output_file = ../reports/06TWAS/twas_weights_summary.html",
            hsq_p = 0.01  # Example additional parameter
    singularity: config["containers"]["R"]
    log: "../results/00LOG/06TWAS/twas_weights_summary.log" 
    shell:
        """
        Rscript -e "rmarkdown::render('{input.rmd_file}', \
        output_file = '{params.output_file}')" > {log} 2>&1
        """
