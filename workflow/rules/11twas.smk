configfile: "../config/config.yaml"
#localrules: prep_exp_data

rule all:
    input: 
        "reports/11TWAS/11twas_weights_report.html"        

#def weights_input_function(wildcards):
#    """
#    Return the list of all .wgt.RDat files needed for a given cell_type
#    by reading the cell-type-specific permanent coordinate file.
#    """
#    import pandas as pd
#    
#    coord_file = config["twas"]["prep_exp_data"]["coord"].format(cell_type=wildcards.cell_type)

#    df = pd.read_csv(coord_file, sep=r'\s+', header=0)
#    # Your file has header: chr start end gene_id
#    genes = df["gene_id"].tolist()

#    return expand(
#        "../results/11TWAS/weights/{cell_type}/{gene}.wgt.RDat",
#        cell_type=wildcards.cell_type,
#        gene=genes
#    )

#def pos_input_function(wildcards):
#    """
#    Return list of successfully computed .wgt.RDat files for a given cell_type.
#    Only includes genes that actually produced a weight file (some may be empty).
#    """
#    import glob, os
#    weight_dir = f"../results/11TWAS/weights/{wildcards.cell_type}"
#    rdats = glob.glob(f"{weight_dir}/*.wgt.RDat")
#    # Filter out zero-byte files (genes that were skipped)
#    rdats = [f for f in rdats if os.path.getsize(f) > 0]
#    return rdats

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
    singularity: config["containers"]["r_eqtl"]
    message: "Creating plink ready expression data and gene coordinate file for FUSION weight calculation"
    benchmark: "reports/benchmarks/11twas.prep_exp_data_{cell_type}.benchmark.txt"    
    log:    config["twas"]["prep_exp_data"]["log"]
    script: "../scripts/twas_prep_expr_for_FUSION.R"

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

rule get_ldref_snplist:
    input:       expand(config["twas"]["get_ldref_snplist"]["input"], chr = range(1, 23))
    output:      config["twas"]["get_ldref_snplist"]["snp_lst"]
    envmodules: "plink/2.0"
    message:    "Restrict genotyped SNPs to LD reference SNPs"
    benchmark:  "reports/benchmarks/11twas.restrict_geno_to_ldref.benchmark.txt"
    log:        config["twas"]["get_ldref_snplist"]["log"]
"../results/00LOG/11TWAS/restrict_genotypes_to_ldref.log"
    shell:      """
                # Extract SNP IDs from LD reference .bim files
                cat {input} | cut -f 2 > {output} 2>&1 | tee {log}                
                """

rule prepare_covar:
    input:  
        lambda w: config["tensorQTL"]["split_covariates"]["output"].format(
            cell_type=w.cell_type, 
            norm_method="quantile",
            geno_pc=config['tensorQTL']['geno_pcs'], 
            exp_pc=config['exp_pc_map'][w.cell_type]
        )
    output: config["twas"]["prepare_covar"]["output"]
    message: "Prepare covariate gene expression data for FUSION"
    singularity: config["containers"]["r_eqtl"]
    benchmark: "reports/benchmarks/11twas.prepare_covariates.benchmark_{cell_type}.txt"
    log:  config["twas"]["prepare_covar"]["log"]
    script:   "../scripts/twas_prep_covar_for_FUSION.R"    

rule compute_weights:
    input:  bed = config["twas"]["convert_vcf"]["bed"],
            bim = config["twas"]["convert_vcf"]["bim"],
            fam = config["twas"]["convert_vcf"]["fam"],
            expr_file = rules.prep_exp_data.output.exp,
            coord_file = rules.prep_exp_data.output.coord,
            ldref_snps = rules.get_ldref_snplist.output,
            covar = rules.prepare_covar.output,
            gemma = config["twas"]["get_gemma"]["output"]
    output: touch(config["twas"]["compute_weights"]["output"])
    params: prefix_in = config["twas"]["convert_vcf"]["prefix"],
            outdir = config["twas"]["convert_vcf"]["outdir"]
    singularity: config["containers"]["twas"]
    message: "Running FUSION weights for ALL genes in {wildcards.cell_type}"
    benchmark: "reports/benchmarks/11twas.compute_all_weights_{cell_type}.benchmark.txt"
    log: config["twas"]["compute_weights"]["log"]
    shell:
        r"""
        echo "Starting weight computation for {wildcards.cell_type}" > {log}

        COORD={input.coord_file}
        EXPR={input.expr_file}
        LDREF={input.ldref_snps}
        COVAR={input.covar}
        GEMMA={input.gemma}
        GENOPREFIX={params.prefix_in}
        OUTDIR={params.outdir}

        mkdir -p $OUTDIR
        
        while read CHR START END GENE; do
            echo "Processing $GENE" >> {log}

            WORKDIR="../results/11TWAS/fusion_input/{wildcards.cell_type}/${{GENE}}"
            mkdir -p $WORKDIR

            CIS_START=$((START - 500000))
            CIS_END=$((END + 500000))
            [ $CIS_START -lt 0 ] && CIS_START=0

            # Extract cis SNPs
            plink2 --bfile $GENOPREFIX \
                   --chr $CHR \
                   --from-bp $CIS_START \
                   --to-bp $CIS_END \
                   --extract $LDREF \
                   --force-intersect \
                   --make-bed \
                   --out $WORKDIR/cis || true

            for ext in bed bim fam; do
                [[ -s "$WORKDIR/cis.$ext" ]] || touch "$WORKDIR/cis.$ext"
            done

            # Extract expression column
            COL=$(head -1 $EXPR | tr '\t' '\n' | grep -n -w $GENE | cut -d: -f1)
            cut -f1,2,$COL $EXPR > $WORKDIR/pheno.txt

            # Run FUSION if SNPs exist
            OUTPREFIX="$OUTDIR/${{GENE}}"
            if [ -s "$WORKDIR/cis.bed" ]; then
                Rscript ../resources/fusion/FUSION.compute_weights.R \
                    --bfile $WORKDIR/cis \
                    --pheno $WORKDIR/pheno.txt \
                    --covar $COVAR \
                    --hsq_p 0.01 \
                    --crossval 5 \
                    --tmp $OUTPREFIX.tmp \
                    --PATH_gemma $GEMMA \
                    --PATH_gcta ../resources/fusion/gcta_nr_robust \
                    --PATH_plink /apps/genomics/plink/1.9/el7/AVX512/intel-2018/serial/plink-1.9/usr/local/bin/plink \
                    --verbose 2 \
                    --models top1,lasso,enet \
                    --out $OUTPREFIX >> {log} 2>&1 || true
            fi

            [ -f "$OUTPREFIX.wgt.RDat" ] || touch "$OUTPREFIX.wgt.RDat"
            
            rm -rf $WORKDIR 
            if [ ! -d "$WORKDIR" ]; then
              echo "$WORKDIR removed" >> {log}
            else
              echo "WARNING: Failed to remove $WORKDIR" >> {log}
            fi 
            
        done < <(tail -n +2 $COORD)

        echo "All genes completed" >> {log}
        """

rule make_pos_file:
    input:  rules.compute_weights.output
    output: config["twas"]["make_pos_file"]["output"]
    params: coord_file = config["twas"]["prep_exp_data"]["coord"],
            cis_window = config["tensorQTL"]["window"]
    message:  "Generate FUSION .pos file for {wildcards.cell_type}"
    benchmark: "reports/benchmarks/11twas.make_pos_file_{cell_type}.benchmark.txt"
    log:      config["twas"]["make_pos_file"]["log"]
    script:   "../scripts/twas_make_fusion_pos_file.py"


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
    input:  pos = expand(rules.make_pos_file.output, cell_type = config['cell_types']),
            rmd_script = config["twas"]["twas_weights_report"]["script"]
    output: config["twas"]["twas_weights_report"]["output"]
    params: cell_types = ','.join(['\'{}\''.format(x) for x in config["cell_types"]]),
            weights_dir = config["twas"]["twas_weights_report"]["weights_dir"],
            bmark_dir = config["twas"]["twas_weights_report"]["bmark_dir"],
            out_file = config["twas"]["twas_weights_report"]["out_file"]
    singularity: config["containers"]["r_eqtl"]
    message:  "Generate TWAS weights report data for report"
    benchmark: "reports/benchmarks/11twas.twas_weights_summary.benchmark.txt"
    log: config["twas"]["twas_weights_report"]["log"]
    shell:
        """
        Rscript -e "rmarkdown::render('{input.rmd_script}', \
            output_file = '{params.output_file}', \
            params = list(cell_types = c({params.cell_types}), log_dir = '{params.log_dir}'))" > {log} 2>&1
        """

