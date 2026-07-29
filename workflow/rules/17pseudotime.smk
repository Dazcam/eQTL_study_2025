configfile: "../config/config.yaml"

def get_pseudotime_tensorqtl(wildcards=None):
    files = []
    for n in [4, 5, 6]:
        files.extend(
            expand(
                config["pseudotime"]["tensorqtl_perm"]["output"],
                trajectory=config["trajectories"],
                quantile_n=n,
                bin=range(1, n + 1),
                exp_pc=config["pseudotime"]["exp_pcs"]
            )
        )
        files.extend(
            expand(
                config["pseudotime"]["tensorqtl_nom"]["output"],
                trajectory=config["trajectories"],
                quantile_n=n,
                bin=range(1, n + 1),
                exp_pc=config["pseudotime"]["exp_pcs"]
            )
        )
    return files

rule all:
    input:
        config["pseudotime"]["report"]["html"]
#        config["pseudotime"]["palantir_all"]["output_h5ad"],        
#        expand(config["pseudotime"]["dynamic_eqtl"]["candidates_tsv"], trajectory = config["trajectories"])
#        get_pseudotime_tensorqtl()
#       expand(config["pseudotime"]["palantair"]["h5ad"], trajectory = config['trajectories']),
#       expand(config["pseudotime"]["pseudobulk"]["sentinel"], trajectory = config['trajectories']),

rule pseudotime_palantir_all:
    input:
        h5ad = config["pseudotime"]["palantir_all"]["input_h5ad"]
    output:
        h5ad = config["pseudotime"]["palantir_all"]["output_h5ad"]
    params:
        exclude_types = config["pseudotime"]["palantir_all"]["exclude_types"],
        root_gene     = config["pseudotime"]["palantir_all"]["root_gene"],
        terminals     = config["pseudotime"]["palantir_all"]["terminals"],
        marker_genes  = config["pseudotime"]["palantir_all"]["marker_genes"],
        n_neighbors   = config["pseudotime"]["palantir_all"]["n_neighbors"],
        n_pcs         = config["pseudotime"]["palantir_all"]["n_pcs"],
        n_dcs         = config["pseudotime"]["palantir_all"]["n_dcs"],
        num_waypoints = config["pseudotime"]["palantir_all"]["num_waypoints"]
    conda:    config["pseudotime"]["env"]
    resources: threads = 8, mem_mb = 200000, time = "02:00:00"
    threads:  8
    log:      config["pseudotime"]["palantir_all"]["log"]
    message:  "Running Palantir on full filtered atlas for fate probability estimation"
    script:   config["pseudotime"]["palantir_all"]["script"]

rule subset_trajectory:
    input:  h5ad   = config["pseudotime"]["subset"]["input_h5ad"]
    output: h5ad   = config["pseudotime"]["subset"]["h5ad"]
    conda:  config["pseudotime"]["env"]
    resources: threads = 2, mem_mb = 180000, time = "00:30:00"
    threads: 2
    log:    config["pseudotime"]["subset"]["log"]
    message: "Subsetting trajectory: {wildcards.trajectory}"
    script: config["pseudotime"]["subset"]["script"]

rule palantair_trajectory:
    input:  h5ad   = config["pseudotime"]["subset"]["h5ad"]
    output: h5ad   = config["pseudotime"]["palantair"]["h5ad"]
    conda:  config["pseudotime"]["env"]
    resources: threads = 4, mem_mb = 96000, time = "01:00:00"
    threads: 4
    log:    config["pseudotime"]["palantair"]["log"]
    message: "Running Palantir pseudotime for trajectory: {wildcards.trajectory}"
    params:
        n_neighbors   = config["pseudotime"]["palantair"].get("n_neighbors", 10),
        n_pcs         = config["pseudotime"]["palantair"].get("n_pcs", 30),
        n_dcs         = config["pseudotime"]["palantair"].get("n_dcs", 15),
        num_waypoints = config["pseudotime"]["palantair"].get("num_waypoints", 1000),
        root_gene     = config["pseudotime"]["palantair"].get("root_gene", "PAX6")
    script: config["pseudotime"]["palantair"]["script"]

rule pseudobulk_trajectory:
    input:  h5ad     = config["pseudotime"]["palantair"]["h5ad"],
            geno     = config["pseudotime"]["pseudobulk"]["geno_sample_file"]
    output: sentinel = config["pseudotime"]["pseudobulk"]["sentinel"]
    conda:  config["pseudotime"]["env"]
    resources: threads = 2, mem_mb = 90000, time = "00:15:00"
    threads: 2
    log:    config["pseudotime"]["pseudobulk"]["log"]
    message: "Generating pseudobulk quantile bins for trajectory: {wildcards.trajectory}"
    script: config["pseudotime"]["pseudobulk"]["script"]

rule pseudotime_prep_tensorqtl:
    input:  pseudobulk  = config["pseudotime"]["pseudobulk"]["sentinel"],
            cov_file    = config["pseudotime"]["prep_tensorqtl"]["cov_file"],
            sex_file    = config["pseudotime"]["prep_tensorqtl"]["sex_file"],
            gene_lookup = config["pseudotime"]["prep_tensorqtl"]["gene_lookup"]
    output: exp_out     = config["pseudotime"]["prep_tensorqtl"]["exp_out"],
            cov_out     = config["pseudotime"]["prep_tensorqtl"]["cov_out"]
    params: pseudoblk_dir = config["pseudotime"]["prep_tensorqtl"]["pseudoblk_dir"],
            out_dir       = config["pseudotime"]["prep_tensorqtl"]["out_dir"]
    singularity: config["containers"]["r_eqtl"]
    resources: threads = 1, mem_mb = 6000, time = "2:00:00"
    log:    config["pseudotime"]["prep_tensorqtl"]["log"]
    message: "Prepping tensorQTL inputs for {wildcards.trajectory} Q{wildcards.quantile_n} bin {wildcards.bin}"
    script: config["pseudotime"]["prep_tensorqtl"]["script"]

rule pseudotime_split_covariates:
    input:  config["pseudotime"]["split_covariates"]["input"]
    output: config["pseudotime"]["split_covariates"]["output"]
    singularity: config["containers"]["r_eqtl"]
    resources: threads = 1, mem_mb = 6000, time = "1:00:00"
    log:    config["pseudotime"]["split_covariates"]["log"]
    message: "Splitting covariates for {wildcards.trajectory} Q{wildcards.quantile_n} bin {wildcards.bin} expPC {wildcards.exp_pc}"
    script: config["pseudotime"]["split_covariates"]["script"]

rule pseudotime_zip_bed:
    input:  config["pseudotime"]["zip_bed"]["input"]
    output: config["pseudotime"]["zip_bed"]["output"]
    envmodules: "HTSlib/1.21-GCC-13.3.0"
    resources: threads = 1, mem_mb = 2000, time = "15:00"
    log:    config["pseudotime"]["zip_bed"]["log"]
    message: "bgzip and tabix BED for {wildcards.trajectory} Q{wildcards.quantile_n} bin {wildcards.bin}"
    shell:
            """
            bgzip -c {input} > {output}
            tabix -p bed {output}
            """

rule pseudotime_tensorqtl_nom:
    input:  genotypes  = config["tensorQTL"]["convert_genotypes"]["output"],
            counts     = config["pseudotime"]["zip_bed"]["output"],
            covariates = config["pseudotime"]["split_covariates"]["output"]
    output: config["pseudotime"]["tensorqtl_nom"]["output"]
    params: prefix_in  = config["tensorQTL"]["tensorqtl_nom"]["prefix_in"],
            prefix_out = config["pseudotime"]["tensorqtl_nom"]["prefix_out"],
            window     = config["tensorQTL"]["window"]
    singularity: config["containers"]["tensorqtl"]
    resources: threads = 4, mem_mb = 20000, time = "30:00"    
    threads: 4
    log:    config["pseudotime"]["tensorqtl_nom"]["log"]
    benchmark: config["pseudotime"]["tensorqtl_nom"]["benchmark"]
    message: "tensorQTL nominal for {wildcards.trajectory} Q{wildcards.quantile_n} bin {wildcards.bin} expPC {wildcards.exp_pc}"
    shell:
            """
            python3 -m tensorqtl {params.prefix_in} {input.counts} {params.prefix_out} \
               --covariates {input.covariates} \
               --window {params.window} \
               --mode cis_nominal >> {log} 2>&1
            """

rule pseudotime_extract_dosages:
    input:
        dosage    = config["pseudotime"]["dynamic_eqtl"]["dosage_file"],
        dosage_idx= config["pseudotime"]["dynamic_eqtl"]["dosage_index"],
        pvar      = config["pseudotime"]["extract_dosages"]["pvar_file"],
        cell_type_perm = lambda wildcards: expand(
            config["tensorQTL"]["tensorqtl_perm"]["output"],
            cell_type   = wildcards.trajectory.split("-to-")[1],
            norm_method = "quantile",
            geno_pc     = 4,
            exp_pc      = config["pseudotime"]["cell_type_exp_pc_map"][wildcards.trajectory]
        )
    output:
        dosage_sub = config["pseudotime"]["extract_dosages"]["output"]
    envmodules: "HTSlib/1.21-GCC-13.3.0"
    resources:  threads = 1, mem_mb = 32000, time = "1:00:00"
    log:        config["pseudotime"]["extract_dosages"]["log"]
    message:    "Extracting genotype dosages for trajectory: {wildcards.trajectory}"
    shell:
        """
        # Extract lead eSNP rsIDs from cell type perm file
        zcat {input.cell_type_perm} | \
            awk -F'\t' 'NR>1 && $NF+0<0.05 {{print $7}}' \
            > {output.dosage_sub}.snps.tmp

        # Get CHROM POS for lead SNPs from pvar
        grep -v '^##' {input.pvar} | \
            awk -F'\t' 'NR==FNR{{snps[$1]=1; next}} FNR>1 && ($3 in snps){{print $1"\t"$2}}' \
            {output.dosage_sub}.snps.tmp - \
            > {output.dosage_sub}.pos.tmp

        # Extract dosages in single pass
        zcat {input.dosage} | \
            awk -F'\t' 'NR==FNR{{pos[$1"\t"$2]=1; next}} FNR==1 || ($1"\t"$2 in pos)' \
            {output.dosage_sub}.pos.tmp - \
            > {output.dosage_sub} 2>> {log}

        # Clean up temp files
        rm {output.dosage_sub}.snps.tmp {output.dosage_sub}.pos.tmp
        """


rule pseudotime_dynamic_eqtl:
    input:
        perm = expand(
            config["pseudotime"]["tensorqtl_perm"]["output"],
            trajectory = "{trajectory}",
            quantile_n = 6,
            bin        = range(1, 7),
            exp_pc     = 30
        ),
        pseudobulk  = config["pseudotime"]["pseudobulk"]["sentinel"],
        dosage_sub  = config["pseudotime"]["extract_dosages"]["output"],
        pvar        = config["pseudotime"]["extract_dosages"]["pvar_file"]
    output:
        candidates_tsv = config["pseudotime"]["dynamic_eqtl"]["candidates_tsv"]
    params:
        perm_dir           = lambda wildcards: f"../results/17PSEUDOTIME/{wildcards.trajectory}/tensorqtl/perm",
        pseudobulk_dir     = lambda wildcards: f"../results/17PSEUDOTIME/{wildcards.trajectory}/pseudobulk",
        cov_dir            = lambda wildcards: f"../results/17PSEUDOTIME/{wildcards.trajectory}/tensorqtl/prep_input",
        results_rds        = lambda wildcards: f"../results/17PSEUDOTIME/{wildcards.trajectory}/dynamic_eqtl/dynamic_eqtl_results.rds",
        results_tsv        = lambda wildcards: f"../results/17PSEUDOTIME/{wildcards.trajectory}/dynamic_eqtl/dynamic_eqtl_significant.tsv",
        cell_type_perm_dir = "../results/05TENSORQTL/tensorqtl_perm"
    singularity: config["containers"]["r_eqtl"]
    resources:   threads = 1, mem_mb = 10000, time = "10:00"
    threads:     8
    log:         config["pseudotime"]["dynamic_eqtl"]["log"]
    message:     "Running dynamic eQTL LMM for trajectory: {wildcards.trajectory}"
    script:      config["pseudotime"]["dynamic_eqtl"]["script"]


rule pseudotime_report:
    input:  perm           = get_pseudotime_tensorqtl(),
            dynamic        = expand(
                                config["pseudotime"]["dynamic_eqtl"]["candidates_tsv"],
                                trajectory = config["trajectories"]
                             ),
            rmd_script     = config["pseudotime"]["report"]["rmd_script"]
    output: config["pseudotime"]["report"]["html"]
    params: in_dir              = config["pseudotime"]["report"]["in_dir"],
            cell_type_perm_dir  = config["pseudotime"]["report"]["cell_type_perm_dir"],
            output_file         = config["pseudotime"]["report"]["output_file"]
    singularity: config["containers"]["r_eqtl"]
    resources: threads = 1, mem_mb = 8000, time = "1:00:00"
    log:    config["pseudotime"]["report"]["log"]
    message: "Generating pseudotime eQTL report"
    shell:
        """
        Rscript -e "rmarkdown::render('{input.rmd_script}', \
            output_file = '{params.output_file}', \
            params = list(
                in_dir             = '{params.in_dir}',
                cell_type_perm_dir = '{params.cell_type_perm_dir}'
            ))" > {log} 2>&1
        """

