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

PT_Q4_EXPORT_TARGETS = [
    f"../results/05TENSORQTL/prep_input/.pt_export_{traj}_Q4_bin{b}.done"
    for traj in config["trajectories"]
    for b in range(1, 5)
]

rule all:
    input:
        config["pseudotime"]["report"]["html"],
#        [
#            f"../results/05TENSORQTL/prep_input/.pt_export_{traj}_Q4_bin{bin}.done"
#            for traj in config["trajectories"]
#            for bin in range(1, 5)
#        ]
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

rule pseudotime_export_smr_twas_inputs:
    input:
        parquet_sentinel = lambda w: (
            "../results/17PSEUDOTIME/{trajectory}/tensorqtl/nominal/"
            "Q4_bin{bin}_genPC4_expPC{pc}/Q4_bin{bin}_nom.cis_qtl_pairs.1.parquet"
        ).format(trajectory=w.trajectory, bin=w.bin,
                  pc=config["pseudotime"]["bin_exp_pc_map"][w.trajectory][int(w.bin)]),
        perm = lambda w: (
            "../results/17PSEUDOTIME/{trajectory}/tensorqtl/perm/"
            "Q4_bin{bin}_genPC4_expPC{pc}/Q4_bin{bin}_perm.cis_qtl.txt.gz"
        ).format(trajectory=w.trajectory, bin=w.bin,
                  pc=config["pseudotime"]["bin_exp_pc_map"][w.trajectory][int(w.bin)]),
        bed = "../results/17PSEUDOTIME/{trajectory}/tensorqtl/prep_input/Q4_bin{bin}_quantile.bed",
        cov = lambda w: (
            "../results/17PSEUDOTIME/{trajectory}/tensorqtl/prep_input/"
            "Q4_bin{bin}_quantile_genPC4_expPC{pc}_split_covariates.txt"
        ).format(trajectory=w.trajectory, bin=w.bin,
                  pc=config["pseudotime"]["bin_exp_pc_map"][w.trajectory][int(w.bin)])
    output:
        touch("../results/05TENSORQTL/prep_input/.pt_export_{trajectory}_Q4_bin{bin}.done")
    log:
        "../results/00LOG/17PSEUDOTIME/{trajectory}/export_smr_twas_Q4_bin{bin}.log"
    benchmark:
        "reports/benchmarks/17pseudotime.export_smr_twas_inputs_{trajectory}_Q4_bin{bin}.txt"
    wildcard_constraints:
        trajectory = "|".join(config["trajectories"]),
        bin = "[1-4]"
    message: "Exporting {wildcards.trajectory} Q4 bin {wildcards.bin} inputs for SMR/TWAS/SuSiE"
    run:
        import glob, os, shutil

        pc = config["pseudotime"]["bin_exp_pc_map"][wildcards.trajectory][int(wildcards.bin)]
        pt_group = f"{wildcards.trajectory}-Q4-Bin{wildcards.bin}"

        src_dir = os.path.dirname(input.parquet_sentinel)
        dest_dir = f"../results/05TENSORQTL/tensorqtl_nom/{pt_group}_quantile_genPC_4_expPC_{pc}"
        os.makedirs(dest_dir, exist_ok=True)
        for src in glob.glob(os.path.join(src_dir, "Q4_bin*_nom.cis_qtl_pairs.*.parquet")):
            chrom_part = os.path.basename(src).split("cis_qtl_pairs.")[1]
            shutil.copy(src, os.path.join(dest_dir, f"{pt_group}_quantile_nom.cis_qtl_pairs.{chrom_part}"))

        perm_dest_dir = f"../results/05TENSORQTL/tensorqtl_perm/{pt_group}_quantile_genPC_4_expPC_{pc}"
        os.makedirs(perm_dest_dir, exist_ok=True)
        shutil.copy(input.perm, os.path.join(perm_dest_dir, f"{pt_group}_quantile_perm.cis_qtl.txt.gz"))

        bed_dest = f"../results/05TENSORQTL/prep_input/{pt_group}_quantile.bed"
        shutil.copy(input.bed, bed_dest)

        cov_dest = f"../results/05TENSORQTL/prep_input/{pt_group}_quantile_genPC_4_expPC_{pc}_split_covariates.txt"
        shutil.copy(input.cov, cov_dest)

#rule pseudotime_export_smr_twas_inputs:
#    input:
#        parquet_sentinel = lambda w: (
#            "../results/17PSEUDOTIME/{trajectory}/tensorqtl/nominal/"
#            "Q4_bin{bin}_genPC4_expPC{pc}/Q4_bin{bin}_nom.cis_qtl_pairs.1.parquet"
#        ).format(trajectory=w.trajectory, bin=w.bin,
#                  pc=config["pseudotime"]["cell_type_exp_pc_map"][w.trajectory]),
#        perm = lambda w: (
#            "../results/17PSEUDOTIME/{trajectory}/tensorqtl/perm/"
#            "Q4_bin{bin}_genPC4_expPC{pc}/Q4_bin{bin}_perm.cis_qtl.txt.gz"
#        ).format(trajectory=w.trajectory, bin=w.bin,
#                  pc=config["pseudotime"]["cell_type_exp_pc_map"][w.trajectory]),
#        bed = "../results/17PSEUDOTIME/{trajectory}/tensorqtl/prep_input/Q4_bin{bin}_quantile.bed",
#        cov = lambda w: (
#            "../results/17PSEUDOTIME/{trajectory}/tensorqtl/prep_input/"
#            "Q4_bin{bin}_quantile_genPC4_expPC{pc}_split_covariates.txt"
#        ).format(trajectory=w.trajectory, bin=w.bin,
#                  pc=config["pseudotime"]["cell_type_exp_pc_map"][w.trajectory])
#    output:
#        touch("../results/05TENSORQTL/prep_input/.pt_export_{trajectory}_Q4_bin{bin}.done")
#    log:
#        "../results/00LOG/17PSEUDOTIME/{trajectory}/export_smr_twas_Q4_bin{bin}.log"
#    benchmark:
#        "reports/benchmarks/17pseudotime.export_smr_twas_inputs_{trajectory}_Q4_bin{bin}.txt"
#    wildcard_constraints:
#        trajectory = "|".join(config["trajectories"]),
#        bin = "[1-4]"
#    message: "Exporting {wildcards.trajectory} Q4 bin {wildcards.bin} inputs for SMR/TWAS/SuSiE"
#    run:
#        import glob, os, shutil

#        pc = config["pseudotime"]["cell_type_exp_pc_map"][wildcards.trajectory]
#        pt_group = f"{wildcards.trajectory}-Q4-Bin{wildcards.bin}"

#        # Nominal parquets
#        src_dir = os.path.dirname(input.parquet_sentinel)
#        dest_dir = f"../results/05TENSORQTL/tensorqtl_nom/{pt_group}_quantile_genPC_4_expPC_{pc}"
#        os.makedirs(dest_dir, exist_ok=True)
#        for src in glob.glob(os.path.join(src_dir, "Q4_bin*_nom.cis_qtl_pairs.*.parquet")):
#            chrom_part = os.path.basename(src).split("cis_qtl_pairs.")[1]
#            shutil.copy(src, os.path.join(dest_dir, f"{pt_group}_quantile_nom.cis_qtl_pairs.{chrom_part}"))

#        # Permutation output (single file — needed by SuSiE's get_sig_eGenes/run_susie)
#        perm_dest_dir = f"../results/05TENSORQTL/tensorqtl_perm/{pt_group}_quantile_genPC_4_expPC_{pc}"
#        os.makedirs(perm_dest_dir, exist_ok=True)
#        shutil.copy(input.perm, os.path.join(perm_dest_dir, f"{pt_group}_quantile_perm.cis_qtl.txt.gz"))

#        # Expression bed
#        bed_dest = f"../results/05TENSORQTL/prep_input/{pt_group}_quantile.bed"
#        shutil.copy(input.bed, bed_dest)

#        # Split covariates
#        cov_dest = f"../results/05TENSORQTL/prep_input/{pt_group}_quantile_genPC_4_expPC_{pc}_split_covariates.txt"
#        shutil.copy(input.cov, cov_dest)

rule pseudotime_report:
    input:  perm           = get_pseudotime_tensorqtl(),
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

#rule pseudotime_report:
#    input:  perm           = get_pseudotime_tensorqtl(),
#            dynamic        = expand(
#                                config["pseudotime"]["dynamic_eqtl"]["candidates_tsv"],
#                                trajectory = config["trajectories"]
#                             ),
#            rmd_script     = config["pseudotime"]["report"]["rmd_script"]
#    output: config["pseudotime"]["report"]["html"]
#    params: in_dir              = config["pseudotime"]["report"]["in_dir"],
#            cell_type_perm_dir  = config["pseudotime"]["report"]["cell_type_perm_dir"],
#            output_file         = config["pseudotime"]["report"]["output_file"]
#    singularity: config["containers"]["r_eqtl"]
#    resources: threads = 1, mem_mb = 8000, time = "1:00:00"
#    log:    config["pseudotime"]["report"]["log"]
#    message: "Generating pseudotime eQTL report"
#    shell:
#        """
#        Rscript -e "rmarkdown::render('{input.rmd_script}', \
#            output_file = '{params.output_file}', \
#            params = list(
#                in_dir             = '{params.in_dir}',
#                cell_type_perm_dir = '{params.cell_type_perm_dir}'
#            ))" > {log} 2>&1
#        """

