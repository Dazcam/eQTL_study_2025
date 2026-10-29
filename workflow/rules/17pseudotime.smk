configfile: "../config/config.yaml"

# Pseudotime quantile schemes: Q4, Q5 and Q6 equal-cell-number bins per trajectory.
# Single source for the pseudobulk script (params), its declared CSV outputs and the
# report inputs below. The pseudobulk QC figure has one panel per Q value (3 panels).
PT_Q_VALUES = [4, 5, 6]

# Genotype PCs are fixed for pseudotime; same value as the tensorQTL pipeline
PT_GENO_PC = config["tensorQTL"]["geno_pcs"]

# Every pseudobulk CSV written per trajectory ({trajectory} left as a wildcard)
PT_PSEUDOBULK_CSVS = [
    csv
    for n in PT_Q_VALUES
    for csv in expand(config["pseudotime"]["pseudobulk"]["csv"],
                      quantile_n=n, bin=range(1, n + 1), allow_missing=True)
]

# Paper eQTL reference for the report: pseudotime bins (NPC-to-*) excluded, as in 05
PT_PAPER_CELL_TYPES = [ct for ct in config["cell_types"] if "-to-" not in ct]

# Q4 bins exported to 05TENSORQTL as pseudo cell types for SuSiE, SMR and TWAS
PT_Q4_EXPORT_TARGETS = expand(config["pseudotime"]["export"]["sentinel"],
                              trajectory=config["trajectories"], bin=range(1, 5))


def get_pseudotime_tensorqtl(wildcards=None):
    files = []
    for n in PT_Q_VALUES:
        for step in ["tensorqtl_perm", "tensorqtl_nom"]:
            files.extend(
                expand(
                    config["pseudotime"][step]["output"],
                    trajectory=config["trajectories"],
                    quantile_n=n,
                    bin=range(1, n + 1),
                    geno_pc=PT_GENO_PC,
                    exp_pc=config["tensorQTL"]["exp_pcs"]
                )
            )
    return files


def pt_export_paths(trajectory, bin):
    """Source (17PSEUDOTIME) and destination (05TENSORQTL) paths for one Q4 bin.

    The bin is exported as the pseudo cell type {trajectory}-Q4-Bin{bin}, using its
    expression PCs from config["exp_pc_map"] and the tensorQTL path patterns, so the
    files land exactly where downstream pipelines look for any other cell type.
    """
    pt_group = f"{trajectory}-Q4-Bin{bin}"
    pc  = config["exp_pc_map"][pt_group]
    src = dict(trajectory=trajectory, quantile_n=4, bin=bin, geno_pc=PT_GENO_PC, exp_pc=pc)
    dst = dict(cell_type=pt_group, norm_method="quantile", geno_pc=PT_GENO_PC, exp_pc=pc)
    pt, tq = config["pseudotime"], config["tensorQTL"]
    return {
        "nom":  (pt["tensorqtl_nom"]["prefix_out"].format(**src),
                 tq["tensorqtl_nom"]["prefix_out"].format(**dst)),
        "perm": (pt["tensorqtl_perm"]["output"].format(**src),
                 tq["tensorqtl_perm"]["output"].format(**dst)),
        "bed":  (pt["prep_tensorqtl"]["exp_out"].format(**src),
                 tq["prep_tensorQTL_input"]["exp_out"].format(**dst)),
        "cov":  (pt["split_covariates"]["output"].format(**src),
                 tq["split_covariates"]["output"].format(**dst)),
    }


def pt_export_inputs(wildcards):
    paths = pt_export_paths(wildcards.trajectory, wildcards.bin)
    src = {key: paths[key][0] for key in ("perm", "bed", "cov")}
    src["nom"] = paths["nom"][0] + ".cis_qtl_pairs.1.parquet"   # stands in for all chromosomes
    return src


rule pseudotime_subset:
    input:  h5ad   = config["pseudotime"]["subset"]["input_h5ad"]
    output: h5ad   = config["pseudotime"]["subset"]["h5ad"]
    conda:  config["pseudotime"]["env"]
    resources: threads = 2, mem_mb = 180000, time = "00:30:00"
    threads: 2
    log:    config["pseudotime"]["subset"]["log"]
    benchmark: "reports/benchmarks/17pseudotime.pseudotime_subset_{trajectory}.benchmark.txt"
    message: "Subsetting trajectory: {wildcards.trajectory}"
    script: config["pseudotime"]["subset"]["script"]

rule pseudotime_palantir:
    input:  h5ad   = rules.pseudotime_subset.output.h5ad
    output: h5ad   = config["pseudotime"]["palantir"]["h5ad"]
    params: n_neighbors   = config["pseudotime"]["palantir"]["n_neighbors"],
            n_pcs         = config["pseudotime"]["palantir"]["n_pcs"],
            n_dcs         = config["pseudotime"]["palantir"]["n_dcs"],
            num_waypoints = config["pseudotime"]["palantir"]["num_waypoints"],
            root_gene     = config["pseudotime"]["palantir"]["root_gene"]
    conda:  config["pseudotime"]["env"]
    resources: threads = 4, mem_mb = 96000, time = "01:00:00"
    threads: 4
    log:    config["pseudotime"]["palantir"]["log"]
    benchmark: "reports/benchmarks/17pseudotime.pseudotime_palantir_{trajectory}.benchmark.txt"
    message: "Running Palantir pseudotime for trajectory: {wildcards.trajectory}"
    script: config["pseudotime"]["palantir"]["script"]

rule pseudotime_pseudobulk:
    input:  h5ad     = rules.pseudotime_palantir.output.h5ad,
            geno     = config["pseudotime"]["pseudobulk"]["geno_sample_file"]
    output: sentinel = config["pseudotime"]["pseudobulk"]["sentinel"],
            csvs     = PT_PSEUDOBULK_CSVS
    params: q_values = PT_Q_VALUES
    conda:  config["pseudotime"]["env"]
    resources: threads = 2, mem_mb = 90000, time = "00:15:00"
    threads: 2
    log:    config["pseudotime"]["pseudobulk"]["log"]
    benchmark: "reports/benchmarks/17pseudotime.pseudotime_pseudobulk_{trajectory}.benchmark.txt"
    message: "Generating pseudobulk quantile bins for trajectory: {wildcards.trajectory}"
    script: config["pseudotime"]["pseudobulk"]["script"]

rule pseudotime_prep_tensorqtl:
    input:  cov_file    = config["pseudotime"]["prep_tensorqtl"]["cov_file"],
            sex_file    = config["pseudotime"]["prep_tensorqtl"]["sex_file"],
            gene_lookup = config["pseudotime"]["prep_tensorqtl"]["gene_lookup"],
            pseudobulk  = config["pseudotime"]["pseudobulk"]["csv"]
    output: exp_out     = config["pseudotime"]["prep_tensorqtl"]["exp_out"],
            cov_out     = config["pseudotime"]["prep_tensorqtl"]["cov_out"],
            report_rds  = config["pseudotime"]["prep_tensorqtl"]["report_rds"]
    params: cell_type   = "{trajectory}-Q{quantile_n}-Bin{bin}",
            norm_method = "quantile"
    singularity: config["containers"]["r_eqtl"]
    resources: threads = 1, mem_mb = 6000, time = "2:00:00"
    log:    config["pseudotime"]["prep_tensorqtl"]["log"]
    benchmark: "reports/benchmarks/17pseudotime.pseudotime_prep_tensorqtl_{trajectory}_Q{quantile_n}_bin{bin}.benchmark.txt"
    message: "Prepping tensorQTL inputs for {wildcards.trajectory} Q{wildcards.quantile_n} bin {wildcards.bin}"
    script: config["pseudotime"]["prep_tensorqtl"]["script"]

rule pseudotime_split_covariates:
    input:  rules.pseudotime_prep_tensorqtl.output.cov_out
    output: config["pseudotime"]["split_covariates"]["output"]
    params: cell_type = "{trajectory}-Q{quantile_n}-Bin{bin}"
    singularity: config["containers"]["r_eqtl"]
    resources: threads = 1, mem_mb = 6000, time = "1:00:00"
    log:    config["pseudotime"]["split_covariates"]["log"]
    benchmark: "reports/benchmarks/17pseudotime.pseudotime_split_covariates_{trajectory}_Q{quantile_n}_bin{bin}_genPC{geno_pc}_expPC{exp_pc}.benchmark.txt"
    message: "Splitting covariates for {wildcards.trajectory} Q{wildcards.quantile_n} bin {wildcards.bin} genPC {wildcards.geno_pc} expPC {wildcards.exp_pc}"
    script: config["pseudotime"]["split_covariates"]["script"]

rule pseudotime_zip_bed:
    input:  rules.pseudotime_prep_tensorqtl.output.exp_out
    output: config["pseudotime"]["zip_bed"]["output"]
    envmodules: "HTSlib/1.21-GCC-13.3.0"
    resources: threads = 1, mem_mb = 2000, time = "15:00"
    log:    config["pseudotime"]["zip_bed"]["log"]
    benchmark: "reports/benchmarks/17pseudotime.pseudotime_zip_bed_{trajectory}_Q{quantile_n}_bin{bin}.benchmark.txt"
    message: "bgzip and tabix BED for {wildcards.trajectory} Q{wildcards.quantile_n} bin {wildcards.bin}"
    shell:
            """
            bgzip -c {input} > {output}
            tabix -p bed {output}
            """

rule pseudotime_tensorqtl_nom:
    input:  genotypes  = config["tensorQTL"]["convert_genotypes"]["output"],
            counts     = rules.pseudotime_zip_bed.output,
            covariates = rules.pseudotime_split_covariates.output
    output: config["pseudotime"]["tensorqtl_nom"]["output"]
    params: prefix_in  = config["tensorQTL"]["tensorqtl_nom"]["prefix_in"],
            prefix_out = config["pseudotime"]["tensorqtl_nom"]["prefix_out"],
            window     = config["tensorQTL"]["window"]
    singularity: config["containers"]["tensorqtl"]
    resources: threads = 4, mem_mb = 20000, time = "30:00"
    threads: 4
    log:    config["pseudotime"]["tensorqtl_nom"]["log"]
    benchmark: "reports/benchmarks/17pseudotime.pseudotime_tensorqtl_nom_{trajectory}_Q{quantile_n}_bin{bin}_genPC{geno_pc}_expPC{exp_pc}.benchmark.txt"
    message: "tensorQTL nominal for {wildcards.trajectory} Q{wildcards.quantile_n} bin {wildcards.bin} expPC {wildcards.exp_pc}"
    shell:
            """
            python3 -m tensorqtl {params.prefix_in} {input.counts} {params.prefix_out} \
               --covariates {input.covariates} \
               --window {params.window} \
               --mode cis_nominal > {log} 2>&1
            """

rule pseudotime_tensorqtl_perm:
    input:  genotypes  = config["tensorQTL"]["convert_genotypes"]["output"],
            counts     = rules.pseudotime_zip_bed.output,
            covariates = rules.pseudotime_split_covariates.output
    output: config["pseudotime"]["tensorqtl_perm"]["output"]
    params: prefix_in  = config["tensorQTL"]["tensorqtl_perm"]["prefix_in"],
            prefix_out = config["pseudotime"]["tensorqtl_perm"]["prefix_out"],
            window     = config["tensorQTL"]["window"]
    singularity: config["containers"]["tensorqtl"]
    resources: threads = 4, mem_mb = 20000, time = "2:00:00"
    threads: 4
    log:    config["pseudotime"]["tensorqtl_perm"]["log"]
    benchmark: "reports/benchmarks/17pseudotime.pseudotime_tensorqtl_perm_{trajectory}_Q{quantile_n}_bin{bin}_genPC{geno_pc}_expPC{exp_pc}.benchmark.txt"
    message: "tensorQTL permutation for {wildcards.trajectory} Q{wildcards.quantile_n} bin {wildcards.bin} expPC {wildcards.exp_pc}"
    shell:
            """
            python3 -m tensorqtl {params.prefix_in} {input.counts} {params.prefix_out} \
               --covariates {input.covariates} \
               --window {params.window} \
               --mode cis > {log} 2>&1
            """

rule pseudotime_export_smr_twas_inputs:
    input:  unpack(pt_export_inputs)
    output: touch(config["pseudotime"]["export"]["sentinel"])
    log:    config["pseudotime"]["export"]["log"]
    benchmark: "reports/benchmarks/17pseudotime.pseudotime_export_smr_twas_inputs_{trajectory}_Q4_bin{bin}.benchmark.txt"
    wildcard_constraints:
        trajectory = "|".join(config["trajectories"]),
        bin = "[1-4]"
    message: "Exporting {wildcards.trajectory} Q4 bin {wildcards.bin} inputs for SMR/TWAS/SuSiE"
    run:
        import glob, os, shutil

        paths  = pt_export_paths(wildcards.trajectory, wildcards.bin)
        copies = []

        # Nominal parquets: one per chromosome, renamed to the 05 prefix
        src_prefix, dest_prefix = paths["nom"]
        for src in sorted(glob.glob(f"{src_prefix}.cis_qtl_pairs.*.parquet")):
            chrom_part = src.split(".cis_qtl_pairs.")[1]
            copies.append((src, f"{dest_prefix}.cis_qtl_pairs.{chrom_part}"))

        # Permutation output, expression BED and split covariates
        copies += [paths[key] for key in ("perm", "bed", "cov")]

        with open(log[0], "w") as log_fh:
            for src, dest in copies:
                os.makedirs(os.path.dirname(dest), exist_ok=True)
                shutil.copy(src, dest)
                log_fh.write(f"{src} -> {dest}\n")

rule pseudotime_report:
    input:  perm           = get_pseudotime_tensorqtl(),
            rmd_script     = config["pseudotime"]["report"]["rmd_script"]
    output: config["pseudotime"]["report"]["html"]
    params: in_dir              = config["pseudotime"]["report"]["in_dir"],
            cell_type_perm_dir  = config["pseudotime"]["report"]["cell_type_perm_dir"],
            output_file         = config["pseudotime"]["report"]["output_file"],
            bmark_dir           = config["pseudotime"]["report"]["bmark_dir"],
            gene_lookup         = config["pseudotime"]["report"]["gene_lookup"],
            trajectories        = ",".join(config["trajectories"]),
            q_values            = ",".join(map(str, PT_Q_VALUES)),
            exp_pcs             = ",".join(map(str, config["tensorQTL"]["exp_pcs"])),
            geno_pc             = PT_GENO_PC,
            window              = config["tensorQTL"]["window"],
            paper_cell_types    = ",".join(PT_PAPER_CELL_TYPES),
            paper_exp_pcs       = ",".join(str(config["exp_pc_map"][ct]) for ct in PT_PAPER_CELL_TYPES),
            n_neighbors         = config["pseudotime"]["palantir"]["n_neighbors"],
            n_pcs               = config["pseudotime"]["palantir"]["n_pcs"],
            n_dcs               = config["pseudotime"]["palantir"]["n_dcs"],
            num_waypoints       = config["pseudotime"]["palantir"]["num_waypoints"],
            root_gene           = config["pseudotime"]["palantir"]["root_gene"]
    singularity: config["containers"]["r_eqtl"]
    resources: threads = 1, mem_mb = 8000, time = "1:00:00"
    log:    config["pseudotime"]["report"]["log"]
    benchmark: "reports/benchmarks/17pseudotime.pseudotime_report.benchmark.txt"
    message: "Generating pseudotime eQTL report"
    shell:
        """
        Rscript -e "rmarkdown::render('{input.rmd_script}', \
            output_file = '{params.output_file}', \
            params = list(
                in_dir             = '{params.in_dir}',
                cell_type_perm_dir = '{params.cell_type_perm_dir}',
                bmark_dir          = '{params.bmark_dir}',
                gene_lookup        = '{params.gene_lookup}',
                trajectories       = '{params.trajectories}',
                q_values           = '{params.q_values}',
                exp_pcs            = '{params.exp_pcs}',
                geno_pc            = '{params.geno_pc}',
                window             = '{params.window}',
                paper_cell_types   = '{params.paper_cell_types}',
                paper_exp_pcs      = '{params.paper_exp_pcs}',
                n_neighbors        = '{params.n_neighbors}',
                n_pcs              = '{params.n_pcs}',
                n_dcs              = '{params.n_dcs}',
                num_waypoints      = '{params.num_waypoints}',
                root_gene          = '{params.root_gene}'
            ))" > {log} 2>&1
        """
