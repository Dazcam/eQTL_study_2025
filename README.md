## Repo for eQTL study 2024

+ We have 96 samples, on MEGA kit v2

TODO


FOR MARK:

- Relevant pipeline is in workflow/rules/01parse.smk
  - rules cat_fq and run_parse are the main rules to consider
  - resource directive sets resources for individual jobs per rule
  - cat_fq is run locally in order to access fqs on neurocluster
- Individual fastqs are pulled and pooled from a json file: samples_plate1.json
- Slurm paramters are set in: config/profile/config.yaml
- Jobs sbatched via snakemake.sh
- rule all is in Snakefile
- Some logs and benchmarks in workflow/reports dir
