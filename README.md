# sc-eQTL Atlas | Developing Human Brain | Beta

> **A single-nucleus eQTL atlas of the prenatal human cerebral cortex**  
> Cardiff University · Division of Psychological Medicine and Clinical Neurosciences  
> Manuscript in preparation · 2026

[![Docs](https://github.com/Dazcam/eQTL_study_2025/actions/workflows/publish.yml/badge.svg)](https://dazcam.github.io/eQTL_study_2025/)
[![License: CC BY 4.0](https://img.shields.io/badge/License-CC%20BY%204.0-lightgrey.svg)](https://creativecommons.org/licenses/by/4.0/)
[![Snakemake](https://img.shields.io/badge/snakemake-≥8.0-brightgreen)](https://snakemake.readthedocs.io)
[![Platform](https://img.shields.io/badge/platform-SLURM%20HPC-blue)](https://slurm.schedmd.com/)

## What this is

We performed **single-nucleus RNA sequencing and genome-wide genotyping** on cerebral cortex from 134 unrelated samples (second trimester) to generate the first cell-type-resolved eQTL atlas of the prenatal human brain.

This repository is an end-to-end computational genomics platform to process ~3 TB of raw single-nucleus RNA sequencing and genome-wide genotyping data. The pipeline identifies genetic variants that influence gene expression in specific brain cell types during development, and links those variants to neuropsychiatric disease risk.

## Pipeline overview

The full analysis runs as 13 sequential Snakemake pipelines on a SLURM HPC cluster, with up to 500 concurrent jobs. Each stage is containerised via Conda environments or Singularity containers for reproducibility.

| # | Stage | Tools |
| :---: | :--- | :--- |
| 01 | scRNA-seq Alignment | Parse Biosciences split-pipe |
| 02 | scRNA-seq QC, clustering & pseudobulk | Scanpy, Papermill |
| 03 | Genotype QC — pre-imputation | PLINK, GenotypeQCtoHRC |
| 04 | Genotype QC — post-imputation | bcftools, PLINK, dbSNP |
| 05 | cis-eQTL mapping | TensorQTL (GPU-accelerated) |
| 06 | eQTL replication & ATAC enrichment | π₁ statistic, custom R |
| 07 | GWAS summary statistic preparation | GWAS standarisation and liftover |
| 08 | Bayesian fine-mapping | SuSiE |
| 09 | Heritability enrichment | S-LDSR (stratified LD score regression) |
| 10 | Causal inference | SMR + HEIDI |
| 11 | TWAS weight computation | FUSION |
| 12 | Causal TWAS | cTWAS |
| 13 | Manuscript figures, tables & data sharing | R (ggplot2, tidyverse) |

Full documentation for each stage — including DAG diagrams, rule descriptions, resource profiles, and technical notes — is at the **[pipeline documentation site](https://dazcam.github.io/eQTL_study_2025/)**.

## App

A fully client-side interactive app to query the eQTL permutation results is available at the **[eQTL Browser](https://dazcam.github.io/eQTL_study_2025/eqtl-browser/)**.

## Repository structure

```
.
├── config/                   # Pipeline configuration (YAML, sample lists, JSON)
├── workflow/
│   ├── Snakefile             # Master workflow entry point
│   ├── rules/                # 13 Snakemake rule files (one per pipeline stage)
│   ├── scripts/              # R, Python, and shell analysis scripts
│   └── envs/                 # Conda environment definitions
├── pipelines/                # Quarto documentation source (one .qmd per stage)
├── resources/                # Reference data, public datasets, sample metadata
├── results/                  # Pipeline outputs (not tracked in git)
└── eqtl-browser/             # Interactive eQTL data browser (static HTML app)
```

## Quickstart

### Prerequisites

| Requirement | Version |
| :--- | :--- |
| Snakemake | ≥ 8.0 |
| Conda / Mamba | any |
| Singularity | ≥ 3.x |
| SLURM | any |
| Python | 3.12 |

See the **[Setup page](https://dazcam.github.io/eQTL_study_2025/setup.html)** for full installation instructions, data access, and environment configuration.

### Run

```bash
# Clone the repo
git clone https://github.com/Dazcam/eQTL_study_2025.git
cd eQTL_study_2025

# Dry run — check the full job graph without executing
snakemake -n --profile config/profile

# Submit to SLURM
bash workflow/snakemake.sh
```

## Data access

| Resource | Location |
| :--- | :--- |
| Raw snRNA-seq data | *Available on request pending publication* |
| Imputed genotypes | *Available on request pending publication* |
| eQTL summary statistics | *To be deposited on publication* |
| TWAS weights | *To be deposited on publication* |

## Citation

> Manuscript in preparation. Details will be updated on preprint release.

## Licence & copyright

Code in this repository is released under [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/).  
© 2026 Cardiff University. See [LICENCE.md](LICENCE.md) for details.
