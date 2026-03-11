#!/bin/bash
#SBATCH --job-name=tar_eqtl
#SBATCH --output=tar_eqtl.%j.out
#SBATCH --error=tar_eqtl.%j.err
#SBATCH --time=10:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --partition=c_compute_neuro1
#SBATCH --account=scw1641

set -euo pipefail

# Directories and output file
outdir=/scratch/c.c1477909/data_sharing
nom_dir=/scratch/c.c1477909/eQTL_study_2025/results/10SMR/smr_input
perm_dir=/scratch/c.c1477909/eQTL_study_2025/results/05TENSORQTL/tensorqtl_perm
out_tar=${outdir}/L2_nom_and_perm_eqtl.tar.gz

#cell_types="ExN-UL ExN-DL InN RG MG Endo-Peri OPC"

cell_types="ExN-UL-0 ExN-UL-1 ExN-UL-2 ExN-DL-0 ExN-DL-1 ExN-DL-2 ExN-DL-3 InN-0 InN-1 RG-0 RG-1 RG-2 RG-3 "

echo "[INFO] Creating output directory: ${outdir}"
mkdir -p "${outdir}"

module load pigz/2.4

# Build tar arguments array
echo "[INFO] Preparing file list for archiving..."
args=()
for ct in ${cell_types}; do
  args+=("-C" "${nom_dir}/${ct}" "${ct}_nom.cis_qtl_pairs.tsv")
  args+=("-C" "${perm_dir}/${ct}_quantile_genPC_4_expPC_40" "${ct}_quantile_perm.cis_qtl.txt.gz")
done

echo "[INFO] Creating tarball and compressing with pigz..."
tar -cf - "${args[@]}" | pigz -6 -p "$SLURM_CPUS_PER_TASK" > "${out_tar}"

echo "[INFO] Tarball created at ${out_tar}"

# Verification step: test gzip integrity with pigz
echo "[INFO] Verifying integrity of the compressed tarball..."
if pigz -t "${out_tar}"; then
  echo "[SUCCESS] Integrity check passed. The archive is valid."
else
  echo "[ERROR] Integrity check failed! The archive may be corrupt."
  exit 1
fi

echo "[INFO] Job completed successfully."
