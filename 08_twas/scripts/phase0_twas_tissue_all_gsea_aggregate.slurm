#!/bin/bash
#SBATCH -J twas_tissue31_ag
#SBATCH -o results_end/twas/pathway_overlap_tissue_all31/logs/%A.aggregate.out
#SBATCH -e results_end/twas/pathway_overlap_tissue_all31/logs/%A.aggregate.err
#SBATCH -p compute_nodes
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH -t 01:00:00

set -euo pipefail

: "${PROJECT_ROOT:?Set PROJECT_ROOT}"
: "${CONDA_ROOT:?Set CONDA_ROOT}"
: "${PHASE0_CONDA_ENV:?Set PHASE0_CONDA_ENV}"

source "${CONDA_ROOT}/etc/profile.d/conda.sh"
conda activate "${PHASE0_CONDA_ENV}"

cd "${PROJECT_ROOT}"
python results_end/twas/scripts/phase0_twas_tissue_all_gsea_overlap.py aggregate
