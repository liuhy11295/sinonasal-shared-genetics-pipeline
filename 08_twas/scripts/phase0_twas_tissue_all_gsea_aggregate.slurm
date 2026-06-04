#!/bin/bash
#SBATCH -J twas_tissue31_ag
#SBATCH -o results_end/twas/pathway_overlap_tissue_all31/logs/%A.aggregate.out
#SBATCH -e results_end/twas/pathway_overlap_tissue_all31/logs/%A.aggregate.err
#SBATCH -D /platform_data/p_user/p010/phase0
#SBATCH -p compute_nodes
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH -t 01:00:00

set -euo pipefail

export PROJECT_ROOT=/platform_data/p_user/p010/phase0
export CONDA_ROOT=/platform_data/p_user/p010/miniconda3
export PHASE0_CONDA_ENV=/platform_data/p_user/p010/conda_envs/nasal_pairwise_env

source "${CONDA_ROOT}/etc/profile.d/conda.sh"
conda activate "${PHASE0_CONDA_ENV}"

cd "${PROJECT_ROOT}"
python results_end/twas/scripts/phase0_twas_tissue_all_gsea_overlap.py aggregate
