#!/bin/bash
#SBATCH -J step9_twas_agg
#SBATCH -o %j.step9_twas_agg.out
#SBATCH -e %j.step9_twas_agg.err
#SBATCH -D /platform_data/p_user/p010/phase0
#SBATCH -p compute_nodes
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=32G
#SBATCH -t 06:00:00

set -euo pipefail
export PROJECT_ROOT=/platform_data/p_user/p010/phase0
export CONDA_ROOT=/platform_data/p_user/p010/miniconda3
export PHASE0_CONDA_ENV=/platform_data/p_user/p010/conda_envs/nasal_pairwise_env
source "${CONDA_ROOT}/etc/profile.d/conda.sh"
conda activate "${PHASE0_CONDA_ENV}"
cd "$PROJECT_ROOT"
python results_end/twas/scripts/step9_twas_pipeline.py aggregate
