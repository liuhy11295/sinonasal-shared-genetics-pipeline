#!/bin/bash
#SBATCH -J step9_twas_full
#SBATCH -o %A_%a.step9_twas_full.out
#SBATCH -e %A_%a.step9_twas_full.err
#SBATCH -D /platform_data/p_user/p010/phase0
#SBATCH -p compute_nodes
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=8G
#SBATCH -t 1-00:00

set -euo pipefail
export PROJECT_ROOT=/platform_data/p_user/p010/phase0
export CONDA_ROOT=/platform_data/p_user/p010/miniconda3
export PHASE0_CONDA_ENV=/platform_data/p_user/p010/conda_envs/nasal_pairwise_env
source "${CONDA_ROOT}/etc/profile.d/conda.sh"
conda activate "${PHASE0_CONDA_ENV}"
cd "$PROJECT_ROOT"
python results_end/twas/scripts/step9_twas_pipeline.py array_task
