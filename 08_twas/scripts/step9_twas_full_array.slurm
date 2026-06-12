#!/bin/bash
#SBATCH -J step9_twas_full
#SBATCH -o %A_%a.step9_twas_full.out
#SBATCH -e %A_%a.step9_twas_full.err
#SBATCH -p compute_nodes
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=8G
#SBATCH -t 1-00:00

set -euo pipefail
: "${PROJECT_ROOT:?Set PROJECT_ROOT}"
: "${CONDA_ROOT:?Set CONDA_ROOT}"
: "${PHASE0_CONDA_ENV:?Set PHASE0_CONDA_ENV}"
source "${CONDA_ROOT}/etc/profile.d/conda.sh"
conda activate "${PHASE0_CONDA_ENV}"
cd "$PROJECT_ROOT"
python results_end/twas/scripts/step9_twas_pipeline.py array_task
