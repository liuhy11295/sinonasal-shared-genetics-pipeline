#!/bin/bash
#SBATCH -J step9_twas_agg
#SBATCH -o %j.step9_twas_agg.out
#SBATCH -e %j.step9_twas_agg.err
#SBATCH -p compute_nodes
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=32G
#SBATCH -t 06:00:00

set -euo pipefail
: "${PROJECT_ROOT:?Set PROJECT_ROOT}"
: "${CONDA_ROOT:?Set CONDA_ROOT}"
: "${PHASE0_CONDA_ENV:?Set PHASE0_CONDA_ENV}"
source "${CONDA_ROOT}/etc/profile.d/conda.sh"
conda activate "${PHASE0_CONDA_ENV}"
cd "$PROJECT_ROOT"
python results_end/twas/scripts/step9_twas_pipeline.py aggregate
