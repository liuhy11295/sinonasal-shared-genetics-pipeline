#!/bin/bash
#SBATCH -J step9_rdeps
#SBATCH -o %j.step9_rdeps.out
#SBATCH -e %j.step9_rdeps.err
#SBATCH -D /platform_data/p_user/p010/phase0
#SBATCH -p compute_nodes
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem=16G
#SBATCH -t 0-02:00

set -euo pipefail

CONDA_ROOT=/platform_data/p_user/p010/miniconda3
ENV=/platform_data/p_user/p010/conda_envs/nasal_pairwise_env
source "${CONDA_ROOT}/etc/profile.d/conda.sh"
conda activate "${ENV}"

conda install -y -p "${ENV}" -c conda-forge r-here r-optparse

Rscript -e 'cat("R=", R.version.string, "\n"); library(here); library(optparse); cat("FUSION R deps OK\n")'
