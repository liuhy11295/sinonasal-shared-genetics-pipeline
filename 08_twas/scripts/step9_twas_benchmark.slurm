#!/bin/bash
#SBATCH -J step9_twas_bench
#SBATCH -o %j.step9_twas_bench.out
#SBATCH -e %j.step9_twas_bench.err
#SBATCH -p compute_nodes
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=16G
#SBATCH -t 1-00:00

set -euo pipefail

: "${PROJECT_ROOT:?Set PROJECT_ROOT}"
: "${CONDA_ROOT:?Set CONDA_ROOT}"
: "${PHASE0_CONDA_ENV:?Set PHASE0_CONDA_ENV}"
source "${CONDA_ROOT}/etc/profile.d/conda.sh"
conda activate "${PHASE0_CONDA_ENV}"

cd "$PROJECT_ROOT"
mkdir -p results_end/twas/logs

if [ ! -s results_end/twas/twas_input_qc.tsv ]; then
  python results_end/twas/scripts/step9_twas_pipeline.py input_qc
fi
BENCH_PAIR=$(python results_end/twas/scripts/step9_twas_pipeline.py benchmark_pair)
echo "$BENCH_PAIR" | tee results_end/twas/logs/benchmark_pair.txt
rm -rf "results_end/twas/per_pair/$BENCH_PAIR"
python results_end/twas/scripts/step9_twas_pipeline.py benchmark

cat results_end/twas/resource_benchmark.tsv
cat results_end/twas/parallelization_plan.tsv
