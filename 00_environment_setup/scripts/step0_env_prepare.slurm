#!/bin/bash
#SBATCH -J nasal_step0_env
#SBATCH -o /platform_data/p_user/p010/phase0/results/phase0_extension/step0_1/logs/%j.step0_env.out
#SBATCH -e /platform_data/p_user/p010/phase0/results/phase0_extension/step0_1/logs/%j.step0_env.err
#SBATCH -D /platform_data/p_user/p010/phase0
#SBATCH -p compute_nodes
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem=32G
#SBATCH -t 0-08:00:00

set -euo pipefail
BASE=/platform_data/p_user/p010
WORK=$BASE/phase0/results/phase0_extension/step0_1
ENV_ROOT=$BASE/conda_envs
MAIN_ENV=$ENV_ROOT/nasal_pairwise_env
CONDA_ROOT=$BASE/miniconda3
CONDA=$CONDA_ROOT/bin/conda
CHECK_DIR=$WORK/checks
LOG_DIR=$WORK/logs
mkdir -p "$CHECK_DIR" "$LOG_DIR" "$ENV_ROOT"

exec > >(tee -a "$LOG_DIR/step0_env_prepare.runtime.log") 2>&1

echo "STEP0_START=$(date -Is)"
echo "HOST=$(hostname)"
echo "BASE=$BASE"
echo "MAIN_ENV=$MAIN_ENV"

if [ ! -x "$CONDA" ]; then
  echo "ERROR: p010 conda not found at $CONDA"
  exit 10
fi

source "$CONDA_ROOT/etc/profile.d/conda.sh"
conda config --show channels > "$CHECK_DIR/conda_channels.before.txt" || true

CORE_PKGS=(
  python=3.10
  r-base=4.3
  plink
  plink2
  pandas
  numpy
  scipy
  statsmodels
  pyarrow
  polars
  tqdm
  pyyaml
  joblib
  r-data.table
  r-tidyverse
  r-dplyr
  r-readr
  r-stringr
  r-tidyr
  r-ggplot2
  r-optparse
  r-argparse
  r-coloc
  r-susier
  r-matrix
  r-foreach
  r-doparallel
)

if [ ! -d "$MAIN_ENV" ]; then
  echo "Creating $MAIN_ENV"
  conda create -y -p "$MAIN_ENV" --override-channels -c conda-forge -c bioconda "${CORE_PKGS[@]}"
else
  echo "Updating existing $MAIN_ENV"
  conda install -y -p "$MAIN_ENV" --override-channels -c conda-forge -c bioconda "${CORE_PKGS[@]}"
fi

conda activate "$MAIN_ENV"

echo "## conda list" > "$CHECK_DIR/conda_list.txt"
conda list >> "$CHECK_DIR/conda_list.txt"

python - <<'PY' > "$CHECK_DIR/python_package_versions.tsv"
import importlib, sys
mods = ['pandas','numpy','scipy','statsmodels','pyarrow','polars','tqdm','yaml','joblib']
print('package\tstatus\tversion')
for m in mods:
    try:
        mod = importlib.import_module(m)
        print(f'{m}\tOK\t{getattr(mod, "__version__", "NA")}')
    except Exception as e:
        print(f'{m}\tMISSING\t{type(e).__name__}: {e}')
print(f'python\tOK\t{sys.version.split()[0]}')
PY

Rscript - <<'RS' > "$CHECK_DIR/r_package_versions.tsv"
pkgs <- c('data.table','tidyverse','dplyr','readr','stringr','tidyr','ggplot2','optparse','argparse','coloc','susieR','Matrix','foreach','doParallel','TwoSampleMR','MRPRESSO','PLACO')
cat('package\tstatus\tversion\n')
for (p in pkgs) {
  ok <- suppressWarnings(requireNamespace(p, quietly=TRUE))
  if (ok) {
    cat(p, 'OK', as.character(utils::packageVersion(p)), sep='\t'); cat('\n')
  } else {
    cat(p, 'MISSING', 'NA', sep='\t'); cat('\n')
  }
}
RS

Rscript -e 'sessionInfo()' > "$CHECK_DIR/R_sessionInfo.txt"

{
  echo "## plink"
  if command -v plink >/dev/null 2>&1; then plink --version || true; else echo "MISSING"; fi
  echo "## plink2"
  if command -v plink2 >/dev/null 2>&1; then plink2 --version || true; else echo "MISSING"; fi
} > "$CHECK_DIR/plink_versions.txt" 2>&1

cat > "$CHECK_DIR/special_tool_status.tsv" <<'EOF'
tool	status	note
LCV	MISSING	No forced install in Step 0; needs existing script/package discovery before use.
CPASSOC	MISSING	No forced install in Step 0; needs existing script/package discovery before use.
EOF

# Discover existing user-provided scripts/packages for hard-to-install methods without modifying them.
{
  echo "## Existing PLACO/LCV/CPASSOC/MR scripts under p010"
  find "$BASE" -path "$BASE/.cache" -prune -o -type f \( -iname '*placo*' -o -iname '*lcv*' -o -iname '*cpassoc*' -o -iname '*TwoSampleMR*' -o -iname '*MRPRESSO*' -o -iname '*MR-PRESSO*' \) -print 2>/dev/null | sort | head -n 300
} > "$CHECK_DIR/existing_special_tools_scan.txt"

echo "STEP0_DONE=$(date -Is)"
