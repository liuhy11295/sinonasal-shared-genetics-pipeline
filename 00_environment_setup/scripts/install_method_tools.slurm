#!/bin/bash
#SBATCH -J nasal_method_tools
#SBATCH -o %j.install_method_tools.out
#SBATCH -e %j.install_method_tools.err
#SBATCH -p compute_nodes
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem=64G
#SBATCH -t 0-12:00:00

set -euo pipefail
PROJECT="${PROJECT_ROOT:?Set PROJECT_ROOT}"
WORK=$PROJECT/results/phase0_extension/step0_1
ENV_ROOT="${ENV_ROOT:?Set ENV_ROOT}"
CONDA_ROOT="${CONDA_ROOT:?Set CONDA_ROOT}"
TOOLS=$WORK/tools
CHECK=$WORK/checks
LOG=$WORK/logs
mkdir -p "$ENV_ROOT" "$TOOLS" "$CHECK" "$LOG"
exec > >(tee -a "$LOG/install_method_tools.runtime.log") 2>&1

echo "STEP0B_START=$(date -Is)"
echo "HOST=$(hostname)"
echo "CONDA_ROOT=$CONDA_ROOT"
source "$CONDA_ROOT/etc/profile.d/conda.sh"

STATUS=$CHECK/special_tools_versions.tsv
MISSING=$CHECK/missing_special_tools.tsv
: > "$STATUS"
: > "$MISSING"
printf 'tool\tstatus\tversion_or_path\tnote\n' >> "$STATUS"
printf 'tool\treason\tnote\n' >> "$MISSING"

record_ok() { printf '%s\tOK\t%s\t%s\n' "$1" "$2" "$3" >> "$STATUS"; }
record_missing() { printf '%s\tMISSING\t%s\t%s\n' "$1" "$2" "$3" >> "$STATUS"; printf '%s\t%s\t%s\n' "$1" "$2" "$3" >> "$MISSING"; }

create_or_update_env() {
  local env_path="$1"; shift
  if [ ! -d "$env_path" ]; then
    conda create -y -p "$env_path" --override-channels -c conda-forge -c bioconda "$@"
  else
    conda install -y -p "$env_path" --override-channels -c conda-forge -c bioconda "$@"
  fi
}

# MR env: install into p010-owned miniconda env only.
MR_ENV=$ENV_ROOT/nasal_mr_env
echo "## Preparing MR env: $MR_ENV"
create_or_update_env "$MR_ENV" \
  r-base=4.3 r-remotes r-devtools r-data.table r-dplyr r-tidyr r-readr r-stringr r-ggplot2 \
  r-optparse r-argparse r-jsonlite r-httr r-curl r-openssl r-xml2 r-vroom r-matrix \
  r-rcpp r-rcppeigen r-glmnet r-ape r-psych r-metafor git curl unzip make gcc_linux-64 gxx_linux-64
conda activate "$MR_ENV"
Rscript - <<'RS'
options(repos = c(MRCIEU = 'https://mrcieu.r-universe.dev', CRAN = 'https://cloud.r-project.org'))
install_if_missing <- function(pkg) {
  if (!requireNamespace(pkg, quietly=TRUE)) install.packages(pkg, dependencies=TRUE)
}
install_if_missing('remotes')
if (!requireNamespace('TwoSampleMR', quietly=TRUE)) {
  tryCatch(install.packages('TwoSampleMR', dependencies=TRUE), error=function(e) message('TwoSampleMR install.packages failed: ', conditionMessage(e)))
}
if (!requireNamespace('TwoSampleMR', quietly=TRUE)) {
  tryCatch(remotes::install_github('MRCIEU/TwoSampleMR', upgrade='never', dependencies=TRUE), error=function(e) message('TwoSampleMR github failed: ', conditionMessage(e)))
}
if (!requireNamespace('MRPRESSO', quietly=TRUE)) {
  tryCatch(remotes::install_github('rondolab/MR-PRESSO', upgrade='never', dependencies=TRUE), error=function(e) message('MRPRESSO github failed: ', conditionMessage(e)))
}
RS
if Rscript -e "cat(as.character(packageVersion('TwoSampleMR')))" > "$CHECK/TwoSampleMR.version" 2>/dev/null; then record_ok TwoSampleMR "$(cat "$CHECK/TwoSampleMR.version")" "$MR_ENV"; else record_missing TwoSampleMR package_not_available "$MR_ENV"; fi
if Rscript -e "cat(as.character(packageVersion('MRPRESSO')))" > "$CHECK/MRPRESSO.version" 2>/dev/null; then record_ok MRPRESSO "$(cat "$CHECK/MRPRESSO.version")" "$MR_ENV"; else record_missing MRPRESSO package_not_available "$MR_ENV"; fi
conda deactivate || true

# PLACO / CPASSOC env and source/script downloads.
PC_ENV=$ENV_ROOT/nasal_placo_cpassoc_env
echo "## Preparing PLACO/CPASSOC env: $PC_ENV"
create_or_update_env "$PC_ENV" r-base=4.3 r-data.table r-dplyr r-readr r-stringr r-optparse r-argparse r-matrix git curl unzip
conda activate "$PC_ENV"
PLACODIR=$TOOLS/PLACO
if [ ! -d "$PLACODIR/.git" ]; then rm -rf "$PLACODIR"; git clone --depth 1 https://github.com/RayDebashree/PLACO.git "$PLACODIR" || true; fi
if [ -f "$PLACODIR/PLACO_v0.2.0.R" ]; then
  record_ok PLACO "$PLACODIR/PLACO_v0.2.0.R" source_script
else
  record_missing PLACO download_failed "https://github.com/RayDebashree/PLACO"
fi
CPDIR=$TOOLS/CPASSOC
mkdir -p "$CPDIR"
if [ ! -f "$CPDIR/CPASSOC.zip" ]; then
  (curl -L --fail --connect-timeout 30 -o "$CPDIR/CPASSOC.zip" "https://hal.case.edu/~xxz10/zhu-web/CPASSOC/CPASSOC.zip" || \
   curl -L --fail --connect-timeout 30 -o "$CPDIR/CPASSOC.zip" "http://hal.case.edu/~xxz10/zhu-web/CPASSOC/CPASSOC.zip" || true)
fi
if [ -s "$CPDIR/CPASSOC.zip" ]; then
  unzip -o "$CPDIR/CPASSOC.zip" -d "$CPDIR/unpacked" >/dev/null 2>&1 || true
fi
CP_FOUND=$(find "$CPDIR" -type f \( -iname '*cpassoc*' -o -iname '*.R' -o -perm -111 \) 2>/dev/null | head -n 1 || true)
if [ -n "$CP_FOUND" ]; then record_ok CPASSOC "$CP_FOUND" downloaded; else record_missing CPASSOC download_failed "CPASSOC.zip not found/unusable"; fi
conda deactivate || true

# LCV env and source download.
LCV_ENV=$ENV_ROOT/nasal_lcv_env
echo "## Preparing LCV env: $LCV_ENV"
create_or_update_env "$LCV_ENV" r-base=4.3 r-data.table r-matrix r-optparse r-argparse git curl unzip
conda activate "$LCV_ENV"
LCVDIR=$TOOLS/LCV
if [ ! -d "$LCVDIR/.git" ]; then rm -rf "$LCVDIR"; git clone --depth 1 https://github.com/lukejoconnor/LCV.git "$LCVDIR" || true; fi
LCV_FOUND=$(find "$LCVDIR" -type f \( -iname '*lcv*' -o -iname '*.R' \) 2>/dev/null | head -n 1 || true)
if [ -n "$LCV_FOUND" ]; then record_ok LCV "$LCVDIR" source_repo; else record_missing LCV download_failed "https://github.com/lukejoconnor/LCV"; fi
conda deactivate || true

# Environment snapshots.
for env in "$MR_ENV" "$PC_ENV" "$LCV_ENV"; do
  b=$(basename "$env")
  conda list -p "$env" > "$CHECK/${b}.conda_list.txt" || true
done

echo "STEP0B_DONE=$(date -Is)"
