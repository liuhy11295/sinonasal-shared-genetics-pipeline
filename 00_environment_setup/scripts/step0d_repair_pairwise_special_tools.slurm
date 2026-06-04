#!/bin/bash
#SBATCH -J nasal_step0d_repair
#SBATCH -o /platform_data/p_user/p010/phase0/results/phase0_extension/step0_1/logs/%j.step0d_repair.out
#SBATCH -e /platform_data/p_user/p010/phase0/results/phase0_extension/step0_1/logs/%j.step0d_repair.err
#SBATCH -D /platform_data/p_user/p010/phase0
#SBATCH -p compute_nodes
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem=64G
#SBATCH -t 0-08:00:00

set -euo pipefail
BASE=/platform_data/p_user/p010
PROJECT=$BASE/phase0
WORK=$PROJECT/results/phase0_extension/step0_1
CONDA_ROOT=$BASE/miniconda3
ENV=$BASE/conda_envs/nasal_pairwise_env
TOOLS=$WORK/tools
CHECK=$WORK/checks
LOG=$WORK/logs
mkdir -p "$TOOLS" "$CHECK" "$LOG"
exec > >(tee -a "$LOG/step0d_repair_pairwise_special_tools.runtime.log") 2>&1

echo "STEP0D_START=$(date -Is)"
echo "HOST=$(hostname)"
echo "ENV=$ENV"
source "$CONDA_ROOT/etc/profile.d/conda.sh"

# Keep everything in nasal_pairwise_env unless a true version conflict appears.
conda install -y -p "$ENV" --override-channels -c conda-forge -c bioconda \
  r-gmp r-arrangements r-iterpc r-cairo r-ieugwasr r-rsnps r-mendelianrandomization \
  r-remotes r-devtools git curl unzip ca-certificates certifi openssl || true
conda activate "$ENV"

STATUS=$CHECK/special_tools_versions.tsv
MISSING=$CHECK/missing_special_tools.tsv
: > "$STATUS"
: > "$MISSING"
printf 'tool\tstatus\tversion_or_path\tnote\n' >> "$STATUS"
printf 'tool\treason\tnote\n' >> "$MISSING"
record_ok() { printf '%s\tOK\t%s\t%s\n' "$1" "$2" "$3" >> "$STATUS"; }
record_missing() { printf '%s\tMISSING\t%s\t%s\n' "$1" "$2" "$3" >> "$STATUS"; printf '%s\t%s\t%s\n' "$1" "$2" "$3" >> "$MISSING"; }

Rscript - <<'RS'
options(repos = c(CRAN = 'https://cloud.r-project.org'))
for (p in c('gmp','arrangements','iterpc','ieugwasr')) {
  if (!requireNamespace(p, quietly=TRUE)) {
    tryCatch(install.packages(p, dependencies=TRUE), error=function(e) message(p, ' CRAN install failed: ', conditionMessage(e)))
  }
}
options(repos = c(MRCIEU = 'https://mrcieu.r-universe.dev', CRAN = 'https://cloud.r-project.org'))
if (!requireNamespace('remotes', quietly=TRUE)) install.packages('remotes')
if (!requireNamespace('TwoSampleMR', quietly=TRUE)) {
  tryCatch(install.packages('TwoSampleMR', dependencies=TRUE), error=function(e) message('TwoSampleMR r-universe failed: ', conditionMessage(e)))
}
if (!requireNamespace('TwoSampleMR', quietly=TRUE)) {
  tryCatch(remotes::install_github('MRCIEU/TwoSampleMR', upgrade='never', dependencies=TRUE), error=function(e) message('TwoSampleMR github failed: ', conditionMessage(e)))
}
if (!requireNamespace('MRPRESSO', quietly=TRUE)) {
  tryCatch(remotes::install_github('rondolab/MR-PRESSO', upgrade='never', dependencies=TRUE), error=function(e) message('MRPRESSO github failed: ', conditionMessage(e)))
}
RS

if Rscript -e "cat(as.character(packageVersion('TwoSampleMR')))" > "$CHECK/TwoSampleMR.version" 2>/dev/null; then record_ok TwoSampleMR "$(cat "$CHECK/TwoSampleMR.version")" "$ENV"; else record_missing TwoSampleMR package_not_available "$ENV"; fi
if Rscript -e "cat(as.character(packageVersion('MRPRESSO')))" > "$CHECK/MRPRESSO.version" 2>/dev/null; then record_ok MRPRESSO "$(cat "$CHECK/MRPRESSO.version")" "$ENV"; else record_missing MRPRESSO package_not_available "$ENV"; fi

# Try lightweight source/script retrieval again, with SSL fallback. These are scripts/tools under p010, not system installs.
PLACODIR=$TOOLS/PLACO
if [ ! -d "$PLACODIR/.git" ]; then
  rm -rf "$PLACODIR"
  git -c http.sslVerify=false clone --depth 1 https://github.com/RayDebashree/PLACO.git "$PLACODIR" || true
fi
if [ -f "$PLACODIR/PLACO_v0.2.0.R" ]; then
  record_ok PLACO "$PLACODIR/PLACO_v0.2.0.R" source_script
else
  # If git is blocked, try raw script URL directly.
  mkdir -p "$PLACODIR"
  curl -k -L --connect-timeout 60 --max-time 180 -o "$PLACODIR/PLACO_v0.2.0.R" https://raw.githubusercontent.com/RayDebashree/PLACO/master/PLACO_v0.2.0.R || true
  if [ -s "$PLACODIR/PLACO_v0.2.0.R" ]; then record_ok PLACO "$PLACODIR/PLACO_v0.2.0.R" source_script_raw; else record_missing PLACO download_failed "GitHub/raw unreachable from compute node"; fi
fi

LCVDIR=$TOOLS/LCV
if [ ! -d "$LCVDIR/.git" ]; then
  rm -rf "$LCVDIR"
  git -c http.sslVerify=false clone --depth 1 https://github.com/lukejoconnor/LCV.git "$LCVDIR" || true
fi
LCV_FOUND=$(find "$LCVDIR" -type f \( -iname '*lcv*' -o -iname '*.R' \) 2>/dev/null | head -n 1 || true)
if [ -n "$LCV_FOUND" ]; then
  record_ok LCV "$LCVDIR" source_repo
else
  mkdir -p "$LCVDIR/R"
  for f in RunLCV.R MomentFunctions.R ExampleRealdataScript.R; do
    curl -k -L --connect-timeout 60 --max-time 180 -o "$LCVDIR/R/$f" "https://raw.githubusercontent.com/lukejoconnor/LCV/master/R/$f" || true
  done
  LCV_FOUND=$(find "$LCVDIR" -type f -name 'RunLCV.R' -size +0c 2>/dev/null | head -n 1 || true)
  if [ -n "$LCV_FOUND" ]; then record_ok LCV "$LCVDIR" source_scripts_raw; else record_missing LCV download_failed "GitHub/raw unreachable from compute node"; fi
fi

CPDIR=$TOOLS/CPASSOC
CP_FOUND=$(find "$CPDIR" -type f \( -iname '*cpassoc*' -o -iname '*.R' -o -perm -111 \) 2>/dev/null | head -n 1 || true)
if [ -n "$CP_FOUND" ]; then record_ok CPASSOC "$CP_FOUND" downloaded; else record_missing CPASSOC download_failed "$CPDIR"; fi

conda list -p "$ENV" > "$CHECK/nasal_pairwise_env.conda_list.after_special_tools.txt" || true
Rscript - <<'RS' > "$CHECK/r_package_versions.after_special_tools.tsv"
pkgs <- c('TwoSampleMR','ieugwasr','MRPRESSO','gmp','arrangements','iterpc','data.table','coloc','susieR')
cat('package\tstatus\tversion\n')
for (p in pkgs) {
  ok <- suppressWarnings(requireNamespace(p, quietly=TRUE))
  if (ok) cat(p, 'OK', as.character(utils::packageVersion(p)), sep='\t') else cat(p, 'MISSING', 'NA', sep='\t')
  cat('\n')
}
RS

echo "STEP0D_DONE=$(date -Is)"
