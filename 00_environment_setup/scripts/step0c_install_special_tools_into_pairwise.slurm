#!/bin/bash
#SBATCH -J nasal_step0c_pairtools
#SBATCH -o %j.step0c_pairtools.out
#SBATCH -e %j.step0c_pairtools.err
#SBATCH -p compute_nodes
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem=64G
#SBATCH -t 0-12:00:00

set -euo pipefail
PROJECT="${PROJECT_ROOT:?Set PROJECT_ROOT}"
WORK=$PROJECT/results/phase0_extension/step0_1
CONDA_ROOT="${CONDA_ROOT:?Set CONDA_ROOT}"
ENV="${PHASE0_CONDA_ENV:?Set PHASE0_CONDA_ENV}"
TOOLS=$WORK/tools
CHECK=$WORK/checks
LOG=$WORK/logs
mkdir -p "$TOOLS" "$CHECK" "$LOG"
exec > >(tee -a "$LOG/step0c_install_special_tools_into_pairwise.runtime.log") 2>&1

echo "STEP0C_START=$(date -Is)"
echo "HOST=$(hostname)"
echo "ENV=$ENV"

source "$CONDA_ROOT/etc/profile.d/conda.sh"
conda activate "$ENV"

STATUS=$CHECK/special_tools_versions.tsv
MISSING=$CHECK/missing_special_tools.tsv
: > "$STATUS"
: > "$MISSING"
printf 'tool\tstatus\tversion_or_path\tnote\n' >> "$STATUS"
printf 'tool\treason\tnote\n' >> "$MISSING"
record_ok() { printf '%s\tOK\t%s\t%s\n' "$1" "$2" "$3" >> "$STATUS"; }
record_missing() { printf '%s\tMISSING\t%s\t%s\n' "$1" "$2" "$3" >> "$STATUS"; printf '%s\t%s\t%s\n' "$1" "$2" "$3" >> "$MISSING"; }

# Add install helpers inside the existing pairwise env only.
conda install -y -p "$ENV" --override-channels -c conda-forge -c bioconda \
  r-remotes r-devtools r-jsonlite r-httr r-curl r-openssl r-xml2 r-vroom \
  r-rcpp r-rcppeigen r-glmnet r-ape r-psych r-metafor git curl unzip make gcc_linux-64 gxx_linux-64 || true
conda activate "$ENV"

Rscript - <<'RS'
options(repos = c(MRCIEU = 'https://mrcieu.r-universe.dev', CRAN = 'https://cloud.r-project.org'))
if (!requireNamespace('remotes', quietly=TRUE)) install.packages('remotes', dependencies=TRUE)
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

if Rscript -e "cat(as.character(packageVersion('TwoSampleMR')))" > "$CHECK/TwoSampleMR.version" 2>/dev/null; then record_ok TwoSampleMR "$(cat "$CHECK/TwoSampleMR.version")" "$ENV"; else record_missing TwoSampleMR package_not_available "$ENV"; fi
if Rscript -e "cat(as.character(packageVersion('MRPRESSO')))" > "$CHECK/MRPRESSO.version" 2>/dev/null; then record_ok MRPRESSO "$(cat "$CHECK/MRPRESSO.version")" "$ENV"; else record_missing MRPRESSO package_not_available "$ENV"; fi

PLACODIR=$TOOLS/PLACO
if [ ! -d "$PLACODIR/.git" ]; then rm -rf "$PLACODIR"; git clone --depth 1 https://github.com/RayDebashree/PLACO.git "$PLACODIR" || true; fi
if [ -f "$PLACODIR/PLACO_v0.2.0.R" ]; then record_ok PLACO "$PLACODIR/PLACO_v0.2.0.R" source_script; else record_missing PLACO download_failed "https://github.com/RayDebashree/PLACO"; fi

LCVDIR=$TOOLS/LCV
if [ ! -d "$LCVDIR/.git" ]; then rm -rf "$LCVDIR"; git clone --depth 1 https://github.com/lukejoconnor/LCV.git "$LCVDIR" || true; fi
LCV_FOUND=$(find "$LCVDIR" -type f \( -iname '*lcv*' -o -iname '*.R' \) 2>/dev/null | head -n 1 || true)
if [ -n "$LCV_FOUND" ]; then record_ok LCV "$LCVDIR" source_repo; else record_missing LCV download_failed "https://github.com/lukejoconnor/LCV"; fi

CPDIR=$TOOLS/CPASSOC
mkdir -p "$CPDIR"
if [ ! -f "$CPDIR/CPASSOC.zip" ]; then
  (curl -L --fail --connect-timeout 30 -o "$CPDIR/CPASSOC.zip" "https://hal.case.edu/~xxz10/zhu-web/CPASSOC/CPASSOC.zip" || \
   curl -L --fail --connect-timeout 30 -o "$CPDIR/CPASSOC.zip" "http://hal.case.edu/~xxz10/zhu-web/CPASSOC/CPASSOC.zip" || true)
fi
if [ -s "$CPDIR/CPASSOC.zip" ]; then unzip -o "$CPDIR/CPASSOC.zip" -d "$CPDIR/unpacked" >/dev/null 2>&1 || true; fi
CP_FOUND=$(find "$CPDIR" -type f \( -iname '*cpassoc*' -o -iname '*.R' -o -perm -111 \) 2>/dev/null | head -n 1 || true)
if [ -n "$CP_FOUND" ]; then record_ok CPASSOC "$CP_FOUND" downloaded; else record_missing CPASSOC download_failed "CPASSOC.zip not found/unusable"; fi

conda list -p "$ENV" > "$CHECK/nasal_pairwise_env.conda_list.after_special_tools.txt" || true
Rscript - <<'RS' > "$CHECK/r_package_versions.after_special_tools.tsv"
pkgs <- c('TwoSampleMR','MRPRESSO','PLACO','data.table','coloc','susieR')
cat('package\tstatus\tversion\n')
for (p in pkgs) {
  ok <- suppressWarnings(requireNamespace(p, quietly=TRUE))
  if (ok) cat(p, 'OK', as.character(utils::packageVersion(p)), sep='\t') else cat(p, 'MISSING', 'NA', sep='\t')
  cat('\n')
}
RS

echo "STEP0C_DONE=$(date -Is)"
