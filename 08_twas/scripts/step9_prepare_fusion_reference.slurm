#!/bin/bash
#SBATCH -J step9_twas_ref
#SBATCH -o %j.step9_twas_ref.out
#SBATCH -e %j.step9_twas_ref.err
#SBATCH -p compute_nodes
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem=32G
#SBATCH -t 0-06:00

set -euo pipefail

PROJECT_ROOT="${PROJECT_ROOT:?Set PROJECT_ROOT}"
cd "$PROJECT_ROOT"

REF_DIR=results_end/twas/reference
FUSION_DIR="${REF_DIR}/fusion_twas"
WEIGHT_DIR="${REF_DIR}/weights"
LDREF_DIR="${REF_DIR}/LDREF"
LOG_DIR=results_end/twas/logs
mkdir -p "$FUSION_DIR" "$WEIGHT_DIR" "$LDREF_DIR" "$LOG_DIR"

log() { echo "[$(date '+%F %T')] $*"; }

download() {
  local url="$1"
  local out="$2"
  if [ -s "$out" ]; then
    log "exists: $out"
    return 0
  fi
  log "download: $url"
  wget -O "${out}.tmp" "$url"
  mv "${out}.tmp" "$out"
}

log "Preparing FUSION TWAS reference under $REF_DIR"

if [ ! -s "${FUSION_DIR}/FUSION.assoc_test.R" ]; then
  if [ ! -d "${FUSION_DIR}/.git" ]; then
    rm -rf "${FUSION_DIR}.tmp"
    git clone --depth 1 https://github.com/gusevlab/fusion_twas.git "${FUSION_DIR}.tmp"
    cp -a "${FUSION_DIR}.tmp"/. "$FUSION_DIR"/
    rm -rf "${FUSION_DIR}.tmp"
  fi
fi

download \
  "https://data.broadinstitute.org/alkesgroup/FUSION/LDREF.tar.bz2" \
  "${REF_DIR}/LDREF.tar.bz2"
if [ ! -e "${LDREF_DIR}/1000G.EUR.1.bim" ] && [ ! -e "${REF_DIR}/LDREF/1000G.EUR.1.bim" ]; then
  tar -xjf "${REF_DIR}/LDREF.tar.bz2" -C "$REF_DIR"
fi

# GTEx v8 FUSION weight bundles. Keep original archives for provenance.
for tissue in \
  Lung \
  Whole_Blood \
  Spleen \
  Cells_EBV-transformed_lymphocytes \
  Skin_Sun_Exposed_Lower_leg \
  Skin_Not_Sun_Exposed_Suprapubic
do
  archive="${WEIGHT_DIR}/GTEx.${tissue}.tar.bz2"
  url="https://alkesgroup.broadinstitute.org/FUSION/WGT/GTEx.${tissue}.tar.bz2"
  download "$url" "$archive"
  if [ ! -s "${WEIGHT_DIR}/${tissue}.pos" ]; then
    mkdir -p "${WEIGHT_DIR}/${tissue}"
    tar -xjf "$archive" -C "${WEIGHT_DIR}/${tissue}"
    pos=$(find "${WEIGHT_DIR}/${tissue}" -type f -name '*.pos' | head -1)
    if [ -n "${pos:-}" ]; then
      ln -sfn "$pos" "${WEIGHT_DIR}/${tissue}.pos"
    fi
  fi
done

results_end/twas/scripts/check_twas_references.sh | tee "${LOG_DIR}/step9_reference_check_after_prepare.tsv"
log "done"
