#!/usr/bin/env bash
set -euo pipefail

NASAL_PROJECT_ROOT="${NASAL_PROJECT_ROOT:?Set NASAL_PROJECT_ROOT}"
ROOT="${RESULTS_ROOT:-${NASAL_PROJECT_ROOT}/results_end}/twas_fusion_gtexv8_magma_candidates"
WEIGHTS="${FUSION_WEIGHTS_DIR:-${NASAL_PROJECT_ROOT}/resources/fusion_gtexv8/weights/GTEx_v8}"
MANIFEST="${ROOT}/download_manifest.tsv"
LOG_DIR="${ROOT}/logs/downloads"
MAX_JOBS="${MAX_JOBS:-1}"

if [ "${MAX_JOBS}" -gt 2 ]; then
  MAX_JOBS=2
fi

mkdir -p "${WEIGHTS}/archives" "${LOG_DIR}"

download_one() {
  local tissue="$1"
  local url="$2"
  local archive="$3"
  local target="$4"
  local log="${LOG_DIR}/${tissue}.download.log"
  mkdir -p "$(dirname "${archive}")" "${target}"
  {
    date
    hostname
    free -h
    echo "tissue=${tissue}"
    echo "url=${url}"
  } > "${log}"
  wget -c -O "${archive}" "${url}" >> "${log}" 2>&1
  tar -xzf "${archive}" -C "${target}" >> "${log}" 2>&1
  {
    date
    free -h
  } >> "${log}"
}

export -f download_one
export LOG_DIR

if [ ! -s "${MANIFEST}" ]; then
  echo "Missing manifest: ${MANIFEST}" >&2
  exit 2
fi

tail -n +2 "${MANIFEST}" | awk -F'\t' 'NF >= 4 {print $1 "\t" $2 "\t" $3 "\t" $4}' \
  | xargs -r -P "${MAX_JOBS}" -n 4 bash -c 'download_one "$0" "$1" "$2" "$3"'
