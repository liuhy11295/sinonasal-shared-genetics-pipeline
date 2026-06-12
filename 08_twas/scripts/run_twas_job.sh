#!/usr/bin/env bash
set -euo pipefail

if [ "$#" -lt 8 ]; then
  echo "Usage: $0 PAIR_ID TISSUE CHR SUMSTATS WEIGHTS_POS WEIGHTS_DIR LDREF_PREFIX OUT_DAT" >&2
  exit 2
fi

PAIR_ID="$1"
TISSUE="$2"
CHR="$3"
SUMSTATS="$4"
WEIGHTS_POS="$5"
WEIGHTS_DIR="$6"
LDREF_PREFIX="$7"
OUT_DAT="$8"

NASAL_PROJECT_ROOT="${NASAL_PROJECT_ROOT:?Set NASAL_PROJECT_ROOT}"
ROOT="${RESULTS_ROOT:-${NASAL_PROJECT_ROOT}/results_end}/twas_fusion_gtexv8_magma_candidates"
FUSION_DIR="${FUSION_DIR:-${NASAL_PROJECT_ROOT}/resources/fusion_gtexv8/fusion_twas}"
FUSION="${FUSION_SCRIPT:-${FUSION_DIR}/FUSION.assoc_test.R}"
ENV_RUN="${FUSION_ENV_RUN:-conda run -n fusion_twas_gtexv8}"

LOG_DIR="${ROOT}/logs/${PAIR_ID}/${TISSUE}"
mkdir -p "${LOG_DIR}" "$(dirname "${OUT_DAT}")" "${ROOT}/audit"

STDOUT_LOG="${LOG_DIR}/chr${CHR}.stdout.log"
STDERR_LOG="${LOG_DIR}/chr${CHR}.stderr.log"
TIME_LOG="${LOG_DIR}/chr${CHR}.time.log"
PRE_MEM="${LOG_DIR}/chr${CHR}.pre_free_h.log"
POST_MEM="${LOG_DIR}/chr${CHR}.post_free_h.log"
FAILED="${ROOT}/audit/failed_jobs.tsv"
MEMWARN="${ROOT}/audit/memory_warning.tsv"

if [ ! -s "${FAILED}" ]; then
  echo -e "pair_id\ttissue\tchr\texit_status\treason\tstdout_log\tstderr_log\ttime_log" > "${FAILED}"
fi
if [ ! -s "${MEMWARN}" ]; then
  echo -e "pair_id\ttissue\tchr\tmax_rss_kb\ttime_log" > "${MEMWARN}"
fi

{
  date
  hostname
  free -h
} > "${PRE_MEM}"

set +e
(
cd "${FUSION_DIR}"
/usr/bin/time -v -o "${TIME_LOG}" ${ENV_RUN} Rscript "${FUSION}" \
  --sumstats "${SUMSTATS}" \
  --weights "${WEIGHTS_POS}" \
  --weights_dir "${WEIGHTS_DIR}" \
  --ref_ld_chr "${LDREF_PREFIX}" \
  --chr "${CHR}" \
  --out "${OUT_DAT}" \
  > "${STDOUT_LOG}" 2> "${STDERR_LOG}"
)
STATUS=$?
set -e

{
  date
  hostname
  free -h
} > "${POST_MEM}"

MAX_RSS="$(awk -F: '/Maximum resident set size/ {gsub(/^[ \t]+/, "", $2); print $2}' "${TIME_LOG}" 2>/dev/null | tail -1)"
if [ -n "${MAX_RSS}" ] && [ "${MAX_RSS}" -gt 12000000 ]; then
  echo -e "${PAIR_ID}\t${TISSUE}\t${CHR}\t${MAX_RSS}\t${TIME_LOG}" >> "${MEMWARN}"
fi

if [ "${STATUS}" -ne 0 ]; then
  REASON="exit_${STATUS}"
  if grep -Eqi "killed|cannot allocate memory|segfault|segmentation" "${STDERR_LOG}" "${STDOUT_LOG}"; then
    REASON="memory_or_crash"
  fi
  echo -e "${PAIR_ID}\t${TISSUE}\t${CHR}\t${STATUS}\t${REASON}\t${STDOUT_LOG}\t${STDERR_LOG}\t${TIME_LOG}" >> "${FAILED}"
  exit "${STATUS}"
fi

exit 0
