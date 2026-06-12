#!/usr/bin/env bash
set -euo pipefail

NASAL_PROJECT_ROOT="${NASAL_PROJECT_ROOT:?Set NASAL_PROJECT_ROOT}"
ROOT="${RESULTS_ROOT:-${NASAL_PROJECT_ROOT}/results_end}/twas_fusion_gtexv8_magma_candidates"
WEIGHTS="${FUSION_WEIGHTS_DIR:-${NASAL_PROJECT_ROOT}/resources/fusion_gtexv8/weights/GTEx_v8}"
LDREF_PREFIX="${FUSION_LDREF_PREFIX:?Set FUSION_LDREF_PREFIX}"
MAX_JOBS="${MAX_JOBS:-1}"
SWAP_USED_STOP_KB="${SWAP_USED_STOP_KB:-1782579}"
SWAP_CLEAR_MIN_USED_KB="${SWAP_CLEAR_MIN_USED_KB:-1}"
SWAP_CLEAR_SAFETY_KB="${SWAP_CLEAR_SAFETY_KB:-2097152}"

if [ "${MAX_JOBS}" -gt 1 ]; then
  echo "This 16G machine is restricted to MAX_JOBS=1 after trial RSS approached 10.8G." >&2
  MAX_JOBS=1
fi

mkdir -p "${ROOT}/audit" "${ROOT}/results/raw"

FAILED="${ROOT}/audit/failed_jobs.tsv"
MEMWARN="${ROOT}/audit/memory_warning.tsv"
JOB_AUDIT="${ROOT}/audit/fusion_job_plan.tsv"
SWAP_AUDIT="${ROOT}/audit/swap_activity.tsv"
SWAP_CLEANUP="${ROOT}/audit/swap_cleanup.tsv"

if [ ! -s "${FAILED}" ]; then
  echo -e "pair_id\ttissue\tchr\texit_status\treason\tstdout_log\tstderr_log\ttime_log" > "${FAILED}"
fi
if [ ! -s "${MEMWARN}" ]; then
  echo -e "pair_id\ttissue\tchr\tmax_rss_kb\ttime_log" > "${MEMWARN}"
fi
echo -e "pair_id\ttissue\tchr\tn_weight_rows\tstatus\toutput" > "${JOB_AUDIT}"
if [ ! -s "${SWAP_AUDIT}" ]; then
  echo -e "pair_id\ttissue\tchr\tmax_si\tmax_so\tswap_used_kb\tstatus" > "${SWAP_AUDIT}"
fi
if [ ! -s "${SWAP_CLEANUP}" ]; then
  echo -e "pair_id\ttissue\tchr\tswap_total_before_kb\tswap_used_before_kb\tmem_available_before_kb\taction\tswap_total_after_kb\tswap_used_after_kb\tstatus" > "${SWAP_CLEANUP}"
fi

while IFS=$'\t' read -r pair_id tissue pos_path n_rows; do
  if [ "${pair_id}" = "pair_id" ]; then
    continue
  fi
  if [ "${n_rows}" = "0" ]; then
    echo -e "${pair_id}\t${tissue}\tNA\t0\tskip_no_candidate_weights\t" >> "${JOB_AUDIT}"
    continue
  fi
  sumstats="${ROOT}/input/sumstats/${pair_id}.fusion.sumstats"
  if [ ! -s "${sumstats}" ]; then
    echo -e "${pair_id}\t${tissue}\tNA\t${n_rows}\tskip_missing_sumstats\t" >> "${JOB_AUDIT}"
    continue
  fi
  weights_dir="${WEIGHTS}/${tissue}/"
  if [ ! -d "${weights_dir}" ]; then
    echo -e "${pair_id}\t${tissue}\tNA\t${n_rows}\tskip_missing_weights_dir\t" >> "${JOB_AUDIT}"
    continue
  fi
  for chr in $(awk 'NR>1 {print $4}' "${pos_path}" | sort -n | uniq); do
    chr_n="$(awk -v c="${chr}" 'NR>1 && $4==c {n++} END{print n+0}' "${pos_path}")"
    if [ "${chr_n}" = "0" ]; then
      continue
    fi
    out="${ROOT}/results/raw/${pair_id}/${tissue}/chr${chr}.dat"
    if [ -s "${out}" ]; then
      echo -e "${pair_id}\t${tissue}\t${chr}\t${chr_n}\tskip_existing\t${out}" >> "${JOB_AUDIT}"
      continue
    fi
    echo -e "${pair_id}\t${tissue}\t${chr}\t${chr_n}\trun\t${out}" >> "${JOB_AUDIT}"
    "$(dirname "$0")/run_twas_job.sh" "${pair_id}" "${tissue}" "${chr}" "${sumstats}" "${pos_path}" "${weights_dir}" "${LDREF_PREFIX}" "${out}"
    last_warn="$(awk -F'\t' -v p="${pair_id}" -v t="${tissue}" -v c="${chr}" 'NR>1 && $1==p && $2==t && $3==c {print $4}' "${MEMWARN}" | tail -1)"
    if [ -n "${last_warn}" ]; then
      echo "Memory warning for ${pair_id} ${tissue} chr${chr}: ${last_warn} KB. Stopping batch." >&2
      exit 99
    fi
    vm_tmp="$(mktemp)"
    vmstat 1 4 > "${vm_tmp}"
    max_si="$(awk 'NR>2 {if($7>m)m=$7} END{print m+0}' "${vm_tmp}")"
    max_so="$(awk 'NR>2 {if($8>m)m=$8} END{print m+0}' "${vm_tmp}")"
    rm -f "${vm_tmp}"
    swap_total_kb="$(free -k | awk '$1=="Swap:" {print $2+0}')"
    swap_used_kb="$(free -k | awk '$1=="Swap:" {print $3+0}')"
    mem_available_kb="$(awk '$1=="MemAvailable:" {print $2+0}' /proc/meminfo)"
    cleanup_action="none"
    cleanup_status="not_needed"
    swap_total_after_kb="${swap_total_kb}"
    swap_used_after_kb="${swap_used_kb}"
    if [ "${swap_total_kb}" -gt 0 ] && [ "${swap_used_kb}" -ge "${SWAP_CLEAR_MIN_USED_KB}" ]; then
      cleanup_action="swapoff_swapon"
      if [ "${mem_available_kb}" -le $((swap_used_kb + SWAP_CLEAR_SAFETY_KB)) ]; then
        cleanup_status="insufficient_mem_available"
      elif sudo -n /usr/sbin/swapoff -a && sudo -n /usr/sbin/swapon -a; then
        cleanup_status="cleared"
        swap_total_after_kb="$(free -k | awk '$1=="Swap:" {print $2+0}')"
        swap_used_after_kb="$(free -k | awk '$1=="Swap:" {print $3+0}')"
      else
        cleanup_status="clear_failed"
        sudo -n /usr/sbin/swapon -a >/dev/null 2>&1 || true
        swap_total_after_kb="$(free -k | awk '$1=="Swap:" {print $2+0}')"
        swap_used_after_kb="$(free -k | awk '$1=="Swap:" {print $3+0}')"
      fi
    elif [ "${swap_total_kb}" -eq 0 ]; then
      cleanup_action="swapon_a"
      if sudo -n /usr/sbin/swapon -a; then
        cleanup_status="enabled_swap"
        swap_total_after_kb="$(free -k | awk '$1=="Swap:" {print $2+0}')"
        swap_used_after_kb="$(free -k | awk '$1=="Swap:" {print $3+0}')"
      else
        cleanup_status="enable_failed"
      fi
    fi
    echo -e "${pair_id}\t${tissue}\t${chr}\t${swap_total_kb}\t${swap_used_kb}\t${mem_available_kb}\t${cleanup_action}\t${swap_total_after_kb}\t${swap_used_after_kb}\t${cleanup_status}" >> "${SWAP_CLEANUP}"
    swap_used_kb="${swap_used_after_kb}"
    echo -e "${pair_id}\t${tissue}\t${chr}\t${max_si}\t${max_so}\t${swap_used_kb}\tok" >> "${SWAP_AUDIT}"
    if [ "${swap_used_kb}" -gt "${SWAP_USED_STOP_KB}" ]; then
      echo "Swap usage above threshold after ${pair_id} ${tissue} chr${chr}: used=${swap_used_kb} KB, threshold=${SWAP_USED_STOP_KB} KB. Stopping batch." >&2
      exit 98
    fi
    free -h
  done
done < "${ROOT}/audit/restricted_weight_pair_tissue_audit.tsv"
