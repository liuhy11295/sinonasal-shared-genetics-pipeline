#!/usr/bin/env python3
from __future__ import annotations

import os
import re
from pathlib import Path


ROOT = Path(os.environ.get("PROJECT_ROOT", "/platform_data/p_user/p010/phase0"))
PAIR_ID = os.environ.get("PAIR_ID", "")
BENCH = ROOT / "results/resource_benchmark"
TIME_LOG = Path(os.environ.get("TIME_LOG", ""))


def parse_time_log(path: Path):
    text = path.read_text(errors="ignore") if path.exists() else ""
    peak_kb = None
    for line in text.splitlines():
        if "Maximum resident set size" in line:
            m = re.search(r"(\d+)", line)
            if m:
                peak_kb = int(m.group(1))
    return peak_kb


def main():
    BENCH.mkdir(parents=True, exist_ok=True)
    peak_kb = parse_time_log(TIME_LOG)
    peak_gb = (peak_kb or 0) / 1024 / 1024
    runtime_files = sorted(BENCH.glob(f"{PAIR_ID}.step8_runtime.tsv"))
    wall = ""
    n_snps = ""
    n_genes = ""
    if runtime_files:
        lines = runtime_files[-1].read_text().splitlines()
        if len(lines) > 1:
            parts = lines[1].split("\t")
            wall = parts[1]
            n_snps = parts[2]
            n_genes = parts[4]
    profile = BENCH / "resource_profile.tsv"
    profile.write_text(
        "step\tpair_id\tn_snps\tn_genes\twalltime_minutes\tpeak_memory_gb\tavg_memory_gb\tpeak_cpu_threads\n"
        f"MTAG_MAGMA_GSEA\t{PAIR_ID}\t{n_snps}\t{n_genes}\t{wall}\t{peak_gb:.3f}\tNA\t{os.environ.get('SLURM_CPUS_PER_TASK','')}\n"
    )
    safe = peak_gb * 1.5 if peak_gb else 512.0
    mem_jobs = int(512 // safe) if safe else 1
    cpu_jobs = 40
    final = max(1, min(mem_jobs, cpu_jobs))
    plan = BENCH / "parallelization_plan.tsv"
    plan.write_text(
        "peak_memory_per_job\tsafe_memory_per_job\tmemory_based_jobs\tcpu_based_jobs\tfinal_parallel_jobs\n"
        f"{peak_gb:.3f}\t{safe:.3f}\t{mem_jobs}\t{cpu_jobs}\t{final}\n"
    )
    print(f"RESOURCE_PROFILE_WRITTEN\t{profile}")
    print(f"PARALLELIZATION_PLAN_WRITTEN\t{plan}")


if __name__ == "__main__":
    main()
