#!/usr/bin/env python3
"""Select the evidence pair with the largest munged-GWAS SNP overlap.

This is a lightweight preflight script. It reads only SNP identifiers from the
two munged files for each pair and writes an overlap table. It does not run any
analysis and does not write large intermediate files.
"""

from __future__ import annotations

import csv
import gzip
from pathlib import Path
from typing import Dict, Iterable, List


PROJECT_ROOT = Path(__file__).resolve().parents[3]
RESULTS = PROJECT_ROOT / "results"
PAIR_MANIFEST = RESULTS / "phase0_extension" / "step2_input_tables" / "pair_manifest.tsv"
EVIDENCE_MANIFEST = RESULTS / "step8_mtag" / "evidence_pair_manifest.tsv"
OUTDIR = RESULTS / "placo_genomewide" / "benchmark"
EXCLUDE = "CHRONIC_RHINITIS_PANUKB_J31"


def open_text(path: Path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="ignore", newline="")
    return path.open("r", encoding="utf-8", errors="ignore", newline="")


def read_tsv(path: Path) -> List[Dict[str, str]]:
    with path.open("r", encoding="utf-8", errors="ignore", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def count_snps(path: Path) -> tuple[int, set[str]]:
    snps: set[str] = set()
    with open_text(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames or "SNP" not in reader.fieldnames:
            raise ValueError(f"{path} lacks SNP column")
        for row in reader:
            snp = row.get("SNP", "")
            if snp:
                snps.add(snp)
    return len(snps), snps


def main() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    pairs = {
        row["pair_id"]: row
        for row in read_tsv(PAIR_MANIFEST)
        if row.get("pair_id") and EXCLUDE not in "\t".join(row.values())
    }
    evidence = [
        row for row in read_tsv(EVIDENCE_MANIFEST)
        if row.get("pair_id") in pairs and EXCLUDE not in "\t".join(row.values())
    ]
    rows = []
    cache: dict[Path, tuple[int, set[str]]] = {}
    for row in evidence:
        pair_id = row["pair_id"]
        manifest = pairs[pair_id]
        g1 = Path(manifest.get("trait1_munged", ""))
        g2 = Path(manifest.get("trait2_munged", ""))
        status = "PASS"
        reason = ""
        n1 = n2 = overlap = 0
        try:
            if not g1.exists() or not g2.exists():
                raise FileNotFoundError("missing munged GWAS file")
            if g1 not in cache:
                cache[g1] = count_snps(g1)
            if g2 not in cache:
                cache[g2] = count_snps(g2)
            n1, s1 = cache[g1]
            n2, s2 = cache[g2]
            overlap = len(s1 & s2)
        except Exception as exc:
            status = "FAIL"
            reason = str(exc)
        rows.append({
            "pair_id": pair_id,
            "trait1": row.get("trait1", manifest.get("trait1", "")),
            "trait2": row.get("trait2", manifest.get("trait2", "")),
            "gwas1_file": str(g1),
            "gwas2_file": str(g2),
            "n_snps_gwas1": n1,
            "n_snps_gwas2": n2,
            "n_snps_overlap": overlap,
            "status": status,
            "failed_reason": reason,
        })

    rows.sort(key=lambda item: int(item["n_snps_overlap"]), reverse=True)
    out = OUTDIR / "pair_overlap_benchmark_candidates.tsv"
    with out.open("w", encoding="utf-8", newline="") as handle:
        fields = [
            "pair_id", "trait1", "trait2", "gwas1_file", "gwas2_file",
            "n_snps_gwas1", "n_snps_gwas2", "n_snps_overlap",
            "status", "failed_reason",
        ]
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    if not rows or rows[0]["status"] != "PASS":
        raise SystemExit("No valid evidence pair with readable munged GWAS files")
    selected = OUTDIR / "benchmark_pair.tsv"
    with selected.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerow(rows[0])
    print(f"BENCHMARK_PAIR={rows[0]['pair_id']}")
    print(f"N_SNPS_OVERLAP={rows[0]['n_snps_overlap']}")


if __name__ == "__main__":
    main()
