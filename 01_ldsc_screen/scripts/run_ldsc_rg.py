#!/usr/bin/env python3
"""Run LDSC rg for all trait pairs and apply BH-FDR pair screening."""

from __future__ import annotations

import argparse
import csv
import math
import os
import re
import subprocess
from pathlib import Path


def read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.exists() or path.stat().st_size == 0:
        return []
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, str]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields, extrasaction="ignore", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def fnum(value: str | None) -> float | None:
    try:
        if value is None or str(value).strip() == "":
            return None
        parsed = float(value)
        return parsed if math.isfinite(parsed) else None
    except ValueError:
        return None


def safe_id(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", value)


def bh_adjust(p_values: list[float]) -> list[float]:
    """Benjamini-Hochberg adjusted P values in original order."""
    n = len(p_values)
    order = sorted(range(n), key=p_values.__getitem__, reverse=True)
    adjusted = [1.0] * n
    running = 1.0
    for index in order:
        rank = sum(p_values[j] <= p_values[index] for j in range(n))
        running = min(running, p_values[index] * n / rank)
        adjusted[index] = min(running, 1.0)
    return adjusted


def parse_rg_log(path: Path) -> dict[str, str]:
    out = {"rg": "", "se": "", "z": "", "p": ""}
    if not path.exists():
        return out
    text = path.read_text(errors="ignore")
    for line in text.splitlines():
        if "Genetic Correlation:" in line:
            vals = re.findall(r"[-+]?(?:\d*\.\d+|\d+)(?:[eE][-+]?\d+)?", line)
            if vals:
                out["rg"] = vals[0]
        elif line.strip().startswith("P:"):
            vals = re.findall(r"[-+]?(?:\d*\.\d+|\d+)(?:[eE][-+]?\d+)?", line)
            if vals:
                out["p"] = vals[0]
        elif "Z-score:" in line:
            vals = re.findall(r"[-+]?(?:\d*\.\d+|\d+)(?:[eE][-+]?\d+)?", line)
            if vals:
                out["z"] = vals[0]
    return out


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", default=os.environ.get("PROJECT_ROOT", ""))
    parser.add_argument("--ldsc-py", default=os.environ.get("LDSC_PY", "ldsc.py"))
    parser.add_argument("--pair-manifest", default="")
    args = parser.parse_args()

    if not args.project_root:
        raise SystemExit("Set --project-root or PROJECT_ROOT")
    root = Path(args.project_root).expanduser().resolve()
    out = root / "results/phase0_extension/step1_ldsc_screen"
    out.mkdir(parents=True, exist_ok=True)
    pair_manifest = Path(args.pair_manifest) if args.pair_manifest else root / "results/phase0_extension/step2_input_tables/pair_manifest.tsv"
    pairs = read_tsv(pair_manifest)
    if not pairs:
        raise SystemExit(f"No pairs found in {pair_manifest}")

    ref_ld = root / "data/reference/ldsc/1000G_Phase3_ldscores/LDscore"
    w_ld = root / "data/reference/ldsc/1000G_Phase3_weights_hm3_no_MHC/1000G_Phase3_weights_hm3_no_MHC"
    raw_rows: list[dict[str, str]] = []

    for row in pairs:
        pair_id = row.get("pair_id", "")
        s1 = row.get("trait1_munged", "")
        s2 = row.get("trait2_munged", "")
        result = {
            "pair_id": pair_id,
            "trait1": row.get("trait1", ""),
            "trait2": row.get("trait2", ""),
            "trait1_munged": s1,
            "trait2_munged": s2,
            "rg": "",
            "se": "",
            "z": "",
            "p": "",
            "status": "not_run",
            "log": "",
        }
        if not pair_id or not Path(s1).exists() or not Path(s2).exists():
            result["status"] = "missing_munged_sumstats"
            raw_rows.append(result)
            continue
        prefix = out / "logs" / safe_id(pair_id)
        prefix.parent.mkdir(parents=True, exist_ok=True)
        cmd = [
            args.ldsc_py,
            "--rg", f"{s1},{s2}",
            "--ref-ld-chr", str(ref_ld),
            "--w-ld-chr", str(w_ld),
            "--out", str(prefix),
        ]
        proc = subprocess.run(cmd, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        log_path = Path(str(prefix) + ".log")
        log_path.write_text(proc.stdout)
        result["log"] = str(log_path)
        result["status"] = "PASS" if proc.returncode == 0 else "FAIL"
        result.update(parse_rg_log(log_path))
        raw_rows.append(result)

    tested = [r for r in raw_rows if fnum(r.get("p")) is not None]
    adjusted = bh_adjust([fnum(r["p"]) for r in tested]) if tested else []
    q_by_pair = {row["pair_id"]: q for row, q in zip(tested, adjusted)}
    for row in raw_rows:
        q = q_by_pair.get(row["pair_id"])
        row["ldsc_bh_q"] = "" if q is None else f"{q:.12g}"
        rg = fnum(row.get("rg"))
        row["ldsc_bh_fdr_significant"] = str(
            q is not None and q < 0.05 and rg is not None and rg > 0
        )

    fields = [
        "pair_id", "trait1", "trait2", "trait1_munged", "trait2_munged",
        "rg", "se", "z", "p", "ldsc_bh_q",
        "ldsc_bh_fdr_significant", "status", "log",
    ]
    write_tsv(out / "ldsc_rg_summary.tsv", raw_rows, fields)
    write_tsv(
        out / "ldsc_rg_significant.tsv",
        [r for r in raw_rows if r["ldsc_bh_fdr_significant"] == "True"],
        fields,
    )
    write_tsv(out / "ldsc_qc.tsv", [{
        "n_pairs": str(len(raw_rows)),
        "n_tested": str(len(tested)),
        "fdr_threshold": "0.05",
        "n_significant": str(sum(r["ldsc_bh_fdr_significant"] == "True" for r in raw_rows)),
    }], ["n_pairs", "n_tested", "fdr_threshold", "n_significant"])


if __name__ == "__main__":
    main()
