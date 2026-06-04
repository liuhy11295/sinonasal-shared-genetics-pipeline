#!/usr/bin/env python3
"""Run MTAG for pairs positive in either parallel LDSC or parallel LAVA screening."""

from __future__ import annotations

import argparse
import csv
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


def safe_id(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", value)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", default=os.environ.get("PROJECT_ROOT", "/platform_data/p_user/p010/phase0"))
    parser.add_argument("--mtag-py", default=os.environ.get("MTAG_PY", "mtag.py"))
    parser.add_argument("--run", action="store_true", help="Actually execute MTAG. Without this flag only writes the selected pair table and commands.")
    args = parser.parse_args()

    root = Path(args.project_root)
    step1 = root / "results/phase0_extension/step1_ldsc_screen"
    step2_lava = root / "results/phase0_extension/step2_lava_screen"
    step2_inputs = root / "results/phase0_extension/step2_input_tables"
    out = root / "results/phase0_extension/step3_mtag_after_ldsc_lava"
    out.mkdir(parents=True, exist_ok=True)

    ldsc = {r.get("pair_id", ""): r for r in read_tsv(step1 / "ldsc_rg_significant.tsv") if r.get("pair_id")}
    lava = {r.get("pair_id", "") for r in read_tsv(step2_lava / "lava_local_rg_significant.tsv") if r.get("pair_id")}
    manifest = {r.get("pair_id", ""): r for r in read_tsv(step2_inputs / "pair_manifest.tsv") if r.get("pair_id")}
    selected_ids = sorted((set(ldsc) | lava) & set(manifest))

    selected = []
    commands = []
    for pair_id in selected_ids:
        row = manifest[pair_id]
        reasons = []
        if pair_id in ldsc:
            reasons.append("LDSC_bonferroni_positive")
        if pair_id in lava:
            reasons.append("LAVA_bonferroni_positive")
        out_prefix = out / "full" / safe_id(pair_id) / "mtag"
        cmd = [
            args.mtag_py,
            "--sumstats", f"{row.get('trait1_munged', '')},{row.get('trait2_munged', '')}",
            "--out", str(out_prefix),
            "--force",
        ]
        selected.append({
            "pair_id": pair_id,
            "trait1": row.get("trait1", ""),
            "trait2": row.get("trait2", ""),
            "trait1_munged": row.get("trait1_munged", ""),
            "trait2_munged": row.get("trait2_munged", ""),
            "selection_reason": ";".join(reasons),
            "mtag_out_prefix": str(out_prefix),
        })
        commands.append({"pair_id": pair_id, "command": " ".join(cmd)})
        if args.run:
            out_prefix.parent.mkdir(parents=True, exist_ok=True)
            proc = subprocess.run(cmd, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
            (out_prefix.parent / "mtag.log").write_text(proc.stdout)

    write_tsv(out / "mtag_input_pairs.tsv", selected, ["pair_id", "trait1", "trait2", "trait1_munged", "trait2_munged", "selection_reason", "mtag_out_prefix"])
    write_tsv(out / "mtag_commands.tsv", commands, ["pair_id", "command"])
    write_tsv(out / "mtag_selection_qc.tsv", [{
        "n_ldsc_positive_pairs": str(len(ldsc)),
        "n_lava_positive_pairs": str(len(lava)),
        "n_selected_pairs": str(len(selected_ids)),
        "selection_rule": "LDSC_bonferroni_positive OR LAVA_bonferroni_positive",
        "run_mtag": str(args.run),
    }], ["n_ldsc_positive_pairs", "n_lava_positive_pairs", "n_selected_pairs", "selection_rule", "run_mtag"])


if __name__ == "__main__":
    main()
