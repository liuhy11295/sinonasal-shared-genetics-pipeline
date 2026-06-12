#!/usr/bin/env python3
"""Select significant local-rg rows using LAVA-adjusted P values."""

from __future__ import annotations

import argparse
import csv
import math
import os
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


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", default=os.environ.get("PROJECT_ROOT", ""))
    parser.add_argument("--lava-results", default="")
    args = parser.parse_args()

    if not args.project_root:
        raise SystemExit("Set --project-root or PROJECT_ROOT")
    root = Path(args.project_root).expanduser().resolve()
    out = root / "results/phase0_extension/step2_lava_screen"
    out.mkdir(parents=True, exist_ok=True)
    pair_manifest = read_tsv(root / "results/phase0_extension/step2_input_tables/pair_manifest.tsv")
    pair_ids = {r.get("pair_id", "") for r in pair_manifest if r.get("pair_id")}
    write_tsv(out / "lava_input_pairs.tsv", pair_manifest, ["pair_id", "trait1", "trait2", "trait1_munged", "trait2_munged"])

    lava_path = Path(args.lava_results) if args.lava_results else root / "results/phase0_server/nasal4_vs_other_package/figure_local_rg_edges.tsv"
    lava_rows = [
        r for r in read_tsv(lava_path)
        if not pair_ids or r.get("trait_pair", r.get("pair_id", "")) in pair_ids
    ]
    tested = [r for r in lava_rows if fnum(r.get("p_adj")) is not None]
    if lava_rows and not tested:
        raise SystemExit("LAVA result table must contain numeric p_adj values")
    for row in lava_rows:
        p_adj = fnum(row.get("p_adj"))
        local_rg = fnum(row.get("local_rg", row.get("value", row.get("rho"))))
        row["pair_id"] = row.get("trait_pair", row.get("pair_id", ""))
        row["lava_adjusted_significant"] = str(
            p_adj is not None
            and p_adj <= 0.05
            and local_rg is not None
            and local_rg > 0
        )

    fields = sorted(
        {k for row in lava_rows for k in row}
        | {"pair_id", "lava_adjusted_significant"}
    )
    write_tsv(out / "lava_local_rg_all.tsv", lava_rows, fields)
    write_tsv(
        out / "lava_local_rg_significant.tsv",
        [r for r in lava_rows if r.get("lava_adjusted_significant") == "True"],
        fields,
    )
    write_tsv(out / "lava_qc.tsv", [{
        "n_input_pairs": str(len(pair_ids)),
        "n_lava_tests": str(len(tested)),
        "adjusted_p_threshold": "0.05",
        "n_lava_positive_rows": str(sum(r.get("lava_adjusted_significant") == "True" for r in lava_rows)),
        "source_lava_results": str(lava_path),
    }], ["n_input_pairs", "n_lava_tests", "adjusted_p_threshold", "n_lava_positive_rows", "source_lava_results"])


if __name__ == "__main__":
    main()
