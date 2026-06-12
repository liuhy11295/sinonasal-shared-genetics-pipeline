#!/usr/bin/env python3
"""Aggregate per-pair LAVA local-rg outputs using LAVA-adjusted P values."""

from __future__ import annotations

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
    project_root = os.environ.get("PROJECT_ROOT", "")
    if not project_root:
        raise SystemExit("Set PROJECT_ROOT")
    root = Path(project_root).expanduser().resolve()
    out = root / "results/phase0_extension/step2_lava_screen"
    rows: list[dict[str, str]] = []
    for path in sorted((out / "per_pair").glob("*.lava_local_rg.tsv")):
        rows.extend(read_tsv(path))
    tested = [r for r in rows if fnum(r.get("p_adj")) is not None]
    if rows and not tested:
        raise SystemExit("LAVA output must contain numeric p_adj values")
    for row in rows:
        p_adj = fnum(row.get("p_adj"))
        local_rg = fnum(row.get("local_rg", row.get("value", row.get("rho"))))
        row["pair_id"] = row.get("pair_id") or row.get("trait_pair", "")
        row["lava_adjusted_significant"] = str(
            p_adj is not None
            and p_adj <= 0.05
            and local_rg is not None
            and local_rg > 0
        )
    fields = sorted(
        {k for row in rows for k in row}
        | {"pair_id", "lava_adjusted_significant"}
    )
    write_tsv(out / "lava_local_rg_all.tsv", rows, fields)
    write_tsv(
        out / "lava_local_rg_significant.tsv",
        [r for r in rows if r.get("lava_adjusted_significant") == "True"],
        fields,
    )
    write_tsv(out / "lava_qc.tsv", [{
        "n_lava_tests": str(len(tested)),
        "adjusted_p_threshold": "0.05",
        "n_lava_positive_rows": str(sum(r.get("lava_adjusted_significant") == "True" for r in rows)),
    }], ["n_lava_tests", "adjusted_p_threshold", "n_lava_positive_rows"])


if __name__ == "__main__":
    main()
