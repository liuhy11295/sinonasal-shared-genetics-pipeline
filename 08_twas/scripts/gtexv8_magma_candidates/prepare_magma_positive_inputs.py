#!/usr/bin/env python3
"""Prepare MAGMA Bonferroni-positive genes as candidate inputs for 49-tissue TWAS."""

from __future__ import annotations

import csv
from pathlib import Path


BASE = Path("/home/lhy/nasal/results_end")
OUT = BASE / "twas_fusion_gtexv8_magma_candidates"
MAGMA_SIG = BASE / "magma/significant_genes.tsv"
MAGMA_ALL = BASE / "magma/gene_results.tsv"


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


def is_true(value: str | None) -> bool:
    return str(value).strip().lower() in {"true", "t", "1", "yes", "y"}


def fnum(value: str | None) -> float | None:
    try:
        if value is None or str(value).strip() == "":
            return None
        return float(value)
    except ValueError:
        return None


def main() -> None:
    source = MAGMA_SIG if MAGMA_SIG.exists() else MAGMA_ALL
    rows = read_tsv(source)
    candidates = []
    for row in rows:
        bonf = fnum(row.get("gene_bonf_p"))
        significant = is_true(row.get("significant_gene")) or (bonf is not None and bonf <= 0.05)
        if not significant:
            continue
        candidates.append({
            "pair_id": row.get("pair_id", ""),
            "trait1": row.get("trait1", ""),
            "trait2": row.get("trait2", ""),
            "gene_symbol": row.get("gene_symbol", row.get("gene", "")),
            "gene_id": row.get("gene_id", ""),
            "chr": row.get("chr", row.get("CHR", "")),
            "gene_start": row.get("start", row.get("START", "")),
            "gene_end": row.get("stop", row.get("END", "")),
            "magma_p": row.get("p", row.get("P", "")),
            "magma_bonf_p": row.get("gene_bonf_p", ""),
            "magma_bonf_sig": "True",
            "magma_source": str(source),
        })
    fields = [
        "pair_id", "trait1", "trait2", "gene_symbol", "gene_id", "chr",
        "gene_start", "gene_end", "magma_p", "magma_bonf_p",
        "magma_bonf_sig", "magma_source",
    ]
    write_tsv(OUT / "input/magma_bonf_candidate_pair_genes.tsv", candidates, fields)
    unique_genes = sorted({r["gene_symbol"] for r in candidates if r["gene_symbol"]})
    (OUT / "input/magma_bonf_candidate_unique_genes.txt").write_text("\n".join(unique_genes) + "\n")
    write_tsv(OUT / "audit/magma_candidate_input_audit.tsv", [{
        "source": str(source),
        "n_pair_gene_records": str(len(candidates)),
        "n_unique_genes": str(len(unique_genes)),
        "n_pairs": str(len({r["pair_id"] for r in candidates if r["pair_id"]})),
        "rule": "MAGMA gene_bonf_p <= 0.05 or significant_gene=True",
    }], ["source", "n_pair_gene_records", "n_unique_genes", "n_pairs", "rule"])


if __name__ == "__main__":
    main()
