#!/usr/bin/env python3
"""Aggregate coloc.abf locus outputs and hand off positive loci to coloc-SuSiE."""

from __future__ import annotations

import csv
import os
from pathlib import Path


project_root = os.environ.get("PROJECT_ROOT", "")
if not project_root:
    raise SystemExit("Set PROJECT_ROOT")
PROJECT_ROOT = Path(project_root).expanduser().resolve()
STEP2 = PROJECT_ROOT / "results/phase0_extension/step2_input_tables"
STEP3 = PROJECT_ROOT / "results/phase0_extension/step3_coloc_abf"
PER = STEP3 / "per_locus"
QC = STEP3 / "qc"
PPH4_THRESHOLD = 0.5
MAIN_THRESHOLD = 0.8


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
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fields})


def fnum(value: str | None) -> float | None:
    if value is None:
        return None
    text = str(value).strip()
    if not text or text.lower() in {"na", "nan", "none", "null"}:
        return None
    try:
        return float(text)
    except ValueError:
        return None


def main() -> None:
    result_fields = [
        "pair_id", "trait1", "trait2", "locus_id", "CHR", "START", "END",
        "nsnps", "PP.H0.abf", "PP.H1.abf", "PP.H2.abf", "PP.H3.abf", "PP.H4.abf",
        "lead_SNP", "lead_min_p", "status", "reason", "source",
    ]
    qc_fields = [
        "pair_id", "trait1", "trait2", "locus_id", "status", "reason",
        "n_trait1_region", "n_trait2_region", "n_harmonized_snps",
    ]
    positive_fields = ["pair_id", "trait1", "trait2", "locus_id", "CHR", "START", "END", "PP.H3", "PP.H4", "analysis_group"]
    candidate_fields = [
        "pair_id", "trait1", "trait2", "locus_id", "SNP", "CHR", "BP", "EA", "OA",
        "BETA_trait1", "SE_trait1", "P_trait1", "BETA_trait2", "SE_trait2", "P_trait2", "source",
    ]

    rows: list[dict[str, str]] = []
    for path in sorted(PER.glob("*.coloc_abf.tsv")):
        rows.extend(read_tsv(path))
    rows.sort(key=lambda row: (row.get("pair_id", ""), row.get("locus_id", "")))
    write_tsv(STEP3 / "coloc_abf_results.tsv", rows, result_fields)

    qc_rows: list[dict[str, str]] = []
    for path in sorted(QC.glob("*.qc.tsv")):
        qc_rows.extend(read_tsv(path))
    qc_rows.sort(key=lambda row: (row.get("pair_id", ""), row.get("locus_id", "")))
    write_tsv(STEP3 / "coloc_abf_qc.tsv", qc_rows, qc_fields)

    positives: list[dict[str, str]] = []
    for row in rows:
        pp4 = fnum(row.get("PP.H4.abf"))
        pp3 = fnum(row.get("PP.H3.abf"))
        if pp4 is None or pp4 < PPH4_THRESHOLD:
            continue
        positives.append({
            "pair_id": row.get("pair_id", ""),
            "trait1": row.get("trait1", ""),
            "trait2": row.get("trait2", ""),
            "locus_id": row.get("locus_id", ""),
            "CHR": row.get("CHR", ""),
            "START": row.get("START", ""),
            "END": row.get("END", ""),
            "PP.H3": row.get("PP.H3.abf", ""),
            "PP.H4": row.get("PP.H4.abf", ""),
            "analysis_group": (
                "main"
                if pp4 >= MAIN_THRESHOLD and (pp3 is None or pp4 > pp3)
                else "secondary"
            ),
        })
    write_tsv(STEP3 / "coloc_positive_loci.tsv", positives, positive_fields)

    # Existing coloc-SuSiE code reads this historical step2 path. Keep a synced
    # handoff file so the old SuSiE script can run without changing its internals.
    write_tsv(STEP2 / "coloc_positive_loci.tsv", positives, positive_fields)

    old_candidates = read_tsv(STEP2 / "candidate_snps.tsv")
    seen = {(r.get("pair_id", ""), r.get("locus_id", ""), r.get("SNP", ""), r.get("source", "")) for r in old_candidates}
    for row in rows:
        pp4 = fnum(row.get("PP.H4.abf"))
        snp = row.get("lead_SNP", "")
        key = (row.get("pair_id", ""), row.get("locus_id", ""), snp, "coloc_ABF_lead_SNP")
        if pp4 is None or pp4 < PPH4_THRESHOLD or not snp or key in seen:
            continue
        old_candidates.append({
            "pair_id": row.get("pair_id", ""),
            "trait1": row.get("trait1", ""),
            "trait2": row.get("trait2", ""),
            "locus_id": row.get("locus_id", ""),
            "SNP": snp,
            "CHR": row.get("CHR", ""),
            "BP": "",
            "EA": "",
            "OA": "",
            "BETA_trait1": "",
            "SE_trait1": "",
            "P_trait1": "",
            "BETA_trait2": "",
            "SE_trait2": "",
            "P_trait2": "",
            "source": "coloc_ABF_lead_SNP",
        })
        seen.add(key)
    write_tsv(STEP3 / "candidate_snps_with_coloc_abf.tsv", old_candidates, candidate_fields)
    write_tsv(STEP2 / "candidate_snps.tsv", old_candidates, candidate_fields)

    summary = [{
        "n_coloc_abf_loci": str(len(rows)),
        "n_coloc_abf_pass": str(sum(1 for r in rows if r.get("status") == "ok")),
        "n_coloc_h4_positive_loci": str(len(positives)),
        "main_positive_rule": "PP.H4 >= 0.80 and PP.H4 > PP.H3",
        "pph4_threshold": str(PPH4_THRESHOLD),
    }]
    write_tsv(STEP3 / "coloc_abf_summary.tsv", summary, ["n_coloc_abf_loci", "n_coloc_abf_pass", "n_coloc_h4_positive_loci", "main_positive_rule", "pph4_threshold"])


if __name__ == "__main__":
    main()
