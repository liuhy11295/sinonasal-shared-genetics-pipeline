#!/usr/bin/env python3
"""Aggregate genome-wide PLACO/CPASSOC sensitivity outputs.

This script does not rerun analyses. It combines per-pair significant-only
outputs and QC logs, updates integrated pair/SNP/locus summary tables with
genome-wide support flags, and writes inventory/readme/update log files.
"""

from __future__ import annotations

import csv
from datetime import datetime
from pathlib import Path
from typing import Dict, Iterable, List

import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[3]
RESULTS = PROJECT_ROOT / "results"
PLACO = RESULTS / "placo_genomewide"
CPASSOC = RESULTS / "cpassoc_genomewide"
INTEGRATED = RESULTS / "integrated"
QC_DIR = RESULTS / "qc"
EXCLUDE = "CHRONIC_RHINITIS_PANUKB_J31"
STEP = "genomewide_placo_cpassoc_sensitivity"


def read(path: Path) -> pd.DataFrame:
    if not path.exists() or path.stat().st_size == 0:
        return pd.DataFrame()
    return pd.read_csv(path, sep="\t", dtype=str, low_memory=False)


def write(df: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, sep="\t", index=False)


def numeric(series: pd.Series) -> pd.Series:
    return pd.to_numeric(series, errors="coerce")


def combine_tables(paths: Iterable[Path]) -> pd.DataFrame:
    frames = []
    for path in sorted(paths):
        df = read(path)
        if not df.empty:
            frames.append(df)
    return pd.concat(frames, ignore_index=True, sort=False) if frames else pd.DataFrame()


def filter_exclude(df: pd.DataFrame) -> pd.DataFrame:
    if df.empty:
        return df
    mask = pd.Series(False, index=df.index)
    for col in ["pair_id", "trait1", "trait2", "nasal_trait", "partner_trait"]:
        if col in df.columns:
            mask |= df[col].fillna("").astype(str).str.contains(EXCLUDE, regex=False)
    return df.loc[~mask].copy()


def method_support(df: pd.DataFrame, pcols: List[str], qcols: List[str]) -> pd.Series:
    if df.empty or "pair_id" not in df.columns:
        return pd.Series(dtype=bool)
    support = pd.Series(False, index=df.index)
    for col in pcols:
        if col in df.columns:
            support |= numeric(df[col]) < 5e-8
    for col in qcols:
        if col in df.columns:
            support |= numeric(df[col]) < 0.05
    return df.loc[support].groupby("pair_id").size().gt(0)


def count_by_pair(df: pd.DataFrame, key: str = "SNP") -> pd.DataFrame:
    if df.empty or "pair_id" not in df.columns:
        return pd.DataFrame(columns=["pair_id", "n"])
    count_col = key if key in df.columns else df.columns[0]
    return df.groupby("pair_id")[count_col].nunique().rename("n").reset_index()


def write_readmes() -> None:
    PLACO.mkdir(parents=True, exist_ok=True)
    CPASSOC.mkdir(parents=True, exist_ok=True)
    (PLACO / "README.md").write_text(
        "# Genome-wide PLACO sensitivity analysis\n\n"
        "Input: trait-pair munged GWAS summary statistics only. MTAG outputs and "
        "candidate SNP tables were not used as PLACO inputs.\n\n"
        "Only genome-wide significant rows are retained: PLACO_p < 5e-8. "
        "Full genome-wide PLACO result rows are not saved.\n",
        encoding="utf-8",
    )
    (CPASSOC / "README.md").write_text(
        "# Genome-wide CPASSOC sensitivity analysis\n\n"
        "Input: trait-pair munged GWAS summary statistics only. MTAG outputs and "
        "candidate SNP tables were not used as CPASSOC inputs.\n\n"
        "Only genome-wide significant rows are retained: SHet_p/SHom_p < 5e-8. "
        "Full genome-wide CPASSOC result rows are not saved.\n",
        encoding="utf-8",
    )


def inventory_row(path: Path, description: str) -> Dict[str, str]:
    rel = path.relative_to(RESULTS).as_posix()
    if path.suffix.lower() == ".tsv" and path.exists() and path.stat().st_size:
        try:
            with path.open("r", encoding="utf-8", errors="ignore") as handle:
                header = handle.readline().rstrip("\n").split("\t")
                n_rows = sum(1 for _ in handle)
            n_cols = len(header)
        except Exception:
            n_rows = ""
            n_cols = ""
    else:
        n_rows = ""
        n_cols = ""
    return {
        "file_path": rel,
        "file_type": path.suffix.lstrip(".") or "dir",
        "n_rows": n_rows,
        "n_cols": n_cols,
        "file_size": path.stat().st_size if path.exists() else 0,
        "description": description,
        "last_updated_step": STEP,
    }


def main() -> None:
    for directory in [PLACO / "significant", CPASSOC / "significant", INTEGRATED, QC_DIR]:
        directory.mkdir(parents=True, exist_ok=True)

    placo_sig = filter_exclude(combine_tables((PLACO / "significant").glob("*.placo.sig.tsv")))
    cpassoc_sig = filter_exclude(combine_tables((CPASSOC / "significant").glob("*.cpassoc.sig.tsv")))
    placo_qc = filter_exclude(combine_tables((PLACO / "logs").glob("*.placo.qc.tsv")))
    cpassoc_qc = filter_exclude(combine_tables((CPASSOC / "logs").glob("*.cpassoc.qc.tsv")))

    write(placo_sig, PLACO / "placo_genomewide_significant_snps.tsv")
    write(cpassoc_sig, CPASSOC / "cpassoc_genomewide_significant_snps.tsv")
    write(placo_qc, PLACO / "placo_qc.tsv")
    write(cpassoc_qc, CPASSOC / "cpassoc_qc.tsv")

    placo_summary = count_by_pair(placo_sig).rename(columns={"n": "placo_genomewide_snp_count"})
    cpassoc_summary = count_by_pair(cpassoc_sig).rename(columns={"n": "cpassoc_genomewide_snp_count"})
    write(placo_summary, PLACO / "placo_genomewide_summary.tsv")
    write(cpassoc_summary, CPASSOC / "cpassoc_genomewide_summary.tsv")

    master_path = INTEGRATED / "pair_master_summary.tsv"
    master = filter_exclude(read(master_path))
    if master.empty:
        evidence = filter_exclude(read(RESULTS / "step8_mtag" / "evidence_pair_manifest.tsv"))
        master = evidence[[c for c in ["pair_id", "trait1", "trait2"] if c in evidence.columns]].drop_duplicates()
    if "pair_id" not in master.columns:
        raise SystemExit("pair_master_summary.tsv lacks pair_id and cannot be updated")

    for col in ["placo_candidate_support", "cpassoc_candidate_support"]:
        if col not in master.columns:
            old_col = "placo_snp_count" if col.startswith("placo") else "cpassoc_snp_count"
            master[col] = numeric(master.get(old_col, pd.Series(0, index=master.index))).fillna(0).gt(0)

    placo_support = method_support(placo_sig, ["PLACO_p"], []).rename("placo_genomewide_support")
    cpassoc_support = method_support(cpassoc_sig, ["SHet_p", "SHom_p"], []).rename("cpassoc_genomewide_support")
    master = master.drop(columns=[c for c in ["placo_genomewide_support", "cpassoc_genomewide_support", "placo_support_type", "cpassoc_support_type"] if c in master.columns])
    master = master.merge(placo_support, left_on="pair_id", right_index=True, how="left")
    master = master.merge(cpassoc_support, left_on="pair_id", right_index=True, how="left")
    master["placo_genomewide_support"] = master["placo_genomewide_support"].fillna(False)
    master["cpassoc_genomewide_support"] = master["cpassoc_genomewide_support"].fillna(False)
    master["placo_support_type"] = master.apply(
        lambda r: "genomewide" if r["placo_genomewide_support"] else ("candidate_only" if bool(r.get("placo_candidate_support", False)) else "none"),
        axis=1,
    )
    master["cpassoc_support_type"] = master.apply(
        lambda r: "genomewide" if r["cpassoc_genomewide_support"] else ("candidate_only" if bool(r.get("cpassoc_candidate_support", False)) else "none"),
        axis=1,
    )
    write(master, master_path)

    snp_tables = []
    if not placo_sig.empty:
        tmp = placo_sig.copy()
        tmp["method"] = "PLACO_genomewide"
        tmp["p_value"] = tmp.get("PLACO_p", "")
        snp_tables.append(tmp)
    if not cpassoc_sig.empty:
        tmp = cpassoc_sig.copy()
        tmp["method"] = "CPASSOC_genomewide"
        if "SHet_p" in tmp.columns and "SHom_p" in tmp.columns:
            tmp["p_value"] = numeric(tmp["SHet_p"]).combine(numeric(tmp["SHom_p"]), min)
        snp_tables.append(tmp)
    snp_master = pd.concat(snp_tables, ignore_index=True, sort=False) if snp_tables else pd.DataFrame()
    write(snp_master, INTEGRATED / "snp_master_table.tsv")

    locus_cols = ["pair_id", "trait1", "trait2", "SNP", "CHR", "BP", "EA", "OA", "method", "p_value"]
    locus = snp_master[[c for c in locus_cols if c in snp_master.columns]].copy() if not snp_master.empty else pd.DataFrame(columns=locus_cols)
    if not locus.empty:
        locus["locus_id"] = locus["CHR"].astype(str) + ":" + locus["BP"].astype(str)
    write(locus, INTEGRATED / "locus_master_table.tsv")

    pair_counts = master[["pair_id", "trait1", "trait2", "placo_genomewide_support", "cpassoc_genomewide_support", "placo_candidate_support", "cpassoc_candidate_support", "placo_support_type", "cpassoc_support_type"]].copy()
    pair_counts = pair_counts.merge(placo_summary, on="pair_id", how="left").merge(cpassoc_summary, on="pair_id", how="left")
    for col in ["placo_genomewide_snp_count", "cpassoc_genomewide_snp_count"]:
        pair_counts[col] = pd.to_numeric(pair_counts[col], errors="coerce").fillna(0).astype(int)
    write(pair_counts, INTEGRATED / "pair_snp_locus_summary.tsv")

    global_qc = read(QC_DIR / "global_qc_report.tsv")
    remove = {"PLACO_genomewide_SNPs", "CPASSOC_genomewide_SNPs", "PLACO_genomewide_pairs", "CPASSOC_genomewide_pairs"}
    if not global_qc.empty and "metric" in global_qc.columns:
        global_qc = global_qc.loc[~global_qc["metric"].isin(remove)]
    appended = pd.DataFrame([
        {"metric": "PLACO_genomewide_SNPs", "value": placo_sig["SNP"].nunique() if "SNP" in placo_sig.columns else 0},
        {"metric": "CPASSOC_genomewide_SNPs", "value": cpassoc_sig["SNP"].nunique() if "SNP" in cpassoc_sig.columns else 0},
        {"metric": "PLACO_genomewide_pairs", "value": placo_sig["pair_id"].nunique() if "pair_id" in placo_sig.columns else 0},
        {"metric": "CPASSOC_genomewide_pairs", "value": cpassoc_sig["pair_id"].nunique() if "pair_id" in cpassoc_sig.columns else 0},
    ])
    global_qc = pd.concat([global_qc, appended], ignore_index=True, sort=False) if not global_qc.empty else appended
    write(global_qc, QC_DIR / "global_qc_report.tsv")

    write_readmes()
    (RESULTS / "README.md").write_text(
        "# Final project results\n\n"
        "This results directory distinguishes candidate-based PLACO/CPASSOC "
        "from genome-wide PLACO/CPASSOC sensitivity analyses. Genome-wide "
        "PLACO/CPASSOC used munged GWAS summary statistics as input and retained "
        "only significant SNP rows.\n",
        encoding="utf-8",
    )
    (RESULTS / "method_notes.md").write_text(
        "# Method notes\n\n"
        "Candidate-based PLACO/CPASSOC was run on variants selected from coloc, MTAG, "
        "and LAVA evidence. Genome-wide PLACO/CPASSOC sensitivity analysis was run "
        "from paired munged GWAS files, excluding candidate SNP tables and MTAG as inputs. "
        "Only P < 5e-8 rows were retained; full genome-wide outputs "
        "were not saved.\n",
        encoding="utf-8",
    )

    tracked = [
        (PLACO / "placo_genomewide_significant_snps.tsv", "Genome-wide PLACO significant SNPs only."),
        (PLACO / "placo_genomewide_summary.tsv", "Genome-wide PLACO pair-level summary."),
        (PLACO / "placo_qc.tsv", "Genome-wide PLACO QC."),
        (PLACO / "README.md", "Genome-wide PLACO README."),
        (CPASSOC / "cpassoc_genomewide_significant_snps.tsv", "Genome-wide CPASSOC significant SNPs only."),
        (CPASSOC / "cpassoc_genomewide_summary.tsv", "Genome-wide CPASSOC pair-level summary."),
        (CPASSOC / "cpassoc_qc.tsv", "Genome-wide CPASSOC QC."),
        (CPASSOC / "README.md", "Genome-wide CPASSOC README."),
        (INTEGRATED / "snp_master_table.tsv", "Integrated SNP master table without evidence tier scoring."),
        (INTEGRATED / "locus_master_table.tsv", "Integrated locus master table without evidence tier scoring."),
        (INTEGRATED / "pair_snp_locus_summary.tsv", "Pair-level SNP/locus support summary."),
        (INTEGRATED / "pair_master_summary.tsv", "Pair master summary with candidate and genome-wide support flags."),
        (QC_DIR / "global_qc_report.tsv", "Global QC with genome-wide PLACO/CPASSOC counts."),
        (RESULTS / "README.md", "Top-level results README."),
        (RESULTS / "method_notes.md", "Top-level method notes."),
    ]
    inventory = pd.DataFrame([inventory_row(path, desc) for path, desc in tracked if path.exists()])
    write(inventory, RESULTS / "final_project_inventory.tsv")
    write(inventory, RESULTS / "project_inventory.tsv")
    md_lines = ["# Project inventory", ""]
    for row in inventory.to_dict("records"):
        md_lines.append(
            f"- `{row['file_path']}`: {row['description']} "
            f"(rows={row['n_rows']}, cols={row['n_cols']}, size={row['file_size']})"
        )
    (RESULTS / "project_inventory.md").write_text("\n".join(md_lines) + "\n", encoding="utf-8")

    update_log = RESULTS / "update_log.md"
    update_log.write_text(
        "# Update log: genome-wide PLACO/CPASSOC sensitivity\n\n"
        f"Updated: {datetime.now().isoformat(timespec='seconds')}\n\n"
        "## Added/updated files\n\n"
        + "\n".join(f"- `{row['file_path']}`" for row in inventory.to_dict("records"))
        + "\n\n## Key statistics after update\n\n"
        f"- Genome-wide PLACO SNPs: {placo_sig['SNP'].nunique() if 'SNP' in placo_sig.columns else 0}\n"
        f"- Genome-wide CPASSOC SNPs: {cpassoc_sig['SNP'].nunique() if 'SNP' in cpassoc_sig.columns else 0}\n"
        f"- Genome-wide PLACO pairs: {placo_sig['pair_id'].nunique() if 'pair_id' in placo_sig.columns else 0}\n"
        f"- Genome-wide CPASSOC pairs: {cpassoc_sig['pair_id'].nunique() if 'pair_id' in cpassoc_sig.columns else 0}\n\n"
        "## Temporary full results\n\n"
        "- Full genome-wide PLACO/CPASSOC result rows were not saved.\n"
        "- Scratch directories are deleted by the per-pair SLURM script after extracting significant rows.\n\n"
        "## Evidence tier status\n\n"
        "- Evidence tiers, high-confidence tiers, and evidence scores were not generated or updated in this step.\n",
        encoding="utf-8",
    )
    print("AGGREGATE_GENOMEWIDE_PLACO_CPASSOC_DONE")


if __name__ == "__main__":
    main()
