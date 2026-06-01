#!/usr/bin/env python3
from __future__ import annotations

import os
import shutil
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path("/platform_data/p_user/p010/phase0")
EXCLUDE = "CHRONIC_RHINITIS_PANUKB_J31"
STAMP = datetime.now().strftime("%Y%m%d_%H%M%S")
ARCHIVE = ROOT / "results" / f"archive_chronic_rhinitis_excluded_{STAMP}"
STEP5 = ROOT / "results/phase0_extension/step5_placo_cpassoc"


def read(path: Path) -> pd.DataFrame:
    if not path.exists() or path.stat().st_size == 0:
        return pd.DataFrame()
    return pd.read_csv(path, sep="\t", dtype=str, low_memory=False)


def write(df: pd.DataFrame, path: Path):
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, sep="\t", index=False)


def backup(path: Path):
    if path.exists():
        rel = path.relative_to(ROOT)
        dest = ARCHIVE / rel
        dest.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(path, dest)


def is_excluded(df: pd.DataFrame) -> pd.Series:
    mask = pd.Series(False, index=df.index)
    for col in ["pair_id", "trait1", "trait2", "nasal_trait", "partner_trait"]:
        if col in df.columns:
            mask = mask | df[col].fillna("").astype(str).str.contains(EXCLUDE, regex=False)
    return mask


def bh(p):
    p = pd.to_numeric(p, errors="coerce")
    out = pd.Series(np.nan, index=p.index, dtype=float)
    mask = p.notna()
    vals = p[mask].astype(float)
    n = len(vals)
    if n == 0:
        return out
    order = vals.sort_values().index
    ranked = vals.loc[order].to_numpy()
    q = ranked * n / np.arange(1, n + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    q = np.clip(q, 0, 1)
    out.loc[order] = q
    return out


def bonf(p):
    p = pd.to_numeric(p, errors="coerce")
    n = int(p.notna().sum())
    return (p * n).clip(upper=1.0)


def filter_table(path: Path) -> tuple[int, int]:
    df = read(path)
    if df.empty:
        return (0, 0)
    backup(path)
    before = len(df)
    df = df.loc[~is_excluded(df)].copy()
    write(df, path)
    return before, len(df)


def step5():
    report = []
    placo_path = STEP5 / "placo_results.tsv"
    cpassoc_path = STEP5 / "cpassoc_results.tsv"
    sig_path = STEP5 / "cross_trait_significant_snps.tsv"
    cand_path = STEP5 / "candidate_snps.tsv"
    summary_path = STEP5 / "step5_summary.tsv"

    placo = read(placo_path)
    backup(placo_path)
    placo = placo.loc[~is_excluded(placo)].copy()
    if not placo.empty:
        placo["PLACO_bh_q"] = bh(placo["PLACO_p"])
        placo["PLACO_bonf_p"] = bonf(placo["PLACO_p"])
        placo["PLACO_bh_sig"] = pd.to_numeric(placo["PLACO_bh_q"], errors="coerce") <= 0.05
        placo["PLACO_bonf_sig"] = pd.to_numeric(placo["PLACO_bonf_p"], errors="coerce") <= 0.05
    write(placo, placo_path)
    report.append(("placo_results.tsv", len(placo), int(pd.to_numeric(placo.get("PLACO_bh_q"), errors="coerce").le(0.05).sum()) if not placo.empty else 0, int(pd.to_numeric(placo.get("PLACO_bonf_p"), errors="coerce").le(0.05).sum()) if not placo.empty else 0))

    cpassoc = read(cpassoc_path)
    backup(cpassoc_path)
    cpassoc = cpassoc.loc[~is_excluded(cpassoc)].copy()
    if not cpassoc.empty:
        for col in ["SHet_p", "SHom_p"]:
            prefix = col.replace("_p", "")
            cpassoc[f"{prefix}_bh_q"] = bh(cpassoc[col])
            cpassoc[f"{prefix}_bonf_p"] = bonf(cpassoc[col])
            cpassoc[f"{prefix}_bh_sig"] = pd.to_numeric(cpassoc[f"{prefix}_bh_q"], errors="coerce") <= 0.05
            cpassoc[f"{prefix}_bonf_sig"] = pd.to_numeric(cpassoc[f"{prefix}_bonf_p"], errors="coerce") <= 0.05
    write(cpassoc, cpassoc_path)
    cpassoc_bh = int((pd.to_numeric(cpassoc.get("SHet_bh_q"), errors="coerce").le(0.05) | pd.to_numeric(cpassoc.get("SHom_bh_q"), errors="coerce").le(0.05)).sum()) if not cpassoc.empty else 0
    cpassoc_bonf = int((pd.to_numeric(cpassoc.get("SHet_bonf_p"), errors="coerce").le(0.05) | pd.to_numeric(cpassoc.get("SHom_bonf_p"), errors="coerce").le(0.05)).sum()) if not cpassoc.empty else 0
    report.append(("cpassoc_results.tsv", len(cpassoc), cpassoc_bh, cpassoc_bonf))

    backup(sig_path)
    rows = []
    if not placo.empty:
        p = pd.to_numeric(placo["PLACO_bonf_p"], errors="coerce")
        for _, r in placo.loc[p <= 0.05].iterrows():
            rows.append({**{k: r.get(k, "") for k in ["pair_id", "trait1", "trait2", "SNP", "CHR", "BP", "EA", "OA"]}, "method": "PLACO", "p_value": r.get("PLACO_p", ""), "source": "PLACO_bonferroni_significant_SNP"})
    if not cpassoc.empty:
        for method, pcol, rawcol in [("CPASSOC_SHet", "SHet_bonf_p", "SHet_p"), ("CPASSOC_SHom", "SHom_bonf_p", "SHom_p")]:
            p = pd.to_numeric(cpassoc[pcol], errors="coerce")
            for _, r in cpassoc.loc[p <= 0.05].iterrows():
                rows.append({**{k: r.get(k, "") for k in ["pair_id", "trait1", "trait2", "SNP", "CHR", "BP", "EA", "OA"]}, "method": method, "p_value": r.get(rawcol, ""), "source": "CPASSOC_bonferroni_significant_SNP"})
    sig = pd.DataFrame(rows, columns=["pair_id", "trait1", "trait2", "SNP", "CHR", "BP", "EA", "OA", "method", "p_value", "source"]).drop_duplicates()
    write(sig, sig_path)

    cand = read(cand_path)
    backup(cand_path)
    if not cand.empty:
        cand = cand.loc[~is_excluded(cand)].copy()
        cand = cand.loc[~cand["source"].fillna("").str.contains("PLACO|CPASSOC", regex=True)].copy()
        add = []
        for _, r in sig.iterrows():
            src = "PLACO_significant_SNP" if r["method"] == "PLACO" else "CPASSOC_significant_SNP"
            add.append({"pair_id": r["pair_id"], "trait1": r["trait1"], "trait2": r["trait2"], "locus_id": "", "SNP": r["SNP"], "CHR": r["CHR"], "BP": r["BP"], "EA": r["EA"], "OA": r["OA"], "BETA_trait1": "", "SE_trait1": "", "P_trait1": "", "BETA_trait2": "", "SE_trait2": "", "P_trait2": "", "source": src})
        if add:
            cand = pd.concat([cand, pd.DataFrame(add)], ignore_index=True).drop_duplicates()
    write(cand, cand_path)

    backup(summary_path)
    summary = pd.DataFrame([{
        "placo_analyzed_pairs": placo["pair_id"].nunique() if not placo.empty else 0,
        "cpassoc_analyzed_pairs": cpassoc["pair_id"].nunique() if not cpassoc.empty else 0,
        "placo_significant_snps": int((sig["method"] == "PLACO").sum()) if not sig.empty else 0,
        "cpassoc_significant_snps": int(sig["method"].astype(str).str.startswith("CPASSOC").sum()) if not sig.empty else 0,
        "cross_trait_significant_snps": len(sig),
        "new_candidate_snps": len(sig),
        "failed_records": "",
        "correction": "Bonferroni<=0.05 after excluding CHRONIC_RHINITIS_PANUKB_J31",
    }])
    write(summary, summary_path)
    return report


def filter_step8():
    files = [
        ROOT / "results/step8_mtag/evidence_pair_manifest.tsv",
        ROOT / "results/step8_mtag/mtag_full_qc_report.tsv",
        ROOT / "results/step8_magma/pair_specific_gene_results.tsv",
        ROOT / "results/step8_magma/magma_qc_report.tsv",
        ROOT / "results/step8_magma/snp_evidence_tiers.tsv",
        ROOT / "results/step8_magma/gene_evidence_tiers.tsv",
        ROOT / "results/step8_pathway/pair_specific_gsea_go_bp.tsv",
        ROOT / "results/step8_pathway/pair_specific_gsea_kegg.tsv",
        ROOT / "results/step8_pathway/pair_specific_gsea_reactome.tsv",
        ROOT / "results/step8_pathway/pair_specific_ora_go_bp.tsv",
        ROOT / "results/step8_pathway/pair_specific_ora_kegg.tsv",
        ROOT / "results/step8_pathway/pair_specific_ora_reactome.tsv",
        ROOT / "results/step8_pathway/pathway_qc_report.tsv",
        ROOT / "results/step8_pathway/top3_pathways_by_pair_method.tsv",
    ]
    counts = []
    for path in files:
        if path.exists():
            counts.append((str(path.relative_to(ROOT)), *filter_table(path)))
    for d in [
        ROOT / "results/step8_mtag/full",
        ROOT / "results/step8_mtag/input",
        ROOT / "results/step8_magma/input",
        ROOT / "results/step8_magma/per_pair",
        ROOT / "results/step8_pathway/per_pair",
    ]:
        if not d.exists():
            continue
        for p in d.glob(f"*{EXCLUDE}*"):
            rel = p.relative_to(ROOT)
            dest = ARCHIVE / rel
            dest.parent.mkdir(parents=True, exist_ok=True)
            shutil.move(str(p), str(dest))
    return counts


def main():
    ARCHIVE.mkdir(parents=True, exist_ok=True)
    step5_report = step5()
    step8_counts = filter_step8()
    report = ARCHIVE / "exclusion_recorrection_report.tsv"
    with open(report, "w") as fh:
        fh.write("section\titem\trows_after\tbh_sig\tbonf_sig_or_rows_before\trows_after_filter\n")
        for item, rows, bh_sig, bonf_sig in step5_report:
            fh.write(f"step5\t{item}\t{rows}\t{bh_sig}\t{bonf_sig}\t\n")
        for item, before, after in step8_counts:
            fh.write(f"step8\t{item}\t\t\t{before}\t{after}\n")
    print(f"ARCHIVE={ARCHIVE}")
    print(f"REPORT={report}")


if __name__ == "__main__":
    main()
