#!/usr/bin/env python3
from __future__ import annotations

import math
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from joblib import Parallel, delayed
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests

project_root = os.environ.get("PROJECT_ROOT", "").strip()
if not project_root:
    raise SystemExit("Set PROJECT_ROOT to the analysis project root.")
ROOT = Path(project_root)
TWAS = ROOT / "results_end/twas"
OUT = TWAS / "pathway_overlap_tissue_all31"
GMT = ROOT / "results/phase0_extension/step8_magma_gene_pathway/resources/gmt"

DBS = [
    ("GO_BP", "GO_Biological_Process_2023.gmt", ROOT / "results/pathway/gsea_go_bp.tsv"),
    ("KEGG", "KEGG_2021_Human.gmt", ROOT / "results/pathway/gsea_kegg.tsv"),
    ("Reactome", "Reactome_2022.gmt", ROOT / "results/pathway/gsea_reactome.tsv"),
]

FIELDS = [
    "pair_id",
    "trait1",
    "trait2",
    "tissue",
    "database",
    "pathway_id",
    "pathway_description",
    "setSize",
    "NES",
    "pvalue",
    "p.adjust",
    "qvalue",
    "core_enrichment",
]


def parse_gmt(path: Path):
    with path.open(encoding="utf-8", errors="ignore") as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            yield parts[0], parts[1], [g for g in parts[2:] if g]


def read_tasks():
    qc = pd.read_csv(TWAS / "twas_qc_summary.tsv", sep="\t", dtype=str).fillna("")
    qc = qc.loc[qc["status"].eq("PASS")].copy()
    qc = qc.sort_values(["pair_id", "tissue"])
    return qc[["pair_id", "tissue"]].drop_duplicates().to_dict("records")


def load_pair_tissue(pair_id: str, tissue: str):
    usecols = ["pair_id", "trait1", "trait2", "tissue", "gene_symbol", "TWAS.Z", "TWAS.P"]
    chunks = []
    for chunk in pd.read_csv(TWAS / "twas_all_results.tsv", sep="\t", dtype=str, usecols=usecols, chunksize=200000):
        sub = chunk.loc[(chunk["pair_id"] == pair_id) & (chunk["tissue"] == tissue)].copy()
        if not sub.empty:
            chunks.append(sub)
    if not chunks:
        return pd.DataFrame(columns=usecols)
    df = pd.concat(chunks, ignore_index=True).fillna("")
    df["P_num"] = pd.to_numeric(df["TWAS.P"], errors="coerce")
    df["Z_num"] = pd.to_numeric(df["TWAS.Z"], errors="coerce")
    df = df.sort_values("P_num").groupby("gene_symbol", as_index=False, dropna=False).head(1)
    return df


def test_gene_set(pair_id, trait1, trait2, tissue, db, ranks, universe, pathway_id, desc, genes):
    geneset = set(genes) & set(universe)
    if len(geneset) < 10 or len(geneset) > 500:
        return None
    set_vals = np.array([ranks[g] for g in geneset], dtype=float)
    bg_vals = np.array([v for g, v in ranks.items() if g not in geneset], dtype=float)
    if len(set_vals) == 0 or len(bg_vals) == 0:
        return None
    try:
        pval = float(mannwhitneyu(set_vals, bg_vals, alternative="two-sided").pvalue)
    except Exception:
        pval = 1.0
    all_vals = np.array(list(ranks.values()), dtype=float)
    sd = float(np.nanstd(all_vals)) or 1.0
    nes = float((np.nanmean(set_vals) - np.nanmean(all_vals)) / sd * math.sqrt(len(set_vals)))
    ordered = sorted(((g, ranks[g]) for g in geneset), key=lambda x: x[1], reverse=nes >= 0)
    return {
        "pair_id": pair_id,
        "trait1": trait1,
        "trait2": trait2,
        "tissue": tissue,
        "database": db,
        "pathway_id": pathway_id,
        "pathway_description": desc,
        "setSize": len(geneset),
        "NES": nes,
        "pvalue": pval,
        "core_enrichment": "/".join([g for g, _ in ordered[:200]]),
    }


def run_task(task_id: int):
    tasks = read_tasks()
    if task_id < 1 or task_id > len(tasks):
        raise SystemExit(f"task_id {task_id} outside 1..{len(tasks)}")
    task = tasks[task_id - 1]
    pair_id = task["pair_id"]
    tissue = task["tissue"]
    out_dir = OUT / "per_task"
    out_dir.mkdir(parents=True, exist_ok=True)
    outfile = out_dir / f"{task_id}.{pair_id}.{tissue}.twas_tissue_gsea.tsv"

    gene_df = load_pair_tissue(pair_id, tissue)
    if gene_df.empty:
        pd.DataFrame(columns=FIELDS).to_csv(outfile, sep="\t", index=False)
        return
    trait1 = gene_df["trait1"].iloc[0]
    trait2 = gene_df["trait2"].iloc[0]
    vals = gene_df["Z_num"]
    if vals.notna().sum() == 0:
        vals = -np.log10(gene_df["P_num"].clip(lower=1e-300))
    ranks = dict(zip(gene_df["gene_symbol"].astype(str), vals.fillna(0.0)))
    universe = list(ranks.keys())
    n_jobs = max(1, min(int(os.environ.get("SLURM_CPUS_PER_TASK", "1")), 8))

    frames = []
    for db, gmt_name, _ in DBS:
        terms = list(parse_gmt(GMT / gmt_name))
        rows = Parallel(n_jobs=n_jobs, backend="loky", batch_size=16)(
            delayed(test_gene_set)(pair_id, trait1, trait2, tissue, db, ranks, universe, pathway_id, desc, genes)
            for pathway_id, desc, genes in terms
        )
        rows = [r for r in rows if r is not None]
        if rows:
            df = pd.DataFrame(rows)
            q = multipletests(pd.to_numeric(df["pvalue"], errors="coerce").fillna(1.0), method="fdr_bh")[1]
            df["p.adjust"] = q
            df["qvalue"] = q
            frames.append(df[FIELDS])

    result = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame(columns=FIELDS)
    result.to_csv(outfile, sep="\t", index=False)


def aggregate():
    OUT.mkdir(parents=True, exist_ok=True)
    files = sorted((OUT / "per_task").glob("*.twas_tissue_gsea.tsv"))
    frames = []
    for path in files:
        try:
            df = pd.read_csv(path, sep="\t", dtype=str)
        except pd.errors.EmptyDataError:
            continue
        if not df.empty:
            frames.append(df)
    all_df = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame(columns=FIELDS)
    all_df.to_csv(OUT / "twas_tissue_gsea_all.tsv", sep="\t", index=False)

    sig = all_df.loc[pd.to_numeric(all_df.get("p.adjust"), errors="coerce") <= 0.05].copy()
    sig.to_csv(OUT / "twas_tissue_gsea_significant.tsv", sep="\t", index=False)

    overlaps = []
    for db, _, magma_path in DBS:
        if not magma_path.exists():
            continue
        magma = pd.read_csv(magma_path, sep="\t", dtype=str).fillna("")
        if "p.adjust" not in magma.columns:
            continue
        magma_sig = magma.loc[pd.to_numeric(magma["p.adjust"], errors="coerce") <= 0.05].copy()
        tdb = sig.loc[sig["database"].eq(db)].copy()
        if tdb.empty or magma_sig.empty:
            continue
        keep_cols = ["pair_id", "database", "pathway_id", "pathway_description", "NES", "pvalue", "p.adjust", "qvalue", "core_enrichment"]
        keep_cols = [c for c in keep_cols if c in magma_sig.columns]
        merged = tdb.merge(
            magma_sig[keep_cols],
            on=["pair_id", "database", "pathway_id"],
            how="inner",
            suffixes=("_twas", "_magma"),
        )
        if not merged.empty:
            overlaps.append(merged)
    ov = pd.concat(overlaps, ignore_index=True) if overlaps else pd.DataFrame()
    ov.to_csv(OUT / "twas_magma_tissue_gsea_overlap_by_tissue.tsv", sep="\t", index=False)

    if ov.empty:
        merged = pd.DataFrame(
            columns=[
                "pair_id",
                "trait1",
                "trait2",
                "database",
                "pathway_id",
                "pathway_description",
                "n_supporting_tissues",
                "supporting_tissues",
                "min_twas_fdr",
                "best_twas_p",
                "best_twas_nes",
                "magma_fdr",
                "magma_nes",
                "representative_core_enrichment",
            ]
        )
    else:
        ov["twas_fdr_num"] = pd.to_numeric(ov["p.adjust_twas"], errors="coerce")
        ov["twas_p_num"] = pd.to_numeric(ov["pvalue_twas"], errors="coerce")
        ov["abs_nes"] = pd.to_numeric(ov["NES_twas"], errors="coerce").abs()
        rows = []
        for keys, sub in ov.groupby(["pair_id", "trait1", "trait2", "database", "pathway_id"], dropna=False):
            best = sub.sort_values(["twas_fdr_num", "twas_p_num", "abs_nes"], ascending=[True, True, False]).iloc[0]
            pathway_desc = best.get("pathway_description_twas", "") or best.get("pathway_description_magma", "")
            rows.append(
                {
                    "pair_id": keys[0],
                    "trait1": keys[1],
                    "trait2": keys[2],
                    "database": keys[3],
                    "pathway_id": keys[4],
                    "pathway_description": pathway_desc,
                    "n_supporting_tissues": sub["tissue"].nunique(),
                    "supporting_tissues": ",".join(sorted(sub["tissue"].dropna().unique())),
                    "min_twas_fdr": best.get("p.adjust_twas", ""),
                    "best_twas_p": best.get("pvalue_twas", ""),
                    "best_twas_nes": best.get("NES_twas", ""),
                    "magma_fdr": best.get("p.adjust_magma", ""),
                    "magma_nes": best.get("NES_magma", ""),
                    "representative_core_enrichment": best.get("core_enrichment_twas", ""),
                }
            )
        merged = pd.DataFrame(rows).sort_values(["pair_id", "database", "min_twas_fdr"])
    merged.to_csv(OUT / "twas_magma_tissue_gsea_overlap_merged.tsv", sep="\t", index=False)

    if not sig.empty:
        sig_summary = (
            sig.groupby(["tissue", "database"], dropna=False)
            .agg(n_sig_rows=("pathway_id", "size"), n_pairs=("pair_id", "nunique"), n_pathways=("pathway_id", "nunique"))
            .reset_index()
        )
    else:
        sig_summary = pd.DataFrame(columns=["tissue", "database", "n_sig_rows", "n_pairs", "n_pathways"])
    sig_summary.to_csv(OUT / "twas_tissue_gsea_significant_summary.tsv", sep="\t", index=False)

    readme = OUT / "README.md"
    readme.write_text(
        "# TWAS tissue-specific GSEA for all 31 pairs\n\n"
        "This directory contains tissue-specific TWAS GSEA for all 31 disease pairs and six FUSION tissues. "
        "Each pair-tissue was ranked by TWAS.Z using all tested genes. BH FDR was calculated within each "
        "pair-tissue-database. `twas_magma_tissue_gsea_overlap_by_tissue.tsv` keeps tissue-specific overlaps "
        "with MAGMA GSEA positives. `twas_magma_tissue_gsea_overlap_merged.tsv` collapses tissue-specific "
        "overlaps into one row per pair/database/pathway with supporting tissues listed.\n",
        encoding="utf-8",
    )


def main():
    cmd = sys.argv[1] if len(sys.argv) > 1 else ""
    if cmd == "count":
        print(len(read_tasks()))
    elif cmd == "task":
        run_task(int(os.environ.get("SLURM_ARRAY_TASK_ID", "1")))
    elif cmd == "aggregate":
        aggregate()
    else:
        raise SystemExit("usage: phase0_twas_tissue_all_gsea_overlap.py count|task|aggregate")


if __name__ == "__main__":
    main()
