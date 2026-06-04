#!/usr/bin/env python3
import csv
from pathlib import Path

import pandas as pd
from statsmodels.stats.multitest import multipletests


ROOT = Path("/home/lhy/nasal/results_end/twas_fusion_gtexv8_magma_candidates")
RAW = ROOT / "results/raw"
TISSUE_OUT = ROOT / "results/tissue_level"
GENE_OUT = ROOT / "results/gene_level"
AUDIT = ROOT / "audit"


def read_candidates():
    return pd.read_csv(ROOT / "input/magma_bonf_candidate_pair_genes.tsv", sep="\t", dtype=str)


def read_raw_pair_tissue(pair_id, tissue):
    frames = []
    for path in sorted((RAW / pair_id / tissue).glob("chr*.dat")):
        if not path.exists() or path.stat().st_size == 0:
            continue
        try:
            df = pd.read_csv(path, sep=r"\s+|\t", engine="python", dtype=str)
        except Exception:
            continue
        if len(df):
            df["source_file"] = str(path)
            frames.append(df)
    if not frames:
        return pd.DataFrame()
    return pd.concat(frames, ignore_index=True)


def bh_fdr(series):
    p = pd.to_numeric(series, errors="coerce")
    out = pd.Series([pd.NA] * len(p), index=series.index, dtype="object")
    ok = p.notna()
    if ok.any():
        out.loc[ok] = multipletests(p.loc[ok].astype(float), method="fdr_bh")[1]
    return out


def normalize_gene_id(x):
    return str(x).split(".")[0] if pd.notna(x) else ""


def main():
    TISSUE_OUT.mkdir(parents=True, exist_ok=True)
    GENE_OUT.mkdir(parents=True, exist_ok=True)
    AUDIT.mkdir(parents=True, exist_ok=True)
    cand = read_candidates()
    cand["chr"] = cand["chr"].astype(str)

    restricted = sorted((ROOT / "input/restricted_weights").glob("*.magma_candidates.pos"))
    tissues = [p.name.replace(".magma_candidates.pos", "") for p in restricted]
    pair_ids = sorted(cand["pair_id"].unique())
    tissue_frames = []
    audit_rows = []

    for pair_id in pair_ids:
        cand_pair = cand[cand["pair_id"] == pair_id].copy()
        for tissue in tissues:
            raw = read_raw_pair_tissue(pair_id, tissue)
            n_requested = len(cand_pair)
            if raw.empty:
                audit_rows.append({
                    "pair_id": pair_id, "tissue": tissue,
                    "n_candidate_genes_requested": n_requested,
                    "n_genes_tested": 0, "n_genes_skipped": n_requested,
                    "n_TWAS_FDR_0.05": 0, "output_file_exists": "False",
                    "warnings_errors": "no_raw_output",
                })
                continue
            id_col = "ID" if "ID" in raw.columns else ("FILE" if "FILE" in raw.columns else None)
            if id_col is None or "TWAS.P" not in raw.columns:
                audit_rows.append({
                    "pair_id": pair_id, "tissue": tissue,
                    "n_candidate_genes_requested": n_requested,
                    "n_genes_tested": 0, "n_genes_skipped": n_requested,
                    "n_TWAS_FDR_0.05": 0, "output_file_exists": "True",
                    "warnings_errors": "missing_ID_or_TWAS.P",
                })
                continue
            raw["ensembl_id_base"] = raw[id_col].map(normalize_gene_id)
            map_audit = pd.read_csv(AUDIT / "weight_mapping_audit.tsv", sep="\t", dtype=str) if (AUDIT / "weight_mapping_audit.tsv").exists() else pd.DataFrame()
            if not map_audit.empty:
                raw = raw.merge(map_audit[["gene_symbol", "ensembl_id"]], left_on="ensembl_id_base", right_on="ensembl_id", how="left")
            raw["pair_id"] = pair_id
            raw["tissue"] = tissue
            raw["TWAS.FDR.within_pair_tissue"] = bh_fdr(raw["TWAS.P"])
            merged = raw.merge(
                cand_pair[["pair_id", "gene_symbol", "magma_p", "locus_grade"]].drop_duplicates(),
                on=["pair_id", "gene_symbol"],
                how="inner",
                suffixes=("", "_candidate"),
            )
            merged["TWAS.FDR.within_pair_tissue"] = bh_fdr(merged["TWAS.P"])
            tissue_path = TISSUE_OUT / f"{pair_id}.{tissue}.twas.tsv"
            merged.to_csv(tissue_path, sep="\t", index=False)
            tissue_frames.append(merged)
            tested = merged["gene_symbol"].nunique(dropna=True)
            sig = pd.to_numeric(merged["TWAS.FDR.within_pair_tissue"], errors="coerce").le(0.05).sum()
            audit_rows.append({
                "pair_id": pair_id, "tissue": tissue,
                "n_candidate_genes_requested": n_requested,
                "n_genes_tested": int(tested),
                "n_genes_skipped": int(max(n_requested - tested, 0)),
                "n_TWAS_FDR_0.05": int(sig),
                "output_file_exists": "True",
                "warnings_errors": "",
            })

    pd.DataFrame(audit_rows).to_csv(AUDIT / "twas_output_audit.tsv", sep="\t", index=False)
    if not tissue_frames:
        (GENE_OUT / "twas_candidate_gtexv8_49tissue_gene_level.tsv").write_text("")
        print("no_twas_outputs_found")
        return

    all_twas = pd.concat(tissue_frames, ignore_index=True)
    all_twas["TWAS.P.num"] = pd.to_numeric(all_twas["TWAS.P"], errors="coerce")
    all_twas["TWAS.FDR.num"] = pd.to_numeric(all_twas["TWAS.FDR.within_pair_tissue"], errors="coerce")

    gene_rows = []
    for (pair_id, gene_symbol), g in all_twas.dropna(subset=["gene_symbol"]).groupby(["pair_id", "gene_symbol"]):
        valid = g.dropna(subset=["TWAS.P.num"])
        if valid.empty:
            continue
        best = valid.loc[valid["TWAS.P.num"].idxmin()]
        sig_tissues = sorted(valid.loc[valid["TWAS.FDR.num"].le(0.05), "tissue"].astype(str).unique())
        gene_rows.append({
            "pair_id": pair_id,
            "gene_symbol": gene_symbol,
            "best_twas_p": best["TWAS.P"],
            "best_twas_tissue": best["tissue"],
            "best_twas_fdr": best["TWAS.FDR.within_pair_tissue"],
            "twas_supported_any_tissue": str(len(sig_tissues) > 0),
            "n_significant_tissues": len(sig_tissues),
            "significant_tissue_list": ",".join(sig_tissues),
            "magma_p": best.get("magma_p", ""),
            "locus_grade": best.get("locus_grade", ""),
        })
    gene_df = pd.DataFrame(gene_rows)
    gene_df.to_csv(GENE_OUT / "twas_candidate_gtexv8_49tissue_gene_level.tsv", sep="\t", index=False)

    strict = all_twas.copy()
    strict["strict_pairwise_all_tissues_FDR"] = pd.NA
    for pair_id, idx in strict.groupby("pair_id").groups.items():
        strict.loc[idx, "strict_pairwise_all_tissues_FDR"] = bh_fdr(strict.loc[idx, "TWAS.P"])
    strict["strict_multi_tissue_FDR_positive"] = pd.to_numeric(strict["strict_pairwise_all_tissues_FDR"], errors="coerce").le(0.05)
    strict.to_csv(GENE_OUT / "twas_candidate_gtexv8_strict_pairwise_all_tissues_FDR.tsv", sep="\t", index=False)

    fdr_audit = pd.DataFrame([{
        "main_fdr_rule": "BH-FDR within each pair_id x tissue over candidate genes tested by FUSION",
        "strict_fdr_rule": "BH-FDR within each pair_id over all candidate gene x tissue TWAS tests",
        "genome_wide_twas_fdr_used": "False",
    }])
    fdr_audit.to_csv(AUDIT / "fdr_audit.tsv", sep="\t", index=False)
    print(f"gene_level_records={len(gene_df)}")


if __name__ == "__main__":
    main()
