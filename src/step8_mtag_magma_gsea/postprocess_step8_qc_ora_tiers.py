#!/usr/bin/env python3
from __future__ import annotations

import csv
import gzip
import math
import os
from pathlib import Path

import pandas as pd
from scipy.stats import hypergeom
from statsmodels.stats.multitest import multipletests


ROOT = Path(os.environ.get("PROJECT_ROOT", "/platform_data/p_user/p010/phase0"))
MTAG = ROOT / "results/step8_mtag"
MAGMA = ROOT / "results/step8_magma"
PATHWAY = ROOT / "results/step8_pathway"
GMT = ROOT / "results/phase0_extension/step8_magma_gene_pathway/resources/gmt"


def read_tsv(path, **kw):
    if not Path(path).exists() or Path(path).stat().st_size == 0:
        return pd.DataFrame()
    return pd.read_csv(path, sep="\t", dtype=str, low_memory=False, **kw)


def fnum(x):
    try:
        if pd.isna(x) or str(x).strip() == "":
            return math.nan
        return float(x)
    except Exception:
        return math.nan


def qvals(vals):
    arr = pd.to_numeric(vals, errors="coerce")
    out = pd.Series([math.nan] * len(arr))
    mask = arr.notna()
    if mask.any():
        out.loc[mask] = multipletests(arr.loc[mask], method="fdr_bh")[1]
    return out


def parse_gmt(path: Path):
    with open(path, errors="ignore") as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 3:
                yield parts[0], parts[1], set(x for x in parts[2:] if x)


def mtag_qc(manifest):
    rows = []
    for _, r in manifest.iterrows():
        pair = r["pair_id"]
        fgz = MTAG / "full" / pair / "full_mtag.tsv.gz"
        native_log = MTAG / "full" / pair / f"{pair}_full.log"
        n = 0
        pvals = []
        chrs = set()
        if fgz.exists():
            with gzip.open(fgz, "rt") as fh:
                reader = csv.DictReader(fh, delimiter="\t")
                for row in reader:
                    n += 1
                    if row.get("CHR"):
                        chrs.add(str(row["CHR"]).replace("chr", ""))
                    p = fnum(row.get("mtag_pval") or row.get("P"))
                    if math.isfinite(p):
                        pvals.append(p)
        arr = pd.Series(pvals, dtype="float64")
        log_text = native_log.read_text(errors="ignore") if native_log.exists() else ""
        rows.append(
            {
                "pair_id": pair,
                "mtag_file": str(fgz) if fgz.exists() else "",
                "n_snps": n,
                "min_p": arr.min() if len(arr) else "",
                "max_p": arr.max() if len(arr) else "",
                "median_p": arr.median() if len(arr) else "",
                "prop_p_lt_5e8": (arr < 5e-8).mean() if len(arr) else "",
                "prop_p_lt_1e5": (arr < 1e-5).mean() if len(arr) else "",
                "n_chr": len(chrs),
                "is_full_mtag": bool(n > 1_000_000 and len(chrs) >= 20 and (len(arr) and arr.max() > 0.9)),
                "need_rerun_mtag": False,
                "reason": "forced_low_mean_chi2" if "--force" in log_text and "mean chi2" in log_text.lower() else "",
            }
        )
    pd.DataFrame(rows).to_csv(MTAG / "mtag_full_qc_report.tsv", sep="\t", index=False)


def ora_for_db(gene_df, db, gmt_name, outname):
    rows = []
    for pair, g in gene_df.groupby("pair_id"):
        bg = set(g["gene_symbol"].dropna().astype(str))
        sig = set(g.loc[pd.to_numeric(g["gene_bonf_p"], errors="coerce") <= 0.05, "gene_symbol"].dropna().astype(str))
        trait1 = g["trait1"].iloc[0]
        trait2 = g["trait2"].iloc[0]
        M, N = len(bg), len(sig)
        for pid, desc, genes in parse_gmt(GMT / gmt_name):
            gs = genes & bg
            K = len(gs)
            x = len(gs & sig)
            if M and N and K and x:
                rows.append(
                    {
                        "pair_id": pair,
                        "trait1": trait1,
                        "trait2": trait2,
                        "database": db,
                        "pathway_id": pid,
                        "pathway_description": desc,
                        "foreground_size": N,
                        "background_size": M,
                        "overlap": x,
                        "pvalue": hypergeom.sf(x - 1, M, K, N),
                        "overlap_genes": "/".join(sorted(gs & sig)),
                    }
                )
    df = pd.DataFrame(rows)
    if not df.empty:
        df["p.adjust"] = qvals(df["pvalue"])
    df.to_csv(PATHWAY / outname, sep="\t", index=False)


def pathway_qc(gene_df):
    rows = []
    for pair, g in gene_df.groupby("pair_id"):
        rows.append(
            {
                "pair_id": pair,
                "tested_genes": len(g),
                "bonferroni_significant_genes": int((pd.to_numeric(g["gene_bonf_p"], errors="coerce") <= 0.05).sum()),
            }
        )
    qc = pd.DataFrame(rows)
    for name, path in [
        ("GO_BP", PATHWAY / "pair_specific_gsea_go_bp.tsv"),
        ("KEGG", PATHWAY / "pair_specific_gsea_kegg.tsv"),
        ("Reactome", PATHWAY / "pair_specific_gsea_reactome.tsv"),
    ]:
        df = read_tsv(path)
        sig = df[pd.to_numeric(df.get("p.adjust"), errors="coerce") <= 0.05] if not df.empty else pd.DataFrame()
        count = sig.groupby("pair_id").size().rename(f"{name}_fdr05_pathways")
        qc = qc.merge(count, left_on="pair_id", right_index=True, how="left")
    qc.fillna(0).to_csv(PATHWAY / "pathway_qc_report.tsv", sep="\t", index=False)


def magma_qc(gene_df):
    qc = (
        gene_df.assign(significant_gene=pd.to_numeric(gene_df["gene_bonf_p"], errors="coerce") <= 0.05)
        .groupby(["pair_id", "trait1", "trait2"], as_index=False)
        .agg(tested_genes=("gene_symbol", "count"), significant_genes=("significant_gene", "sum"))
    )
    qc["status"] = qc["tested_genes"].apply(lambda x: "PASS" if int(x) > 10000 else "WARN")
    qc.to_csv(MAGMA / "magma_qc_report.tsv", sep="\t", index=False)


def tiers(gene_df):
    high = read_tsv(ROOT / "results/phase0_extension/step3_susie_coloc/high_confidence_shared_variants.tsv")
    cross = read_tsv(ROOT / "results/phase0_extension/step5_placo_cpassoc/cross_trait_significant_snps.tsv")
    snp_rows = []
    for src, df in [("high_confidence_shared_variants", high), ("cross_trait_significant_snps", cross)]:
        if df.empty:
            continue
        for _, r in df.iterrows():
            methods = []
            if src.startswith("high"):
                methods.append("coloc_or_susie")
            if str(r.get("method", "")).upper() == "PLACO":
                methods.append("PLACO")
            if str(r.get("method", "")).upper().startswith("CPASSOC"):
                methods.append("CPASSOC")
            tier = "Level 1" if len(set(methods)) >= 2 or "coloc_or_susie" in methods else "Level 2"
            snp_rows.append(
                {
                    "pair_id": r.get("pair_id", ""),
                    "SNP": r.get("SNP", ""),
                    "CHR": r.get("CHR", ""),
                    "BP": r.get("BP", ""),
                    "support_methods": ";".join(sorted(set(methods))) or src,
                    "evidence_level": tier,
                    "source": src,
                }
            )
    pd.DataFrame(snp_rows).to_csv(MAGMA / "snp_evidence_tiers.tsv", sep="\t", index=False)

    gene_rows = []
    sig = gene_df[pd.to_numeric(gene_df["gene_bonf_p"], errors="coerce") <= 0.05]
    for _, r in sig.iterrows():
        gene_rows.append(
            {
                "pair_id": r["pair_id"],
                "gene_symbol": r["gene_symbol"],
                "gene_bonf_p": r["gene_bonf_p"],
                "support_methods": "MAGMA_full_MTAG",
                "evidence_level": "Level 2",
            }
        )
    pd.DataFrame(gene_rows).to_csv(MAGMA / "gene_evidence_tiers.tsv", sep="\t", index=False)


def readme(manifest, gene_df):
    mtag_qc_df = read_tsv(MTAG / "mtag_full_qc_report.tsv")
    magma_qc_df = read_tsv(MAGMA / "magma_qc_report.tsv")
    pathway_qc_df = read_tsv(PATHWAY / "pathway_qc_report.tsv")
    text = []
    text.append("# Step 8 full MTAG -> MAGMA -> GSEA\n")
    text.append(f"- Evidence pairs analyzed: {len(manifest)}")
    text.append("- LAVA inclusion rule: only the 10 pre-specified LAVA-driven coloc follow-up pairs are marked lava_sig.")
    text.append("- Full MTAG: all 33 pairs have full_mtag.tsv.gz outputs.")
    text.append("- MAGMA: official MAGMA v1.09a, 1000G EUR PLINK reference, NCBI37.3 gene locations, gene window +/-5 kb.")
    text.append("- Pathway GSEA: ranked enrichment over all MAGMA tested genes using GO BP, KEGG and Reactome GMT files.")
    text.append("- ORA: foreground is Bonferroni-significant MAGMA genes; background is the actual MAGMA tested genes per pair.")
    text.append("- If no pathway is significant, this means no pathway passed the current data/threshold, not absence of biology.")
    text.append("\n## Output counts\n")
    text.append(f"- MTAG QC rows: {len(mtag_qc_df)}")
    text.append(f"- MAGMA gene rows: {len(gene_df)}")
    text.append(f"- MAGMA QC rows: {len(magma_qc_df)}")
    text.append(f"- Pathway QC rows: {len(pathway_qc_df)}")
    text.append("\n## Important files\n")
    for p in [
        MTAG / "evidence_pair_manifest.tsv",
        MTAG / "mtag_full_qc_report.tsv",
        MAGMA / "pair_specific_gene_results.tsv",
        MAGMA / "magma_qc_report.tsv",
        MAGMA / "snp_evidence_tiers.tsv",
        MAGMA / "gene_evidence_tiers.tsv",
        PATHWAY / "pair_specific_gsea_go_bp.tsv",
        PATHWAY / "pair_specific_gsea_kegg.tsv",
        PATHWAY / "pair_specific_gsea_reactome.tsv",
        PATHWAY / "pair_specific_ora_go_bp.tsv",
        PATHWAY / "pair_specific_ora_kegg.tsv",
        PATHWAY / "pair_specific_ora_reactome.tsv",
        PATHWAY / "pathway_qc_report.tsv",
    ]:
        text.append(f"- {p}")
    (MAGMA / "README.md").write_text("\n".join(text) + "\n")
    (MAGMA / "method_notes.md").write_text(
        "# Method notes\n\n"
        "Full MTAG was rerun without p-value/locus output filters. Two chronic rhinitis pairs required MTAG `--force` because MTAG reported mean chi-square < 1.02; these pairs should be interpreted cautiously.\n\n"
        "MAGMA inputs used SNP, P and N columns from full trait-1/nasal MTAG output. Gene-level significance used Bonferroni p <= 0.05 within each pair.\n\n"
        "GSEA used all MAGMA tested genes ranked by MAGMA Z where available, otherwise by -log10(P). ORA used only Bonferroni-significant MAGMA genes as foreground and tested genes as background.\n"
    )


def main():
    manifest = read_tsv(MTAG / "evidence_pair_manifest.tsv")
    gene_df = read_tsv(MAGMA / "pair_specific_gene_results.tsv")
    if manifest.empty or gene_df.empty:
        raise SystemExit("missing manifest or pair_specific_gene_results.tsv")
    mtag_qc(manifest)
    magma_qc(gene_df)
    pathway_qc(gene_df)
    ora_for_db(gene_df, "GO_BP", "GO_Biological_Process_2023.gmt", "pair_specific_ora_go_bp.tsv")
    ora_for_db(gene_df, "KEGG", "KEGG_2021_Human.gmt", "pair_specific_ora_kegg.tsv")
    ora_for_db(gene_df, "Reactome", "Reactome_2022.gmt", "pair_specific_ora_reactome.tsv")
    tiers(gene_df)
    readme(manifest, gene_df)
    print("STEP8_POSTPROCESS_DONE")


if __name__ == "__main__":
    main()
