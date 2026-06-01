#!/usr/bin/env python3
from __future__ import annotations

import csv
import gzip
import shutil
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path("/platform_data/p_user/p010/phase0")
RESULTS = ROOT / "results"
EXCLUDE = "CHRONIC_RHINITIS_PANUKB_J31"
OUT_DIRS = [
    "ldsc",
    "lava",
    "mtag",
    "placo",
    "cpassoc",
    "coloc",
    "susie_coloc",
    "magma",
    "pathway",
    "integrated",
    "visualization_inputs",
    "qc",
]


def read(path: str | Path) -> pd.DataFrame:
    p = Path(path)
    if not p.exists() or p.stat().st_size == 0:
        return pd.DataFrame()
    if str(p).endswith(".gz"):
        return pd.read_csv(p, sep="\t", dtype=str, low_memory=False, compression="gzip")
    return pd.read_csv(p, sep="\t", dtype=str, low_memory=False)


def write(df: pd.DataFrame, rel: str):
    path = RESULTS / rel
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, sep="\t", index=False)


def write_empty(rel: str, cols: list[str]):
    write(pd.DataFrame(columns=cols), rel)


def copy_table(src: str | Path, rel: str, filter_chronic: bool = True) -> pd.DataFrame:
    df = read(src)
    if filter_chronic and not df.empty:
        df = filter_out_chronic(df)
    write(df, rel)
    return df


def fnum(s):
    return pd.to_numeric(s, errors="coerce")


def filter_out_chronic(df: pd.DataFrame) -> pd.DataFrame:
    mask = pd.Series(False, index=df.index)
    for col in ["pair_id", "trait_pair", "trait1", "trait2", "trait_a", "trait_b", "nasal_trait", "partner_trait"]:
        if col in df.columns:
            mask = mask | df[col].fillna("").astype(str).str.contains(EXCLUDE, regex=False)
    return df.loc[~mask].copy()


def bh(p):
    p = pd.to_numeric(p, errors="coerce")
    out = pd.Series(np.nan, index=p.index, dtype=float)
    mask = p.notna()
    vals = p.loc[mask].astype(float)
    if vals.empty:
        return out
    order = vals.sort_values().index
    ranked = vals.loc[order].to_numpy()
    q = ranked * len(ranked) / np.arange(1, len(ranked) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    out.loc[order] = np.clip(q, 0, 1)
    return out


def simple_readme(folder: str, title: str, files: list[tuple[str, str]]):
    lines = [f"# {title}", "", f"Directory: `results/{folder}/`", "", "All tables are TSV unless noted.", ""]
    lines.append("慢性鼻炎 `CHRONIC_RHINITIS_PANUKB_J31` was excluded from final downstream paper-ready outputs because of low power; archived pre-exclusion files remain under `results/archive_chronic_rhinitis_excluded_*`.")
    lines.append("")
    lines.append("## Files")
    for name, desc in files:
        lines.append(f"- `{name}`: {desc}")
    (RESULTS / folder / "README.md").write_text("\n".join(lines) + "\n")


def row_count(path: Path) -> int:
    if not path.exists() or path.stat().st_size == 0:
        return 0
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt", errors="ignore") as fh:
        n = sum(1 for _ in fh)
    return max(n - 1, 0)


def ensure_dirs():
    for d in OUT_DIRS:
        (RESULTS / d).mkdir(parents=True, exist_ok=True)


def organize_ldsc():
    h2 = copy_table("results/phase0_server/nasal4_vs_other_package/nasal4_vs_other_ldsc_h2_summary.tsv", "ldsc/heritability_summary.tsv")
    rg = copy_table("results/phase0_server/nasal4_vs_other_package/nasal4_vs_other_ldsc_rg_summary.tsv", "ldsc/rg_summary.tsv")
    if not rg.empty:
        p = fnum(rg.get("p"))
        rg["ldsc_bh_q"] = bh(p)
        rg["ldsc_sig_bh05"] = rg["ldsc_bh_q"] <= 0.05
        sig = rg.loc[rg["ldsc_sig_bh05"]].copy()
        write(rg, "ldsc/rg_summary.tsv")
        write(sig, "ldsc/rg_significant.tsv")
        heat = rg.rename(columns={"trait_a": "source", "trait_b": "target", "rg": "value", "p": "pvalue"})
        write(heat[[c for c in ["pair_id", "source", "target", "value", "pvalue", "ldsc_bh_q"] if c in heat.columns]], "ldsc/ldsc_rg_heatmap_input.tsv")
        net = sig.rename(columns={"trait_a": "source", "trait_b": "target", "rg": "weight", "p": "pvalue"})
        write(net[[c for c in ["pair_id", "source", "target", "weight", "pvalue", "ldsc_bh_q"] if c in net.columns]], "ldsc/ldsc_rg_network_input.tsv")
        bub = rg.rename(columns={"trait_a": "nasal_trait", "trait_b": "partner_trait", "rg": "x", "p": "pvalue", "se": "size"})
        write(bub[[c for c in ["pair_id", "nasal_trait", "partner_trait", "x", "size", "pvalue", "ldsc_bh_q"] if c in bub.columns]], "ldsc/ldsc_rg_bubble_input.tsv")
    else:
        write_empty("ldsc/rg_significant.tsv", [])
    qc = pd.DataFrame(
        [
            {"metric": "heritability_rows", "value": len(h2)},
            {"metric": "rg_rows", "value": len(rg)},
            {"metric": "rg_bh05_rows", "value": int((rg.get("ldsc_sig_bh05", pd.Series(dtype=bool)) == True).sum()) if not rg.empty else 0},
        ]
    )
    write(qc, "ldsc/ldsc_qc.tsv")
    simple_readme("ldsc", "LDSC Results", [
        ("heritability_summary.tsv", "Trait-level LDSC heritability summary."),
        ("rg_summary.tsv", "All LDSC genetic correlation results with BH q-values."),
        ("rg_significant.tsv", "BH q <= 0.05 LDSC pairs."),
        ("ldsc_qc.tsv", "Row counts and LDSC result QC."),
    ])


def organize_lava():
    local = copy_table("results/phase0_server/nasal4_vs_other_package/figure_local_rg_edges.tsv", "lava/local_rg_all.tsv")
    if not local.empty:
        pcol = "p_adj" if "p_adj" in local.columns else "p"
        local["lava_significant"] = fnum(local[pcol]) <= 0.05
        write(local, "lava/local_rg_all.tsv")
        sig = local.loc[local["lava_significant"]].copy()
        write(sig, "lava/local_rg_significant.tsv")
        pair_summary = local.groupby("trait_pair", as_index=False).agg(
            n_loci=("region_id", "count"),
            n_significant_loci=("lava_significant", "sum"),
        ).rename(columns={"trait_pair": "pair_id"})
        write(pair_summary, "lava/local_rg_pair_summary.tsv")
        write(local.rename(columns={"trait_pair": "pair_id", "local_rg": "value"}), "lava/lava_local_rg_heatmap_input.tsv")
        write(local, "lava/lava_local_rg_track_input.tsv")
        write(sig.rename(columns={"trait_pair": "pair_id", "local_rg": "weight"}), "lava/lava_pair_network_input.tsv")
    else:
        for rel in ["local_rg_significant.tsv", "local_rg_pair_summary.tsv", "lava_local_rg_heatmap_input.tsv", "lava_local_rg_track_input.tsv", "lava_pair_network_input.tsv"]:
            write_empty(f"lava/{rel}", [])
    simple_readme("lava", "LAVA Local Genetic Correlation Results", [
        ("local_rg_all.tsv", "All LAVA local rg rows available for final reporting."),
        ("local_rg_significant.tsv", "Local rg rows with adjusted p <= 0.05 where available."),
        ("local_rg_pair_summary.tsv", "Per-pair LAVA locus counts."),
    ])


def organize_mtag():
    qc = copy_table("results/step8_mtag/mtag_full_qc_report.tsv", "mtag/pair_mtag_qc.tsv")
    manifest = copy_table("results/step8_mtag/evidence_pair_manifest.tsv", "mtag/full_mtag_manifest.tsv")
    if not manifest.empty:
        rows = []
        for _, r in manifest.iterrows():
            pair = r["pair_id"]
            fgz = RESULTS / "step8_mtag/full" / pair / "full_mtag.tsv.gz"
            rows.append({"pair_id": pair, "trait1": r.get("trait1", ""), "trait2": r.get("trait2", ""), "full_mtag_file": str(fgz), "file_exists": fgz.exists(), "size_bytes": fgz.stat().st_size if fgz.exists() else 0})
        summary = pd.DataFrame(rows)
        write(summary, "mtag/pair_mtag_summary.tsv")
        existing = read("results/phase0_server/nasal4_vs_other_package/figure_mtag_manhattan.tsv")
        existing = filter_out_chronic(existing) if not existing.empty else existing
        if not existing.empty:
            rename = {}
            if "trait_pair" in existing.columns and "pair_id" not in existing.columns:
                rename["trait_pair"] = "pair_id"
            if "rsid" in existing.columns and "SNP" not in existing.columns:
                rename["rsid"] = "SNP"
            if "chr" in existing.columns and "CHR" not in existing.columns:
                rename["chr"] = "CHR"
            if "pos" in existing.columns and "BP" not in existing.columns:
                rename["pos"] = "BP"
            existing = existing.rename(columns=rename)
            pcol = (
                "mtag_pval" if "mtag_pval" in existing.columns
                else ("p_mtag" if "p_mtag" in existing.columns else ("P" if "P" in existing.columns else None))
            )
            if pcol:
                sig = existing.loc[fnum(existing[pcol]) <= 5e-8].copy()
                ranked = existing.assign(_p=fnum(existing[pcol]))
                if "pair_id" in ranked.columns:
                    top = ranked.sort_values(["pair_id", "_p"], na_position="last").groupby("pair_id").head(200)
                else:
                    top = ranked.sort_values("_p", na_position="last").head(5000)
            else:
                sig = existing.copy()
                top = existing.copy()
            write(sig, "mtag/significant_snps.tsv")
            write(top, "mtag/mtag_top_loci_input.tsv")
            write(existing, "mtag/mtag_manhattan_input.tsv")
        else:
            write_empty("mtag/significant_snps.tsv", [])
            write_empty("mtag/mtag_top_loci_input.tsv", [])
            write_empty("mtag/mtag_manhattan_input.tsv", [])
        write(summary, "mtag/mtag_pair_summary_input.tsv")
    else:
        for rel in ["pair_mtag_summary.tsv", "significant_snps.tsv", "mtag_manhattan_input.tsv", "mtag_top_loci_input.tsv", "mtag_pair_summary_input.tsv"]:
            write_empty(f"mtag/{rel}", [])
    simple_readme("mtag", "Full MTAG Results", [
        ("full_mtag_manifest.tsv", "Manifest of full MTAG output files; full SNP-level files remain in step8_mtag/full to avoid duplication."),
        ("significant_snps.tsv", "Genome-wide significant MTAG SNPs extracted from full MTAG."),
        ("pair_mtag_qc.tsv", "MTAG QC summary."),
    ])


def organize_placo_cpassoc():
    placo = copy_table("results/phase0_extension/step5_placo_cpassoc/placo_results.tsv", "placo/source_core_placo_results.tsv")
    if not placo.empty:
        sig = placo.loc[fnum(placo.get("PLACO_bonf_p")) <= 0.05].copy()
        write(sig, "placo/placo_significant_snps.tsv")
        summ = sig.groupby("pair_id", as_index=False).agg(placo_snp_count=("SNP", "nunique"))
        write(summ, "placo/placo_pair_summary.tsv")
        write(pd.DataFrame([{"metric": "rows", "value": len(placo)}, {"metric": "bonf_significant_rows", "value": len(sig)}]), "placo/placo_qc.tsv")
        write(sig, "placo/placo_manhattan_input.tsv")
        write(summ, "placo/placo_pair_count_input.tsv")
    else:
        for rel in ["placo_significant_snps.tsv", "placo_pair_summary.tsv", "placo_qc.tsv", "placo_manhattan_input.tsv", "placo_pair_count_input.tsv"]:
            write_empty(f"placo/{rel}", [])
    simple_readme("placo", "PLACO Results", [("source_core_placo_results.tsv", "All PLACO candidate SNP results with recalculated correction columns."), ("placo_significant_snps.tsv", "Bonferroni-significant PLACO SNPs.")])

    cp = copy_table("results/phase0_extension/step5_placo_cpassoc/cpassoc_results.tsv", "cpassoc/source_core_cpassoc_results.tsv")
    if not cp.empty:
        mask = (fnum(cp.get("SHet_bonf_p")) <= 0.05) | (fnum(cp.get("SHom_bonf_p")) <= 0.05)
        sig = cp.loc[mask].copy()
        write(sig, "cpassoc/cpassoc_significant_snps.tsv")
        summ = sig.groupby("pair_id", as_index=False).agg(cpassoc_snp_count=("SNP", "nunique"))
        write(summ, "cpassoc/cpassoc_pair_summary.tsv")
        write(pd.DataFrame([{"metric": "rows", "value": len(cp)}, {"metric": "bonf_significant_rows", "value": len(sig)}]), "cpassoc/cpassoc_qc.tsv")
        write(sig, "cpassoc/cpassoc_manhattan_input.tsv")
        write(summ, "cpassoc/cpassoc_pair_count_input.tsv")
    else:
        for rel in ["cpassoc_significant_snps.tsv", "cpassoc_pair_summary.tsv", "cpassoc_qc.tsv", "cpassoc_manhattan_input.tsv", "cpassoc_pair_count_input.tsv"]:
            write_empty(f"cpassoc/{rel}", [])
    simple_readme("cpassoc", "CPASSOC Results", [("source_core_cpassoc_results.tsv", "All CPASSOC candidate SNP results with recalculated correction columns."), ("cpassoc_significant_snps.tsv", "Bonferroni-significant CPASSOC SNPs.")])


def organize_coloc_susie():
    coloc = copy_table("results/phase0_extension/step2_input_tables/coloc_positive_loci.tsv", "coloc/coloc_positive_loci.tsv")
    if not coloc.empty:
        summ = coloc.groupby("pair_id", as_index=False).agg(coloc_loci_count=("locus_id", "nunique"), max_pph4=("PP.H4", lambda x: fnum(x).max()))
        write(summ, "coloc/coloc_pair_summary.tsv")
        write(pd.DataFrame([{"metric": "positive_loci", "value": len(coloc)}, {"metric": "pairs", "value": coloc["pair_id"].nunique()}]), "coloc/coloc_qc.tsv")
        write(coloc, "coloc/coloc_locus_plot_input.tsv")
        write(summ, "coloc/coloc_pair_count_input.tsv")
    else:
        for rel in ["coloc_pair_summary.tsv", "coloc_qc.tsv", "coloc_locus_plot_input.tsv", "coloc_pair_count_input.tsv"]:
            write_empty(f"coloc/{rel}", [])
    simple_readme("coloc", "Coloc Results", [("coloc_positive_loci.tsv", "Positive coloc loci PP.H4 >= 0.5 from existing results."), ("coloc_pair_summary.tsv", "Per-pair coloc locus counts.")])

    susie = copy_table("results/phase0_extension/step3_susie_coloc/coloc_susie_results.tsv", "susie_coloc/source_core_coloc_susie_results.tsv")
    qc = copy_table("results/phase0_extension/step3_susie_coloc/susie_qc_report.tsv", "susie_coloc/susie_qc.tsv")
    if not susie.empty:
        pos = susie.loc[fnum(susie.get("PP.H4")) >= 0.5].copy()
        write(pos, "susie_coloc/susie_positive_loci.tsv")
        summ = pos.groupby("pair_id", as_index=False).agg(susie_loci_count=("locus_id", "nunique"), susie_shared_signals=("shared_signal_id", "nunique"))
        write(summ, "susie_coloc/susie_pair_summary.tsv")
        write(pos, "susie_coloc/susie_locus_plot_input.tsv")
        write(summ, "susie_coloc/susie_pair_count_input.tsv")
    else:
        for rel in ["susie_positive_loci.tsv", "susie_pair_summary.tsv", "susie_locus_plot_input.tsv", "susie_pair_count_input.tsv"]:
            write_empty(f"susie_coloc/{rel}", [])
    simple_readme("susie_coloc", "SuSiE Coloc Results", [("source_core_coloc_susie_results.tsv", "Existing coloc.susie shared signal table."), ("susie_positive_loci.tsv", "PP.H4 >= 0.5 SuSiE coloc rows.")])


def organize_magma_pathway():
    genes = copy_table("results/step8_magma/pair_specific_gene_results.tsv", "magma/gene_results.tsv")
    tiers = copy_table("results/step8_magma/gene_evidence_tiers.tsv", "magma/gene_evidence_tiers.tsv")
    qc = copy_table("results/step8_magma/magma_qc_report.tsv", "magma/magma_qc.tsv")
    if not genes.empty:
        sig = genes.loc[fnum(genes.get("gene_bonf_p")) <= 0.05].copy()
        write(sig, "magma/significant_genes.tsv")
        plot = genes.copy()
        plot["minus_log10_p"] = -np.log10(fnum(plot.get("p")).clip(lower=1e-300))
        write(plot[["pair_id", "trait1", "trait2", "gene_id", "gene_symbol", "chr", "start", "stop", "p", "gene_bonf_p", "minus_log10_p"]], "magma/magma_gene_manhattan_input.tsv")
        top = plot.sort_values(["pair_id", "gene_bonf_p", "p"], na_position="last").groupby("pair_id").head(50)
        write(top, "magma/magma_top_gene_input.tsv")
        write(top[["pair_id", "gene_symbol", "gene_bonf_p", "z"]], "magma/magma_gene_heatmap_input.tsv")
    else:
        for rel in ["significant_genes.tsv", "magma_gene_manhattan_input.tsv", "magma_top_gene_input.tsv", "magma_gene_heatmap_input.tsv"]:
            write_empty(f"magma/{rel}", [])
    simple_readme("magma", "MAGMA Gene Results", [("gene_results.tsv", "All MAGMA gene-level results from full MTAG."), ("significant_genes.tsv", "Bonferroni-significant genes.")])

    files = {
        "gsea_go_bp.tsv": "results/step8_pathway/pair_specific_gsea_go_bp.tsv",
        "gsea_kegg.tsv": "results/step8_pathway/pair_specific_gsea_kegg.tsv",
        "gsea_reactome.tsv": "results/step8_pathway/pair_specific_gsea_reactome.tsv",
        "ora_go_bp.tsv": "results/step8_pathway/pair_specific_ora_go_bp.tsv",
        "ora_kegg.tsv": "results/step8_pathway/pair_specific_ora_kegg.tsv",
        "ora_reactome.tsv": "results/step8_pathway/pair_specific_ora_reactome.tsv",
    }
    dfs = []
    for out, src in files.items():
        df = copy_table(src, f"pathway/{out}")
        if not df.empty:
            df["analysis_file"] = out
            dfs.append(df)
    allp = pd.concat(dfs, ignore_index=True) if dfs else pd.DataFrame()
    if not allp.empty:
        allp["significant_fdr05"] = fnum(allp.get("p.adjust")) <= 0.05
        summary = allp.groupby(["pair_id", "analysis_file"], as_index=False).agg(n_pathways=("pathway_id", "count"), n_sig_fdr05=("significant_fdr05", "sum"))
        write(summary, "pathway/pathway_summary.tsv")
        dot = allp.sort_values(["pair_id", "analysis_file", "p.adjust"], na_position="last").groupby(["pair_id", "analysis_file"]).head(20)
        write(dot, "pathway/pathway_dotplot_input.tsv")
        write(dot, "pathway/pathway_heatmap_input.tsv")
        write(dot.loc[dot["significant_fdr05"]], "pathway/pathway_network_input.tsv")
    else:
        for rel in ["pathway_summary.tsv", "pathway_dotplot_input.tsv", "pathway_heatmap_input.tsv", "pathway_network_input.tsv"]:
            write_empty(f"pathway/{rel}", [])
    simple_readme("pathway", "Pathway Results", [("gsea_go_bp.tsv/gsea_kegg.tsv/gsea_reactome.tsv", "GSEA pathway tables."), ("ora_go_bp.tsv/ora_kegg.tsv/ora_reactome.tsv", "ORA pathway tables."), ("pathway_summary.tsv", "Per-pair pathway counts.")])


def organize_integrated_visual_qc():
    pair = read("results/step8_mtag/evidence_pair_manifest.tsv")
    if pair.empty:
        pair = read("results/phase0_extension/step2_input_tables/priority_pairs.tsv")
    pair = filter_out_chronic(pair)
    base_cols = [c for c in ["pair_id", "trait1", "trait2"] if c in pair.columns]
    master = pair[base_cols].drop_duplicates().copy() if base_cols else pd.DataFrame(columns=["pair_id", "trait1", "trait2"])
    for c in ["ldsc_sig", "lava_sig"]:
        master[c] = pair.set_index("pair_id")[c].reindex(master["pair_id"]).fillna(False).values if c in pair.columns and not master.empty else False

    def merge_count(rel, colname, key="pair_id"):
        nonlocal master
        df = read(RESULTS / rel)
        if df.empty or key not in df.columns:
            master[colname] = 0
            return
        df = filter_out_chronic(df)
        count_col = "SNP" if "SNP" in df.columns else ("locus_id" if "locus_id" in df.columns else df.columns[0])
        s = df.groupby(key)[count_col].nunique().rename(colname)
        master = master.merge(s, on="pair_id", how="left")
        master[colname] = master[colname].fillna(0).astype(int)

    merge_count("placo/placo_significant_snps.tsv", "placo_snp_count")
    merge_count("cpassoc/cpassoc_significant_snps.tsv", "cpassoc_snp_count")
    merge_count("coloc/coloc_positive_loci.tsv", "coloc_loci_count")
    merge_count("susie_coloc/susie_positive_loci.tsv", "susie_loci_count")
    merge_count("mtag/significant_snps.tsv", "mtag_snp_count")
    merge_count("magma/significant_genes.tsv", "magma_sig_gene_count")
    for rel, col in [("pathway/gsea_go_bp.tsv", "go_sig_count"), ("pathway/gsea_kegg.tsv", "kegg_sig_count"), ("pathway/gsea_reactome.tsv", "reactome_sig_count")]:
        df = read(RESULTS / rel)
        if not df.empty and "pair_id" in df.columns:
            df = filter_out_chronic(df)
            sig = df.loc[fnum(df.get("p.adjust")) <= 0.05]
            s = sig.groupby("pair_id")["pathway_id"].nunique().rename(col)
            master = master.merge(s, on="pair_id", how="left")
            master[col] = master[col].fillna(0).astype(int)
        else:
            master[col] = 0

    count_cols = ["placo_snp_count", "cpassoc_snp_count", "coloc_loci_count", "susie_loci_count", "mtag_snp_count", "magma_sig_gene_count", "go_sig_count", "kegg_sig_count", "reactome_sig_count"]
    score = master[["ldsc_sig", "lava_sig"]].astype(str).isin(["True", "TRUE", "1", "true"]).sum(axis=1) + (master[count_cols] > 0).sum(axis=1)
    master["overall_evidence_level"] = np.where(score >= 5, "high", np.where(score >= 3, "moderate", "supporting"))
    write(master, "integrated/pair_master_summary.tsv")
    write(master.loc[master["overall_evidence_level"].eq("high")], "integrated/high_confidence_pairs.tsv")
    write(read("results/magma/significant_genes.tsv"), "integrated/high_confidence_genes.tsv")
    write(read("results/coloc/coloc_positive_loci.tsv"), "integrated/high_confidence_loci.tsv")
    path = read("results/pathway/pathway_network_input.tsv")
    write(path, "integrated/high_confidence_pathways.tsv")
    simple_readme("integrated", "Integrated Evidence Summary", [("pair_master_summary.tsv", "One row per disease pair integrating all evidence types."), ("high_confidence_pairs.tsv", "Pairs with high overall evidence level.")])

    mapping = {
        "01_ldsc_heatmap.tsv": "ldsc/ldsc_rg_heatmap_input.tsv",
        "02_ldsc_network.tsv": "ldsc/ldsc_rg_network_input.tsv",
        "03_lava_heatmap.tsv": "lava/lava_local_rg_heatmap_input.tsv",
        "04_lava_track.tsv": "lava/lava_local_rg_track_input.tsv",
        "05_placo_manhattan.tsv": "placo/placo_manhattan_input.tsv",
        "06_cpassoc_manhattan.tsv": "cpassoc/cpassoc_manhattan_input.tsv",
        "07_coloc_locus.tsv": "coloc/coloc_locus_plot_input.tsv",
        "08_susie_locus.tsv": "susie_coloc/susie_locus_plot_input.tsv",
        "09_mtag_manhattan.tsv": "mtag/mtag_manhattan_input.tsv",
        "10_magma_gene_manhattan.tsv": "magma/magma_gene_manhattan_input.tsv",
        "11_magma_gene_dotplot.tsv": "magma/magma_top_gene_input.tsv",
        "12_pathway_dotplot.tsv": "pathway/pathway_dotplot_input.tsv",
        "13_pathway_heatmap.tsv": "pathway/pathway_heatmap_input.tsv",
        "14_pair_network.tsv": "integrated/pair_master_summary.tsv",
        "15_gene_network.tsv": "integrated/high_confidence_genes.tsv",
        "16_pathway_network.tsv": "pathway/pathway_network_input.tsv",
    }
    for dest, src in mapping.items():
        df = read(RESULTS / src)
        write(df, f"visualization_inputs/{dest}")
    simple_readme("visualization_inputs", "Visualization Input Tables", [(k, f"Copied from `{v}`.") for k, v in mapping.items()])

    qc = pd.DataFrame([
        {"metric": "pair_total", "value": len(master)},
        {"metric": "LDSC_positive_pairs", "value": int(master["ldsc_sig"].astype(str).isin(["True", "TRUE", "1", "true"]).sum())},
        {"metric": "LAVA_positive_pairs", "value": int(master["lava_sig"].astype(str).isin(["True", "TRUE", "1", "true"]).sum())},
        {"metric": "PLACO_SNPs", "value": int(master["placo_snp_count"].sum())},
        {"metric": "CPASSOC_SNPs", "value": int(master["cpassoc_snp_count"].sum())},
        {"metric": "coloc_loci", "value": int(master["coloc_loci_count"].sum())},
        {"metric": "SuSiE_loci", "value": int(master["susie_loci_count"].sum())},
        {"metric": "MAGMA_significant_genes", "value": int(master["magma_sig_gene_count"].sum())},
        {"metric": "GO_significant_pathways", "value": int(master["go_sig_count"].sum())},
        {"metric": "KEGG_significant_pathways", "value": int(master["kegg_sig_count"].sum())},
        {"metric": "Reactome_significant_pathways", "value": int(master["reactome_sig_count"].sum())},
    ])
    write(qc, "qc/global_qc_report.tsv")
    simple_readme("qc", "Global QC", [("global_qc_report.tsv", "Project-wide counts after final exclusions.")])


def final_inventory():
    descriptions = {
        "README.md": "Directory-level description and usage notes.",
        "pair_master_summary.tsv": "Integrated pair-level evidence table.",
        "global_qc_report.tsv": "Global final QC metrics.",
    }
    rows = []
    for p in sorted(RESULTS.rglob("*")):
        if not p.is_file():
            continue
        try:
            rel = p.relative_to(RESULTS)
        except Exception:
            continue
        if rel.parts[0] not in OUT_DIRS and p.name != "final_project_inventory.tsv":
            continue
        ftype = p.suffix.lstrip(".") or "file"
        rows.append({
            "file_path": str(rel),
            "file_type": ftype,
            "rows": row_count(p) if ftype in {"tsv", "csv"} else "",
            "size_bytes": p.stat().st_size,
            "description": descriptions.get(p.name, "Paper-ready organized result or plotting input."),
        })
    inv = pd.DataFrame(rows)
    write(inv, "final_project_inventory.tsv")


def main():
    ensure_dirs()
    organize_ldsc()
    organize_lava()
    organize_mtag()
    organize_placo_cpassoc()
    organize_coloc_susie()
    organize_magma_pathway()
    organize_integrated_visual_qc()
    final_inventory()
    print("FINAL_RESULTS_STRUCTURE_DONE")


if __name__ == "__main__":
    main()
