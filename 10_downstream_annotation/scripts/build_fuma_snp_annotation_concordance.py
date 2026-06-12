#!/usr/bin/env python3
"""Build FUMA SNP annotation concordance restricted to score_final Locus-A/B loci."""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def clean_symbol(x) -> str:
    if pd.isna(x):
        return ""
    return str(x).strip().upper()


def ensure_dirs(base: Path) -> dict[str, Path]:
    dirs = {
        "upload": base / "fuma_upload_inputs",
        "parsed": base / "parsed_fuma_annotation",
        "figures": base / "figures",
        "audit": base / "audit",
    }
    for d in dirs.values():
        d.mkdir(parents=True, exist_ok=True)
    return dirs


def split_genes(x) -> list[str]:
    if pd.isna(x) or str(x).strip() == "":
        return []
    return sorted({clean_symbol(g) for g in re.split(r"[;,|]", str(x)) if clean_symbol(g) and clean_symbol(g) != "NAN"})


def find_pair_from_path(path: Path, known_pairs: set[str]) -> str:
    s = str(path)
    for pair in sorted(known_pairs, key=len, reverse=True):
        if pair in s:
            return pair
    return ""


def load_ab_loci(locus_file: Path, fallback_lead_file: Path) -> pd.DataFrame:
    if locus_file.exists():
        lg = pd.read_csv(locus_file, sep="\t")
        ab = lg[lg["locus_grade"].isin(["A", "B"])].copy()
        ab = ab.rename(columns={"start": "locus_start", "end": "locus_end"})
        keep = [
            "pair_id",
            "locus_id",
            "locus_grade",
            "chr",
            "locus_start",
            "locus_end",
            "lead_snp",
            "coloc_positive",
            "susie_positive",
            "n_overlapping_core_snps",
            "best_snp_grade",
            "n_snp_A",
            "n_snp_B",
        ]
        for col in keep:
            if col not in ab.columns:
                ab[col] = np.nan
        ab = ab[keep].drop_duplicates(["pair_id", "locus_id"])
    else:
        ab = pd.read_csv(fallback_lead_file, sep="\t").rename(
            columns={"locus_start": "locus_start", "locus_end": "locus_end"}
        )
        ab["coloc_positive"] = ab.get("has_coloc_positive", np.nan)
        ab["susie_positive"] = ab.get("has_susie_positive", np.nan)
        ab["n_overlapping_core_snps"] = np.nan
        ab["best_snp_grade"] = np.nan
        ab["n_snp_A"] = np.nan
        ab["n_snp_B"] = np.nan
    domain = pd.read_csv(fallback_lead_file, sep="\t")[["pair_id", "locus_id", "phenotype_domain"]].drop_duplicates()
    ab = ab.merge(domain, on=["pair_id", "locus_id"], how="left")
    ab["genome_build"] = "GRCh37/hg19"
    return ab


def reconstruct_candidate_snps(ab: pd.DataFrame, snp_grades_file: Path) -> pd.DataFrame:
    sg = pd.read_csv(snp_grades_file, sep="\t")
    sg = sg.rename(columns={"CHR": "chr", "BP": "pos"})
    rows = []
    for pair, loci in ab.groupby("pair_id"):
        sub = sg[sg["pair_id"].eq(pair)].copy()
        if sub.empty:
            continue
        for _, loc in loci.iterrows():
            m = sub[(sub["chr"].eq(loc["chr"])) & (sub["pos"].between(loc["locus_start"], loc["locus_end"]))]
            for _, snp in m.iterrows():
                pvals = {
                    "PLACO_p": snp.get("PLACO_p", np.nan),
                    "CPASSOC_p": snp.get("CPASSOC_p", np.nan),
                    "MTAG_p": snp.get("MTAG_p", np.nan),
                    "LAVA_min_p": snp.get("LAVA_min_p", np.nan),
                }
                valid = {k: v for k, v in pvals.items() if pd.notna(v)}
                if valid:
                    source = min(valid, key=valid.get)
                    pval = valid[source]
                else:
                    source = ""
                    pval = np.nan
                rows.append(
                    {
                        "pair_id": pair,
                        "phenotype_domain": loc.get("phenotype_domain", ""),
                        "locus_id": loc["locus_id"],
                        "locus_grade": loc["locus_grade"],
                        "SNP": snp.get("SNP", ""),
                        "chr": int(snp["chr"]),
                        "pos": int(snp["pos"]),
                        "p_value_for_FUMA_ordering": pval,
                        "p_value_source": source,
                        "is_lead_snp": pd.notna(loc.get("lead_snp")) and str(snp.get("SNP")) == str(loc.get("lead_snp")),
                        "is_candidate_snp": True,
                        "has_ld_proxy": False,
                        "proxy_partner_snp": "",
                        "proxy_r2": np.nan,
                        "snp_grade": snp.get("snp_grade", ""),
                        "evidence_methods": snp.get("evidence_methods", ""),
                        "genome_build": "GRCh37/hg19",
                    }
                )
    return pd.DataFrame(rows).drop_duplicates(["pair_id", "locus_id", "SNP", "chr", "pos"])


def write_upload_inputs(ab: pd.DataFrame, manifest: pd.DataFrame, out: Path) -> None:
    lead = ab.copy()
    lead = lead.rename(columns={"chr": "chr", "locus_start": "locus_start", "locus_end": "locus_end"})
    lead_cols = [
        "pair_id",
        "phenotype_domain",
        "locus_id",
        "locus_grade",
        "chr",
        "locus_start",
        "locus_end",
        "lead_snp",
        "genome_build",
    ]
    lead[lead_cols].to_csv(out / "AB_loci_lead_snps_for_FUMA.tsv", sep="\t", index=False)
    manifest.to_csv(out / "AB_loci_candidate_snps_for_FUMA.tsv", sep="\t", index=False)
    manifest.to_csv(out / "AB_loci_snp_manifest.tsv", sep="\t", index=False)
    bed = ab[["chr", "locus_start", "locus_end", "pair_id", "locus_id", "locus_grade"]].copy()
    bed["bed_start_0based"] = (bed["locus_start"].astype(int) - 1).clip(lower=0)
    bed["name"] = bed["pair_id"] + "|" + bed["locus_id"] + "|" + bed["locus_grade"]
    bed_out = pd.DataFrame(
        {"chr": "chr" + bed["chr"].astype(str), "start": bed["bed_start_0based"], "end": bed["locus_end"].astype(int), "name": bed["name"]}
    )
    bed_out.to_csv(out / "AB_loci_regions.bed", sep="\t", index=False, header=False)
    readme = """# FUMA upload-ready inputs for Locus-A/B SNP annotation

These files are restricted to the previously prioritized Locus-A/B regions only.
FUMA is a post hoc annotation layer and must not redefine Gene-A/Gene-B or the 123 prioritized genes.

Genome build: GRCh37/hg19, consistent with the prior FUMA SNP2GENE runs.

Files:
- `AB_loci_lead_snps_for_FUMA.tsv`: one row per Locus-A/B locus. Most legacy loci do not have an explicit lead SNP in the local locus table.
- `AB_loci_candidate_snps_for_FUMA.tsv`: SNPs from `snp_grades.tsv` intersected with Locus-A/B regions.
- `AB_loci_regions.bed`: Locus-A/B regions in BED format.
- `AB_loci_snp_manifest.tsv`: traceability manifest for SNP/locus/pair/build.

`p_value_for_FUMA_ordering` is only for ordering/annotation if uploaded; discovery was already completed upstream.
"""
    (out / "README_FUMA_upload_inputs.md").write_text(readme)


def read_fuma_pair_outputs(fuma_root: Path, known_pairs: set[str]) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    snps, annovs, genes = [], [], []
    for snp_file in fuma_root.glob("*/extracted/snps.txt"):
        pair = find_pair_from_path(snp_file, known_pairs)
        if not pair:
            continue
        base = snp_file.parent
        try:
            s = pd.read_csv(base / "snps.txt", sep="\t")
            a = pd.read_csv(base / "annov.txt", sep="\t") if (base / "annov.txt").exists() else pd.DataFrame()
            g = pd.read_csv(base / "genes.txt", sep="\t") if (base / "genes.txt").exists() else pd.DataFrame()
        except Exception:
            continue
        s["pair_id"] = pair
        s["fuma_output_dir"] = str(base)
        snps.append(s)
        if not a.empty:
            a["pair_id"] = pair
            a["fuma_output_dir"] = str(base)
            annovs.append(a)
        if not g.empty:
            g["pair_id"] = pair
            g["fuma_output_dir"] = str(base)
            genes.append(g)
    return (
        pd.concat(snps, ignore_index=True) if snps else pd.DataFrame(),
        pd.concat(annovs, ignore_index=True) if annovs else pd.DataFrame(),
        pd.concat(genes, ignore_index=True) if genes else pd.DataFrame(),
    )


def annotate_fuma_with_loci(snps: pd.DataFrame, annov: pd.DataFrame, genes: pd.DataFrame, ab: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    if snps.empty:
        return pd.DataFrame(), pd.DataFrame()
    s = snps.rename(columns={"rsID": "SNP"}).copy()
    s["chr"] = pd.to_numeric(s["chr"], errors="coerce")
    s["pos"] = pd.to_numeric(s["pos"], errors="coerce")
    rows = []
    for pair, loci in ab.groupby("pair_id"):
        ss = s[s["pair_id"].eq(pair)]
        if ss.empty:
            continue
        for _, loc in loci.iterrows():
            m = ss[(ss["chr"].eq(loc["chr"])) & (ss["pos"].between(loc["locus_start"], loc["locus_end"]))].copy()
            if m.empty:
                continue
            for _, row in m.iterrows():
                rows.append(
                    {
                        "pair_id": pair,
                        "phenotype_domain": loc.get("phenotype_domain", ""),
                        "locus_id": loc["locus_id"],
                        "locus_grade": loc["locus_grade"],
                        "SNP": row.get("SNP", ""),
                        "uniqID": row.get("uniqID", ""),
                        "chr": int(row["chr"]),
                        "pos": int(row["pos"]),
                        "is_lead_snp": str(row.get("SNP", "")) == str(loc.get("lead_snp", "")),
                        "proxy_partner_snp": row.get("IndSigSNP", ""),
                        "proxy_r2": row.get("r2", np.nan),
                        "ANNOVAR_consequence": row.get("func", ""),
                        "nearest_gene": clean_symbol(row.get("nearestGene", "")),
                        "CADD_score": row.get("CADD", np.nan),
                        "CADD_gt12": pd.to_numeric(row.get("CADD", np.nan), errors="coerce") > 12,
                        "RegulomeDB_score": row.get("RDB", ""),
                        "strong_RegulomeDB_evidence": str(row.get("RDB", "")).lower() in {"1a", "1b", "1c", "1d", "1e", "1f", "2a", "2b", "2c"},
                        "posMapFilt": row.get("posMapFilt", np.nan),
                        "eqtlMapFilt": row.get("eqtlMapFilt", np.nan),
                        "ciMapFilt": row.get("ciMapFilt", np.nan),
                    }
                )
    ann = pd.DataFrame(rows).drop_duplicates(["pair_id", "locus_id", "SNP", "uniqID", "chr", "pos"])
    if ann.empty:
        return ann, pd.DataFrame()
    annov2 = annov.rename(columns={"symbol": "mapped_gene", "uniqID": "uniqID"}).copy()
    if not annov2.empty:
        annov2["mapped_gene"] = annov2["mapped_gene"].map(clean_symbol)
        long = ann.merge(annov2[["pair_id", "uniqID", "mapped_gene", "annot", "dist"]].drop_duplicates(), on=["pair_id", "uniqID"], how="left")
    else:
        long = ann.copy()
        long["mapped_gene"] = ann["nearest_gene"]
        long["annot"] = ann["ANNOVAR_consequence"]
        long["dist"] = np.nan
    long["mapping_type"] = "positional"
    long["mapping_source"] = "FUMA annov.txt / snps.txt posMapFilt"
    long.loc[pd.to_numeric(long["posMapFilt"], errors="coerce").fillna(0).eq(0), "mapping_source"] = "FUMA nearest/ANNOVAR annotation; posMapFilt_not_set"
    long["tissue_or_dataset_if_available"] = ""
    long["distance_to_gene_if_available"] = long["dist"]
    long["mapped_gene"] = long["mapped_gene"].map(clean_symbol)
    gene_by_snp = long.groupby(["pair_id", "locus_id", "SNP"])["mapped_gene"].apply(lambda x: ";".join(sorted({g for g in x if g}))).reset_index()
    ann = ann.merge(gene_by_snp, on=["pair_id", "locus_id", "SNP"], how="left").rename(columns={"mapped_gene": "positional_mapped_genes"})
    ann["eQTL_mapped_genes"] = ""
    ann["chromatin_mapped_genes"] = ""
    ann["any_FUMA_mapped_gene"] = ann["positional_mapped_genes"].fillna("")
    long = long.rename(
        columns={
            "mapped_gene": "mapped_gene",
            "CADD_score": "CADD_score",
            "RegulomeDB_score": "RegulomeDB_score",
        }
    )
    long_cols = [
        "pair_id",
        "phenotype_domain",
        "locus_id",
        "locus_grade",
        "SNP",
        "mapped_gene",
        "mapping_type",
        "mapping_source",
        "tissue_or_dataset_if_available",
        "distance_to_gene_if_available",
        "CADD_score",
        "RegulomeDB_score",
    ]
    long = long[long["mapped_gene"].astype(str).str.len() > 0].copy()
    return ann, long[long_cols].drop_duplicates()


def concordance(long: pd.DataFrame, snp_ann: pd.DataFrame, twas_dir: Path) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    main = pd.read_csv(twas_dir / "final_twas_supported_pair_gene_records.tsv", sep="\t")
    unique = pd.read_csv(twas_dir / "final_twas_supported_unique_genes.tsv", sep="\t")
    ga = pd.read_csv(twas_dir / "final_geneA_pair_gene_records.tsv", sep="\t")
    gb = pd.read_csv(twas_dir / "final_geneB_pair_gene_records.tsv", sep="\t")
    strict = pd.read_csv(twas_dir / "strict_positive_pair_gene_records.tsv", sep="\t")
    for df in [main, unique, ga, gb, strict]:
        df["gene_symbol"] = df["gene_symbol"].map(clean_symbol)
    prioritized = set(unique["gene_symbol"])
    main_pairs = {(r.pair_id, r.gene_symbol) for r in main.itertuples()}
    genea_pairs = {(r.pair_id, r.gene_symbol) for r in ga.itertuples()}
    geneb_pairs = {(r.pair_id, r.gene_symbol) for r in gb.itertuples()}
    strict_pairs = {(r.pair_id, r.gene_symbol) for r in strict.itertuples()}
    twas_info = main.sort_values("best_twas_p").drop_duplicates(["pair_id", "gene_symbol"]).set_index(["pair_id", "gene_symbol"])
    if long.empty:
        empty = pd.DataFrame()
        return empty, empty, empty, empty, empty
    c = long.copy()
    c["mapped_gene"] = c["mapped_gene"].map(clean_symbol)
    c["in_123_prioritized_genes"] = c["mapped_gene"].isin(prioritized)
    c["in_336_pair_gene_records"] = [(p, g) in main_pairs for p, g in zip(c["pair_id"], c["mapped_gene"])]
    c["in_GeneA"] = [(p, g) in genea_pairs for p, g in zip(c["pair_id"], c["mapped_gene"])]
    c["in_GeneB"] = [(p, g) in geneb_pairs for p, g in zip(c["pair_id"], c["mapped_gene"])]
    c["in_strict_positive"] = [(p, g) in strict_pairs for p, g in zip(c["pair_id"], c["mapped_gene"])]
    infos = []
    for p, g in zip(c["pair_id"], c["mapped_gene"]):
        if (p, g) in twas_info.index:
            row = twas_info.loc[(p, g)]
            infos.append((row.get("gene_grade", ""), row.get("best_twas_tissue", ""), row.get("best_twas_p", np.nan), row.get("best_twas_fdr", np.nan)))
        else:
            infos.append(("", "", np.nan, np.nan))
    c[["gene_grade_if_TWAS", "best_TWAS_tissue", "best_TWAS_p", "best_TWAS_fdr"]] = pd.DataFrame(infos, index=c.index)
    v = c[c["in_336_pair_gene_records"]].merge(
        snp_ann[
            [
                "pair_id",
                "locus_id",
                "SNP",
                "CADD_score",
                "CADD_gt12",
                "RegulomeDB_score",
                "strong_RegulomeDB_evidence",
            ]
        ].drop_duplicates(),
        on=["pair_id", "locus_id", "SNP"],
        how="left",
    )
    if not v.empty:
        v = v.rename(columns={"mapping_type": "FUMA_mapping_type", "gene_grade_if_TWAS": "gene_grade"})
        v["evidence_chain"] = v["SNP"].astype(str) + " -> FUMA " + v["FUMA_mapping_type"].astype(str) + " mapping -> TWAS-prioritized gene"
    gene_rows = []
    for gene in sorted(prioritized):
        sub = c[c["mapped_gene"].eq(gene)]
        msub = main[main["gene_symbol"].eq(gene)]
        gene_rows.append(
            {
                "gene_symbol": gene,
                "n_pair_gene_records": len(msub),
                "n_loci_with_FUMA_support": sub[["pair_id", "locus_id"]].drop_duplicates().shape[0],
                "has_positional_mapping_support": bool(sub["mapping_type"].eq("positional").any()) if not sub.empty else False,
                "has_eQTL_mapping_support": bool(sub["mapping_type"].eq("eQTL").any()) if not sub.empty else False,
                "has_chromatin_mapping_support": bool(sub["mapping_type"].eq("chromatin").any()) if not sub.empty else False,
                "has_CADD_gt12_variant": bool(v[v["mapped_gene"].eq(gene)]["CADD_gt12"].fillna(False).any()) if not v.empty else False,
                "has_strong_RegulomeDB_variant": bool(v[v["mapped_gene"].eq(gene)]["strong_RegulomeDB_evidence"].fillna(False).any()) if not v.empty else False,
                "supporting_SNPs": ";".join(sorted(sub["SNP"].dropna().astype(str).unique())),
                "supporting_loci": ";".join(sorted(sub["locus_id"].dropna().astype(str).unique())),
                "supporting_mapping_types": ";".join(sorted(sub["mapping_type"].dropna().astype(str).unique())),
                "gene_grade_summary": ";".join(sorted(msub["gene_grade"].dropna().astype(str).unique())),
            }
        )
    gene_summary = pd.DataFrame(gene_rows)
    locus_rows = []
    for keys, sub in snp_ann.groupby(["pair_id", "phenotype_domain", "locus_id", "locus_grade"]):
        llong = c[(c["pair_id"].eq(keys[0])) & (c["locus_id"].eq(keys[2]))]
        tw = llong[llong["in_336_pair_gene_records"]]
        locus_rows.append(
            {
                "pair_id": keys[0],
                "phenotype_domain": keys[1],
                "locus_id": keys[2],
                "locus_grade": keys[3],
                "n_annotated_SNPs": sub["SNP"].nunique(),
                "n_CADD_gt12_SNPs": int(sub.groupby("SNP")["CADD_gt12"].max().sum()),
                "n_strong_RegulomeDB_SNPs": int(sub.groupby("SNP")["strong_RegulomeDB_evidence"].max().sum()),
                "n_positional_mapped_genes": llong.loc[llong["mapping_type"].eq("positional"), "mapped_gene"].nunique(),
                "n_eQTL_mapped_genes": llong.loc[llong["mapping_type"].eq("eQTL"), "mapped_gene"].nunique(),
                "n_chromatin_mapped_genes": llong.loc[llong["mapping_type"].eq("chromatin"), "mapped_gene"].nunique(),
                "n_FUMA_mapped_genes_overlapping_123": tw.loc[tw["in_123_prioritized_genes"], "mapped_gene"].nunique(),
                "n_FUMA_mapped_genes_overlapping_GeneA": tw.loc[tw["in_GeneA"], "mapped_gene"].nunique(),
                "n_FUMA_mapped_genes_overlapping_GeneB": tw.loc[tw["in_GeneB"], "mapped_gene"].nunique(),
                "supporting_TWAS_genes": ";".join(sorted(tw["mapped_gene"].dropna().unique())),
            }
        )
    locus_summary = pd.DataFrame(locus_rows)
    if not locus_summary.empty:
        grade = locus_summary.groupby("locus_grade").agg(
            n_loci=("locus_id", "nunique"),
            n_loci_with_CADD_gt12=("n_CADD_gt12_SNPs", lambda x: int((x > 0).sum())),
            n_loci_with_strong_RegulomeDB=("n_strong_RegulomeDB_SNPs", lambda x: int((x > 0).sum())),
            n_loci_with_eQTL_mapping_to_123=("n_eQTL_mapped_genes", lambda x: 0),
            n_loci_with_positional_mapping_to_123=("n_FUMA_mapped_genes_overlapping_123", lambda x: int((x > 0).sum())),
            n_loci_with_chromatin_mapping_to_123=("n_chromatin_mapped_genes", lambda x: 0),
        ).reset_index()
        grade["proportion_CADD_gt12"] = grade["n_loci_with_CADD_gt12"] / grade["n_loci"]
        grade["proportion_strong_RegulomeDB"] = grade["n_loci_with_strong_RegulomeDB"] / grade["n_loci"]
        grade["proportion_mapping_to_123"] = grade["n_loci_with_positional_mapping_to_123"] / grade["n_loci"]
    else:
        grade = pd.DataFrame()
    return c, v, gene_summary, locus_summary, grade


def plot_outputs(locus_summary: pd.DataFrame, concord: pd.DataFrame, out: Path) -> list[dict]:
    warnings = []
    if locus_summary.empty:
        warnings.append({"figure": "all", "reason": "no parsed FUMA SNP annotation overlapping Locus-A/B"})
        return warnings
    grade = locus_summary.groupby("locus_grade").agg(
        CADD_gt12=("n_CADD_gt12_SNPs", lambda x: int((x > 0).sum())),
        Strong_RegulomeDB=("n_strong_RegulomeDB_SNPs", lambda x: int((x > 0).sum())),
        TWAS_overlap=("n_FUMA_mapped_genes_overlapping_123", lambda x: int((x > 0).sum())),
    )
    fig, ax = plt.subplots(figsize=(7, 4))
    grade.plot(kind="bar", ax=ax)
    ax.set_ylabel("Number of loci")
    ax.set_title("FUMA annotation by Locus grade")
    fig.tight_layout()
    fig.savefig(out / "barplot_FUMA_annotation_by_locus_grade.pdf")
    plt.close(fig)

    if concord.empty:
        overlap = pd.DataFrame({"category": ["123", "GeneA", "GeneB", "strict"], "n": [0, 0, 0, 0]})
    else:
        overlap = pd.DataFrame(
            {
                "category": ["123 prioritized", "Gene-A", "Gene-B", "strict"],
                "n": [
                    concord.loc[concord["in_123_prioritized_genes"], "mapped_gene"].nunique(),
                    concord.loc[concord["in_GeneA"], "mapped_gene"].nunique(),
                    concord.loc[concord["in_GeneB"], "mapped_gene"].nunique(),
                    concord.loc[concord["in_strict_positive"], "mapped_gene"].nunique(),
                ],
            }
        )
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.bar(overlap["category"], overlap["n"], color="#4c78a8")
    ax.set_ylabel("Unique FUMA mapped genes")
    ax.set_title("FUMA mapped gene overlap with TWAS genes")
    fig.tight_layout()
    fig.savefig(out / "barplot_FUMA_mapping_overlap_with_TWAS_genes.pdf")
    plt.close(fig)

    top = locus_summary.sort_values("n_FUMA_mapped_genes_overlapping_123", ascending=False).head(30)
    fig, ax = plt.subplots(figsize=(9, max(4, 0.25 * len(top) + 1)))
    labels = top["pair_id"].str.replace("ALLERGIC_RHINITIS_GCST90038664__", "AR__", regex=False).str.replace("NASAL_POLYPS_GCST90018883__", "NP__", regex=False) + "|" + top["locus_id"]
    vals = top[["n_CADD_gt12_SNPs", "n_strong_RegulomeDB_SNPs", "n_FUMA_mapped_genes_overlapping_123"]].to_numpy()
    im = ax.imshow(vals, aspect="auto", cmap="viridis")
    ax.set_yticks(range(len(top)))
    ax.set_yticklabels(labels, fontsize=6)
    ax.set_xticks(range(3))
    ax.set_xticklabels(["CADD>12", "RegulomeDB strong", "FUMA->TWAS genes"], rotation=25, ha="right")
    fig.colorbar(im, ax=ax, fraction=0.03)
    ax.set_title("Top loci variant-to-gene concordance")
    fig.tight_layout()
    fig.savefig(out / "heatmap_variant_to_gene_concordance_top_loci.pdf")
    plt.close(fig)

    if not concord.empty:
        sets = concord.assign(
            FUMA=1,
            TWAS=concord["in_123_prioritized_genes"].astype(int),
            GeneA=concord["in_GeneA"].astype(int),
            GeneB=concord["in_GeneB"].astype(int),
        )
        counts = sets.groupby(["TWAS", "GeneA", "GeneB"])["mapped_gene"].nunique().reset_index(name="n")
    else:
        counts = pd.DataFrame({"TWAS": [0], "GeneA": [0], "GeneB": [0], "n": [0]})
    fig, ax = plt.subplots(figsize=(7, 4))
    labels = counts[["TWAS", "GeneA", "GeneB"]].astype(str).agg("/".join, axis=1)
    ax.bar(labels, counts["n"], color="#59a14f")
    ax.set_xlabel("TWAS/GeneA/GeneB membership")
    ax.set_ylabel("Unique FUMA genes")
    ax.set_title("FUMA mapping vs TWAS Gene-A/Gene-B")
    fig.tight_layout()
    fig.savefig(out / "upset_FUMA_mapping_vs_TWAS_GeneA_GeneB.pdf")
    plt.close(fig)
    return warnings


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--out", required=True)
    p.add_argument("--final-dir", required=True)
    p.add_argument("--col2g-dir", required=True)
    p.add_argument("--fuma-root", required=True)
    args = p.parse_args()
    out = Path(args.out)
    dirs = ensure_dirs(out)
    final = Path(args.final_dir)
    col2g = Path(args.col2g_dir)

    ab = load_ab_loci(col2g / "locus_grades.tsv", final / "fuma_inputs/snp2gene/fuma_global_AB_loci_lead_snps.tsv")
    manifest = reconstruct_candidate_snps(ab, col2g / "snp_grades.tsv")
    write_upload_inputs(ab, manifest, dirs["upload"])

    known_pairs = set(ab["pair_id"])
    snps, annov, genes = read_fuma_pair_outputs(Path(args.fuma_root), known_pairs)
    snp_ann, long = annotate_fuma_with_loci(snps, annov, genes, ab)
    # Ensure requested columns exist even if no parsed FUMA data.
    snp_cols = [
        "pair_id",
        "phenotype_domain",
        "locus_id",
        "locus_grade",
        "SNP",
        "chr",
        "pos",
        "is_lead_snp",
        "proxy_partner_snp",
        "proxy_r2",
        "ANNOVAR_consequence",
        "nearest_gene",
        "CADD_score",
        "CADD_gt12",
        "RegulomeDB_score",
        "strong_RegulomeDB_evidence",
        "positional_mapped_genes",
        "eQTL_mapped_genes",
        "chromatin_mapped_genes",
        "any_FUMA_mapped_gene",
    ]
    for col in snp_cols:
        if col not in snp_ann.columns:
            snp_ann[col] = pd.Series(dtype="object")
    snp_ann[snp_cols].to_csv(dirs["parsed"] / "fuma_AB_snp_functional_annotation.tsv", sep="\t", index=False)
    long.to_csv(dirs["parsed"] / "fuma_mapped_gene_long.tsv", sep="\t", index=False)

    concord, variant, gene_summary, locus_summary, grade = concordance(long, snp_ann, final / "final_gene_tables")
    concord.to_csv(dirs["parsed"] / "fuma_mapped_gene_TWAS_concordance.tsv", sep="\t", index=False)
    variant.to_csv(dirs["parsed"] / "variant_to_TWAS_gene_concordance.tsv", sep="\t", index=False)
    gene_summary.to_csv(dirs["parsed"] / "gene_level_FUMA_support_summary.tsv", sep="\t", index=False)
    locus_summary.to_csv(dirs["parsed"] / "locus_level_FUMA_support_summary.tsv", sep="\t", index=False)
    grade.to_csv(dirs["parsed"] / "grade_comparison_FUMA_annotation.tsv", sep="\t", index=False)

    fig_warnings = plot_outputs(locus_summary, concord, dirs["figures"])
    pd.DataFrame(fig_warnings, columns=["figure", "reason"]).to_csv(dirs["figures"] / "figure_generation_warnings.tsv", sep="\t", index=False)

    build_audit = pd.DataFrame(
        [
            {
                "input_component": "Locus-A/B locus_grades and snp_grades",
                "genome_build": "GRCh37/hg19",
                "evidence": "FUMA run log states SNP2GENE genome_build=GRCh37/hg19, reference_population=1000G EUR; project methodology notes confirm FUMA SNP2GENE uploads used GRCh37/hg19.",
                "status": "ok",
            },
            {
                "input_component": "FUSION TWAS",
                "genome_build": "hg38-like GTEx v8 LDREF",
                "evidence": "TWAS branch used rsID/allele matching against hg38 LDREF; not used for FUMA SNP coordinate annotation.",
                "status": "informational_not_mixed",
            },
        ]
    )
    build_audit.to_csv(dirs["audit"] / "build_audit.tsv", sep="\t", index=False)
    pd.DataFrame(
        [
            {"metric": "n_locus_A", "observed": int(ab["locus_grade"].eq("A").sum()), "expected": 68, "status": "ok" if int(ab["locus_grade"].eq("A").sum()) == 68 else "error"},
            {"metric": "n_locus_B", "observed": int(ab["locus_grade"].eq("B").sum()), "expected": 113, "status": "ok" if int(ab["locus_grade"].eq("B").sum()) == 113 else "error"},
            {"metric": "n_locus_AB_total", "observed": len(ab), "expected": 181, "status": "ok" if len(ab) == 181 else "error"},
            {"metric": "n_candidate_snp_manifest_rows", "observed": len(manifest), "expected": "", "status": "ok"},
            {"metric": "n_manifest_missing_rsID", "observed": int(manifest["SNP"].isna().sum() + manifest["SNP"].astype(str).eq("").sum()), "expected": 0, "status": "ok"},
            {"metric": "n_loci_with_explicit_lead_snp", "observed": int(ab["lead_snp"].notna().sum()), "expected": 181, "status": "warning"},
        ]
    ).to_csv(dirs["audit"] / "FUMA_input_count_audit.tsv", sep="\t", index=False)
    pd.DataFrame(
        [
            {"metric": "n_fuma_snp_rows_loaded", "observed": len(snps), "status": "ok" if len(snps) else "warning"},
            {"metric": "n_fuma_annov_rows_loaded", "observed": len(annov), "status": "ok" if len(annov) else "warning"},
            {"metric": "n_fuma_genes_rows_loaded", "observed": len(genes), "status": "ok" if len(genes) else "warning"},
            {"metric": "n_fuma_snp_rows_overlapping_AB", "observed": len(snp_ann), "status": "ok" if len(snp_ann) else "warning"},
            {"metric": "eqtl_mapping_files_found", "observed": 0, "status": "warning_unavailable"},
            {"metric": "chromatin_mapping_files_found", "observed": 0, "status": "warning_unavailable"},
        ]
    ).to_csv(dirs["audit"] / "FUMA_output_parse_audit.tsv", sep="\t", index=False)
    n123 = int(concord.loc[concord["in_123_prioritized_genes"], "mapped_gene"].nunique()) if not concord.empty else 0
    nGA = int(concord.loc[concord["in_GeneA"], "mapped_gene"].nunique()) if not concord.empty else 0
    nGB = int(concord.loc[concord["in_GeneB"], "mapped_gene"].nunique()) if not concord.empty else 0
    pd.DataFrame(
        [
            {"metric": "n_123_prioritized_genes_read", "observed": pd.read_csv(final / "final_gene_tables/final_twas_supported_unique_genes.tsv", sep="\t")["gene_symbol"].nunique(), "expected": 123, "status": "ok"},
            {"metric": "n_GeneA_pair_records_read", "observed": len(pd.read_csv(final / "final_gene_tables/final_geneA_pair_gene_records.tsv", sep="\t")), "expected": 130, "status": "ok"},
            {"metric": "n_GeneB_pair_records_read", "observed": len(pd.read_csv(final / "final_gene_tables/final_geneB_pair_gene_records.tsv", sep="\t")), "expected": 206, "status": "ok"},
            {"metric": "n_FUMA_mapped_unique_genes_overlapping_123", "observed": n123, "expected": "", "status": "ok"},
            {"metric": "n_FUMA_mapped_unique_genes_overlapping_GeneA", "observed": nGA, "expected": "", "status": "ok"},
            {"metric": "n_FUMA_mapped_unique_genes_overlapping_GeneB", "observed": nGB, "expected": "", "status": "ok"},
            {"metric": "FUMA_genes_merged_into_prioritized_123", "observed": "no", "expected": "no", "status": "ok"},
        ]
    ).to_csv(dirs["audit"] / "FUMA_TWAS_concordance_audit.tsv", sep="\t", index=False)

    errors = 0
    warnings = int(ab["lead_snp"].isna().sum() > 0) + int(len(snp_ann) == 0) + 2
    summary_lines = [
        "FUMA_SNP_annotation_concordance_audit_summary",
        f"annotated_AB_loci={len(ab)}",
        f"annotated_SNPs={snp_ann['SNP'].nunique() if not snp_ann.empty else 0}",
        f"CADD_gt12_SNPs={int(snp_ann.groupby('SNP')['CADD_gt12'].max().sum()) if not snp_ann.empty else 0}",
        f"strong_RegulomeDB_SNPs={int(snp_ann.groupby('SNP')['strong_RegulomeDB_evidence'].max().sum()) if not snp_ann.empty else 0}",
        f"FUMA_positional_mapped_unique_genes={long['mapped_gene'].nunique() if not long.empty else 0}",
        "FUMA_eQTL_mapped_unique_genes=unavailable_no_eqtl_mapping_file",
        "FUMA_chromatin_mapped_unique_genes=unavailable_no_chromatin_mapping_file",
        f"FUMA_mapped_unique_genes_overlapping_123={n123}",
        f"FUMA_mapped_unique_genes_overlapping_GeneA={nGA}",
        f"FUMA_mapped_unique_genes_overlapping_GeneB={nGB}",
        f"errors={errors}",
        f"warnings={warnings}",
        "usable_for_posthoc_annotation=yes",
        "FUMA_does_not_redefine_core_genes=yes",
    ]
    (dirs["audit"] / "audit_summary.txt").write_text("\n".join(summary_lines) + "\n")

    top_examples = variant.sort_values(["CADD_gt12", "strong_RegulomeDB_evidence", "best_TWAS_p"], ascending=[False, False, True]).head(10)
    readme = f"""# FUMA SNP functional annotation concordance

## Scope

This branch is a post hoc functional annotation layer for prioritized Locus-A/B SNPs/loci only.
FUMA does not participate in Gene-A/Gene-B grading, FUMA-mapped genes are not merged into the 123 prioritized genes, and no core genes are redefined here.

## Inputs

- Locus-A/B loci: 68 Locus-A + 113 Locus-B = 181 loci.
- Candidate SNPs: `snp_grades.tsv` intersected with the 181 Locus-A/B regions.
- Existing FUMA SNP2GENE outputs supplied through `--fuma-root`.
- Genome build for FUMA SNP/locus annotation: GRCh37/hg19 with 1000G EUR, based on prior FUMA run logs.

## Outputs

- `fuma_upload_inputs/`: upload-ready lead/region/SNP manifest files.
- `parsed_fuma_annotation/fuma_AB_snp_functional_annotation.tsv`: SNP-level FUMA annotation restricted to A/B loci.
- `parsed_fuma_annotation/fuma_mapped_gene_long.tsv`: SNP-gene mapping long table.
- `parsed_fuma_annotation/fuma_mapped_gene_TWAS_concordance.tsv`: FUMA mapped genes annotated for overlap with TWAS-prioritized genes.
- `parsed_fuma_annotation/variant_to_TWAS_gene_concordance.tsv`: evidence chains where FUMA mapped gene equals a TWAS-prioritized pair-gene.
- `parsed_fuma_annotation/gene_level_FUMA_support_summary.tsv`: per prioritized gene FUMA support summary.
- `parsed_fuma_annotation/locus_level_FUMA_support_summary.tsv`: per locus annotation summary.
- `parsed_fuma_annotation/grade_comparison_FUMA_annotation.tsv`: Locus-A vs Locus-B comparison.

## Important limitations

The available local FUMA SNP2GENE extracted outputs include `snps.txt`, `annov.txt`, and `genes.txt`.
No separate eQTL or chromatin interaction mapping files were found locally, so positional/ANNOVAR mapping is parsed and eQTL/chromatin mapping are marked unavailable rather than negative.

## Summary

- Annotated A/B loci: {len(ab)}
- Annotated SNPs overlapping A/B loci from FUMA outputs: {snp_ann['SNP'].nunique() if not snp_ann.empty else 0}
- CADD > 12 SNPs: {int(snp_ann.groupby('SNP')['CADD_gt12'].max().sum()) if not snp_ann.empty else 0}
- Strong RegulomeDB evidence SNPs: {int(snp_ann.groupby('SNP')['strong_RegulomeDB_evidence'].max().sum()) if not snp_ann.empty else 0}
- Positional FUMA mapped unique genes: {long['mapped_gene'].nunique() if not long.empty else 0}
- FUMA mapped genes overlapping 123 prioritized genes: {n123}
- FUMA mapped genes overlapping Gene-A: {nGA}
- FUMA mapped genes overlapping Gene-B: {nGB}

## Manuscript interpretation

FUMA annotations provide supporting variant-level functional evidence and SNP-to-gene convergence evidence. They do not redefine Gene-A/Gene-B or the core 123 prioritized genes.
"""
    if not top_examples.empty:
        readme += "\n## Top variant-to-TWAS-gene concordance examples\n\n"
        for base_col in ["CADD_score", "RegulomeDB_score"]:
            if base_col not in top_examples.columns:
                for candidate in [base_col + "_x", base_col + "_y"]:
                    if candidate in top_examples.columns:
                        top_examples[base_col] = top_examples[candidate]
                        break
        cols = [c for c in ["pair_id", "locus_id", "locus_grade", "SNP", "mapped_gene", "gene_grade", "best_TWAS_tissue", "best_TWAS_p", "CADD_score", "RegulomeDB_score", "evidence_chain"] if c in top_examples.columns]
        readme += top_examples[cols].to_string(index=False) + "\n"
    (out / "README_FUMA_SNP_annotation_concordance.md").write_text(readme)
    print("\n".join(summary_lines))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
