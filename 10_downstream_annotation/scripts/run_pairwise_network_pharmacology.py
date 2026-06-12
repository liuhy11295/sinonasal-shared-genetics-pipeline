#!/usr/bin/env python3
"""Pairwise network pharmacology for score_final TWAS-prioritized genes."""

from __future__ import annotations

import argparse
import json
import math
import re
import sys
import time
from collections import Counter, defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
import requests


DGIDB_URL = "https://dgidb.org/api/graphql"
STRING_NETWORK_URL = "https://string-db.org/api/tsv/network"
STRING_ENRICHMENT_URL = "https://string-db.org/api/tsv/enrichment"


DGIDB_QUERY = """
query($genes:[String!], $first:Int) {
  interactions(geneNames:$genes, first:$first) {
    nodes {
      id
      interactionScore
      evidenceScore
      gene { name conceptId }
      drug { name conceptId approved immunotherapy antiNeoplastic }
    }
  }
}
"""


def clean_gene(x: str) -> str:
    return str(x).strip().upper()


def safe_name(x: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", str(x))


def write_gene_list(path: Path, genes: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(genes) + ("\n" if genes else ""))


def empty_table(path: Path, columns: list[str]) -> pd.DataFrame:
    path.parent.mkdir(parents=True, exist_ok=True)
    df = pd.DataFrame(columns=columns)
    df.to_csv(path, sep="\t", index=False)
    return df


def request_with_retry(method: str, url: str, *, max_tries: int = 3, sleep: float = 1.0, **kwargs):
    last = None
    for attempt in range(1, max_tries + 1):
        try:
            resp = requests.request(method, url, timeout=60, **kwargs)
            if resp.status_code == 200:
                return resp
            last = f"HTTP {resp.status_code}: {resp.text[:300]}"
        except Exception as exc:  # noqa: BLE001
            last = repr(exc)
        time.sleep(sleep * attempt)
    raise RuntimeError(last or "request failed")


def fetch_dgidb(genes: list[str], analysis_level: str, pair_id: str) -> pd.DataFrame:
    cols = [
        "analysis_type",
        "pair_id",
        "gene_symbol",
        "drug_name",
        "drug_concept_id",
        "approved",
        "immunotherapy",
        "anti_neoplastic",
        "interaction_score",
        "evidence_score",
        "interaction_id",
    ]
    if not genes:
        return pd.DataFrame(columns=cols)
    payload = {"query": DGIDB_QUERY, "variables": {"genes": genes, "first": 1000}}
    resp = request_with_retry("POST", DGIDB_URL, json=payload)
    data = resp.json()
    if "errors" in data:
        raise RuntimeError(json.dumps(data["errors"])[:800])
    rows = []
    for node in data.get("data", {}).get("interactions", {}).get("nodes", []) or []:
        gene = node.get("gene") or {}
        drug = node.get("drug") or {}
        rows.append(
            {
                "analysis_type": analysis_level,
                "pair_id": pair_id,
                "gene_symbol": clean_gene(gene.get("name", "")),
                "drug_name": str(drug.get("name", "")).upper(),
                "drug_concept_id": drug.get("conceptId"),
                "approved": drug.get("approved"),
                "immunotherapy": drug.get("immunotherapy"),
                "anti_neoplastic": drug.get("antiNeoplastic"),
                "interaction_score": node.get("interactionScore"),
                "evidence_score": node.get("evidenceScore"),
                "interaction_id": node.get("id"),
            }
        )
    return pd.DataFrame(rows, columns=cols).drop_duplicates()


def fetch_string_network(genes: list[str], analysis_type: str, pair_id: str) -> pd.DataFrame:
    cols = [
        "analysis_type",
        "pair_id",
        "preferredName_A",
        "preferredName_B",
        "stringId_A",
        "stringId_B",
        "score",
        "nscore",
        "fscore",
        "pscore",
        "ascore",
        "escore",
        "dscore",
        "tscore",
    ]
    if len(genes) < 2:
        return pd.DataFrame(columns=cols)
    payload = {
        "identifiers": "\r".join(genes),
        "species": 9606,
        "required_score": 400,
        "caller_identity": "nasal_score_final_network_pharmacology",
    }
    resp = request_with_retry("POST", STRING_NETWORK_URL, data=payload)
    lines = [line for line in resp.text.splitlines() if line.strip()]
    if not lines:
        return pd.DataFrame(columns=cols)
    from io import StringIO

    df = pd.read_csv(StringIO(resp.text), sep="\t")
    if df.empty:
        return pd.DataFrame(columns=cols)
    df.insert(0, "pair_id", pair_id)
    df.insert(0, "analysis_type", analysis_type)
    keep = [c for c in cols if c in df.columns]
    return df[keep].drop_duplicates()


def fetch_string_enrichment(genes: list[str], analysis_type: str, pair_id: str) -> pd.DataFrame:
    cols = [
        "analysis_type",
        "pair_id",
        "category",
        "term",
        "description",
        "number_of_genes",
        "number_of_genes_in_background",
        "ncbiTaxonId",
        "inputGenes",
        "preferredNames",
        "p_value",
        "fdr",
    ]
    if len(genes) < 2:
        return pd.DataFrame(columns=cols)
    payload = {
        "identifiers": "\r".join(genes),
        "species": 9606,
        "caller_identity": "nasal_score_final_network_pharmacology",
    }
    resp = request_with_retry("POST", STRING_ENRICHMENT_URL, data=payload)
    if not resp.text.strip():
        return pd.DataFrame(columns=cols)
    from io import StringIO

    df = pd.read_csv(StringIO(resp.text), sep="\t")
    if df.empty:
        return pd.DataFrame(columns=cols)
    df.insert(0, "pair_id", pair_id)
    df.insert(0, "analysis_type", analysis_type)
    keep = [c for c in cols if c in df.columns]
    return df[keep].drop_duplicates()


def summarize_drugs(dgidb: pd.DataFrame) -> pd.DataFrame:
    cols = [
        "analysis_type",
        "pair_id",
        "drug_name",
        "approved",
        "immunotherapy",
        "anti_neoplastic",
        "n_target_genes",
        "target_gene_list",
        "max_interaction_score",
        "sum_interaction_score",
        "max_evidence_score",
        "n_interactions",
    ]
    if dgidb.empty:
        return pd.DataFrame(columns=cols)
    rows = []
    for keys, sub in dgidb.groupby(["analysis_type", "pair_id", "drug_name"], dropna=False):
        rows.append(
            {
                "analysis_type": keys[0],
                "pair_id": keys[1],
                "drug_name": keys[2],
                "approved": bool(sub["approved"].fillna(False).any()),
                "immunotherapy": bool(sub["immunotherapy"].fillna(False).any()),
                "anti_neoplastic": bool(sub["anti_neoplastic"].fillna(False).any()),
                "n_target_genes": sub["gene_symbol"].nunique(),
                "target_gene_list": ";".join(sorted(sub["gene_symbol"].dropna().unique())),
                "max_interaction_score": pd.to_numeric(sub["interaction_score"], errors="coerce").max(),
                "sum_interaction_score": pd.to_numeric(sub["interaction_score"], errors="coerce").sum(),
                "max_evidence_score": pd.to_numeric(sub["evidence_score"], errors="coerce").max(),
                "n_interactions": len(sub),
            }
        )
    out = pd.DataFrame(rows, columns=cols)
    return out.sort_values(["n_target_genes", "sum_interaction_score", "drug_name"], ascending=[False, False, True])


def ppi_metrics(edges: pd.DataFrame, genes: list[str], analysis_type: str, pair_id: str) -> pd.DataFrame:
    cols = [
        "analysis_type",
        "pair_id",
        "gene_symbol",
        "degree",
        "weighted_degree",
        "betweenness",
        "closeness",
        "clustering",
        "component_size",
    ]
    graph = nx.Graph()
    graph.add_nodes_from(genes)
    for _, row in edges.iterrows():
        a, b = clean_gene(row.get("preferredName_A", "")), clean_gene(row.get("preferredName_B", ""))
        if a and b and a != b:
            graph.add_edge(a, b, weight=float(row.get("score", 0) or 0))
    if graph.number_of_nodes() == 0:
        return pd.DataFrame(columns=cols)
    degree = dict(graph.degree())
    weighted_degree = dict(graph.degree(weight="weight"))
    betweenness = nx.betweenness_centrality(graph, weight="weight", normalized=True) if graph.number_of_edges() else {}
    closeness = nx.closeness_centrality(graph) if graph.number_of_edges() else {}
    clustering = nx.clustering(graph, weight="weight") if graph.number_of_edges() else {}
    comp_size = {}
    for comp in nx.connected_components(graph):
        for node in comp:
            comp_size[node] = len(comp)
    rows = []
    for gene in genes:
        rows.append(
            {
                "analysis_type": analysis_type,
                "pair_id": pair_id,
                "gene_symbol": gene,
                "degree": degree.get(gene, 0),
                "weighted_degree": weighted_degree.get(gene, 0.0),
                "betweenness": betweenness.get(gene, 0.0),
                "closeness": closeness.get(gene, 0.0),
                "clustering": clustering.get(gene, 0.0),
                "component_size": comp_size.get(gene, 1),
            }
        )
    return pd.DataFrame(rows, columns=cols).sort_values(["degree", "weighted_degree"], ascending=[False, False])


def plot_placeholder(path: Path, title: str, message: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.axis("off")
    ax.text(0.5, 0.62, title, ha="center", va="center", fontsize=11, fontweight="bold")
    ax.text(0.5, 0.42, message, ha="center", va="center", fontsize=9, wrap=True)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def plot_ppi(edges: pd.DataFrame, genes: list[str], path: Path, title: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if len(genes) < 2 or edges.empty:
        plot_placeholder(path, title, "No STRING PPI edge at required_score >= 400.")
        return
    graph = nx.Graph()
    graph.add_nodes_from(genes)
    for _, row in edges.iterrows():
        a, b = clean_gene(row.get("preferredName_A", "")), clean_gene(row.get("preferredName_B", ""))
        if a and b and a != b:
            graph.add_edge(a, b, weight=float(row.get("score", 0) or 0))
    if graph.number_of_edges() == 0:
        plot_placeholder(path, title, "No STRING PPI edge at required_score >= 400.")
        return
    fig, ax = plt.subplots(figsize=(max(7, min(13, 0.35 * len(graph.nodes) + 5)), 7))
    pos = nx.spring_layout(graph, seed=13, weight="weight", k=None)
    degrees = dict(graph.degree())
    sizes = [80 + 70 * degrees.get(n, 0) for n in graph.nodes()]
    widths = [0.6 + 2.0 * graph[u][v].get("weight", 0.4) for u, v in graph.edges()]
    nx.draw_networkx_edges(graph, pos, ax=ax, width=widths, alpha=0.45, edge_color="#4c566a")
    nx.draw_networkx_nodes(graph, pos, ax=ax, node_size=sizes, node_color="#4daf4a", edgecolors="#1f2933", linewidths=0.6)
    nx.draw_networkx_labels(graph, pos, ax=ax, font_size=7)
    ax.set_title(title, fontsize=11)
    ax.axis("off")
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def plot_drug_gene(dgidb: pd.DataFrame, path: Path, title: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if dgidb.empty:
        plot_placeholder(path, title, "No DGIdb drug-gene interaction returned.")
        return
    top_drugs = (
        dgidb.groupby("drug_name")["gene_symbol"]
        .nunique()
        .sort_values(ascending=False)
        .head(20)
        .index.tolist()
    )
    sub = dgidb[dgidb["drug_name"].isin(top_drugs)].copy()
    graph = nx.Graph()
    for _, row in sub.iterrows():
        gene = row["gene_symbol"]
        drug = "DRUG:" + row["drug_name"]
        graph.add_node(gene, bipartite="gene")
        graph.add_node(drug, bipartite="drug")
        graph.add_edge(gene, drug, weight=float(row.get("interaction_score", 1) or 1))
    if graph.number_of_edges() == 0:
        plot_placeholder(path, title, "No plottable DGIdb drug-gene interaction returned.")
        return
    fig, ax = plt.subplots(figsize=(max(8, min(15, 0.22 * graph.number_of_nodes() + 6)), 8))
    pos = nx.spring_layout(graph, seed=17, k=0.8)
    node_colors = ["#377eb8" if graph.nodes[n].get("bipartite") == "gene" else "#e41a1c" for n in graph.nodes()]
    node_sizes = [80 + 45 * graph.degree(n) for n in graph.nodes()]
    labels = {n: n.replace("DRUG:", "") for n in graph.nodes()}
    nx.draw_networkx_edges(graph, pos, ax=ax, edge_color="#9aa0a6", alpha=0.35, width=0.8)
    nx.draw_networkx_nodes(graph, pos, ax=ax, node_color=node_colors, node_size=node_sizes, edgecolors="#111827", linewidths=0.4)
    nx.draw_networkx_labels(graph, pos, labels=labels, ax=ax, font_size=6)
    ax.set_title(title, fontsize=11)
    ax.axis("off")
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def plot_top_drugs(summary: pd.DataFrame, path: Path, title: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if summary.empty:
        plot_placeholder(path, title, "No DGIdb drug summary available.")
        return
    sub = summary.head(15).iloc[::-1]
    fig, ax = plt.subplots(figsize=(8, max(4, 0.34 * len(sub) + 1.5)))
    colors = ["#d95f02" if x else "#7570b3" for x in sub["approved"]]
    ax.barh(sub["drug_name"], sub["n_target_genes"], color=colors)
    ax.set_xlabel("Target genes")
    ax.set_title(title, fontsize=11)
    ax.tick_params(axis="y", labelsize=7)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def plot_overview_bar(df: pd.DataFrame, path: Path, title: str, xcol: str, ycol: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if df.empty:
        plot_placeholder(path, title, "No data.")
        return
    sub = df.sort_values(ycol, ascending=False).head(30).iloc[::-1]
    fig, ax = plt.subplots(figsize=(9, max(5, 0.25 * len(sub) + 1.5)))
    ax.barh(sub[xcol], sub[ycol], color="#4c78a8")
    ax.set_xlabel(ycol)
    ax.set_title(title, fontsize=11)
    ax.tick_params(axis="y", labelsize=6)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def read_inputs(score_dir: Path) -> pd.DataFrame:
    main = pd.read_csv(score_dir / "final_twas_supported_pair_gene_records.tsv", sep="\t")
    needed = {"pair_id", "gene_symbol", "gene_grade", "phenotype_domain"}
    missing = needed - set(main.columns)
    if missing:
        raise ValueError(f"main table missing columns: {sorted(missing)}")
    main["gene_symbol"] = main["gene_symbol"].map(clean_gene)
    main = main.dropna(subset=["pair_id", "gene_symbol"]).drop_duplicates()
    return main


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--score-dir", required=True, help="Directory containing prioritized pairwise gene tables.")
    parser.add_argument("--out-dir", required=True, help="Output directory for exploratory drug-target annotations.")
    parser.add_argument("--sleep", type=float, default=0.25)
    args = parser.parse_args()

    score_dir = Path(args.score_dir)
    out_dir = Path(args.out_dir)
    dirs = {
        "inputs": out_dir / "input_gene_sets",
        "dgidb": out_dir / "tables" / "dgidb_interactions",
        "ppi": out_dir / "tables" / "string_ppi",
        "enrich": out_dir / "tables" / "string_enrichment",
        "fig_drug": out_dir / "figures" / "drug_gene_network",
        "fig_ppi": out_dir / "figures" / "ppi_network",
        "fig_top_drugs": out_dir / "figures" / "top_drugs",
        "fig_summary": out_dir / "figures" / "summary",
        "audit": out_dir / "audit",
    }
    for d in dirs.values():
        d.mkdir(parents=True, exist_ok=True)

    main_df = read_inputs(score_dir)
    pair_domains = main_df.drop_duplicates("pair_id").set_index("pair_id")["phenotype_domain"].to_dict()
    pair_ids = sorted(main_df["pair_id"].unique())

    all_dgidb, all_drugs, all_ppi, all_metrics, all_hubs, all_enrich = [], [], [], [], [], []
    manifest, audit_rows = [], []

    for i, pair_id in enumerate(pair_ids, start=1):
        pair_df = main_df[main_df["pair_id"] == pair_id]
        analysis_sets = {
            "GeneA": sorted(pair_df.loc[pair_df["gene_grade"].eq("Gene-A"), "gene_symbol"].dropna().unique()),
            "GeneA_or_B": sorted(pair_df["gene_symbol"].dropna().unique()),
        }
        for analysis_type, genes in analysis_sets.items():
            tag = f"{safe_name(pair_id)}.{analysis_type}"
            print(f"[{i}/{len(pair_ids)}] {pair_id} {analysis_type} n={len(genes)}", flush=True)
            write_gene_list(dirs["inputs"] / f"{tag}.genes.txt", genes)
            manifest.append(
                {
                    "pair_id": pair_id,
                    "phenotype_domain": pair_domains.get(pair_id, ""),
                    "analysis_type": analysis_type,
                    "gene_file": str(dirs["inputs"] / f"{tag}.genes.txt"),
                    "n_genes": len(genes),
                    "gene_list": ";".join(genes),
                    "run_dgidb": len(genes) >= 1,
                    "run_string_ppi": len(genes) >= 2,
                    "run_string_enrichment": len(genes) >= 2,
                }
            )

            dgidb_err = ppi_err = enrich_err = ""
            try:
                dgidb = fetch_dgidb(genes, analysis_type, pair_id)
            except Exception as exc:  # noqa: BLE001
                dgidb_err = repr(exc)
                dgidb = pd.DataFrame()
            if dgidb.empty:
                dgidb = empty_table(
                    dirs["dgidb"] / f"{tag}.dgidb_interactions.tsv",
                    [
                        "analysis_type",
                        "pair_id",
                        "gene_symbol",
                        "drug_name",
                        "drug_concept_id",
                        "approved",
                        "immunotherapy",
                        "anti_neoplastic",
                        "interaction_score",
                        "evidence_score",
                        "interaction_id",
                    ],
                )
            else:
                dgidb.to_csv(dirs["dgidb"] / f"{tag}.dgidb_interactions.tsv", sep="\t", index=False)
            drug_sum = summarize_drugs(dgidb)
            drug_sum.to_csv(dirs["dgidb"] / f"{tag}.drug_summary.tsv", sep="\t", index=False)

            time.sleep(args.sleep)
            try:
                ppi = fetch_string_network(genes, analysis_type, pair_id)
            except Exception as exc:  # noqa: BLE001
                ppi_err = repr(exc)
                ppi = pd.DataFrame()
            if ppi.empty:
                ppi = empty_table(
                    dirs["ppi"] / f"{tag}.string_ppi_edges.tsv",
                    [
                        "analysis_type",
                        "pair_id",
                        "preferredName_A",
                        "preferredName_B",
                        "stringId_A",
                        "stringId_B",
                        "score",
                        "nscore",
                        "fscore",
                        "pscore",
                        "ascore",
                        "escore",
                        "dscore",
                        "tscore",
                    ],
                )
            else:
                ppi.to_csv(dirs["ppi"] / f"{tag}.string_ppi_edges.tsv", sep="\t", index=False)
            metrics = ppi_metrics(ppi, genes, analysis_type, pair_id)
            metrics.to_csv(dirs["ppi"] / f"{tag}.node_metrics.tsv", sep="\t", index=False)
            hubs = metrics.sort_values(["degree", "weighted_degree"], ascending=[False, False]).head(10)
            hubs.to_csv(dirs["ppi"] / f"{tag}.hub_genes.tsv", sep="\t", index=False)

            time.sleep(args.sleep)
            try:
                enrich = fetch_string_enrichment(genes, analysis_type, pair_id)
            except Exception as exc:  # noqa: BLE001
                enrich_err = repr(exc)
                enrich = pd.DataFrame()
            if enrich.empty:
                enrich = empty_table(
                    dirs["enrich"] / f"{tag}.string_enrichment.tsv",
                    [
                        "analysis_type",
                        "pair_id",
                        "category",
                        "term",
                        "description",
                        "number_of_genes",
                        "number_of_genes_in_background",
                        "ncbiTaxonId",
                        "inputGenes",
                        "preferredNames",
                        "p_value",
                        "fdr",
                    ],
                )
            else:
                enrich.to_csv(dirs["enrich"] / f"{tag}.string_enrichment.tsv", sep="\t", index=False)

            plot_ppi(ppi, genes, dirs["fig_ppi"] / f"{tag}.ppi_network.pdf", f"{pair_id} | {analysis_type} STRING PPI")
            plot_drug_gene(dgidb, dirs["fig_drug"] / f"{tag}.drug_gene_network.pdf", f"{pair_id} | {analysis_type} DGIdb")
            plot_top_drugs(drug_sum, dirs["fig_top_drugs"] / f"{tag}.top_drugs.pdf", f"{pair_id} | {analysis_type} top drugs")

            all_dgidb.append(dgidb)
            all_drugs.append(drug_sum)
            all_ppi.append(ppi)
            all_metrics.append(metrics)
            all_hubs.append(hubs)
            all_enrich.append(enrich)
            audit_rows.append(
                {
                    "pair_id": pair_id,
                    "analysis_type": analysis_type,
                    "n_input_genes": len(genes),
                    "n_dgidb_interactions": len(dgidb),
                    "n_drugs": drug_sum["drug_name"].nunique() if not drug_sum.empty else 0,
                    "n_string_edges": len(ppi),
                    "n_string_nodes_with_degree_gt0": int((metrics["degree"] > 0).sum()) if not metrics.empty else 0,
                    "n_string_enrichment_terms": len(enrich),
                    "dgidb_error": dgidb_err,
                    "string_ppi_error": ppi_err,
                    "string_enrichment_error": enrich_err,
                }
            )

    manifest_df = pd.DataFrame(manifest)
    audit_df = pd.DataFrame(audit_rows)
    manifest_df.to_csv(out_dir / "network_pharmacology_input_manifest.tsv", sep="\t", index=False)
    audit_df.to_csv(dirs["audit"] / "network_pharmacology_audit.tsv", sep="\t", index=False)

    def concat_or_empty(parts: list[pd.DataFrame]) -> pd.DataFrame:
        parts = [p for p in parts if p is not None and not p.empty]
        return pd.concat(parts, ignore_index=True).drop_duplicates() if parts else pd.DataFrame()

    dgidb_all = concat_or_empty(all_dgidb)
    drug_all = concat_or_empty(all_drugs)
    ppi_all = concat_or_empty(all_ppi)
    metric_all = concat_or_empty(all_metrics)
    hub_all = concat_or_empty(all_hubs)
    enrich_all = concat_or_empty(all_enrich)
    dgidb_all.to_csv(out_dir / "all_dgidb_interactions.tsv", sep="\t", index=False)
    drug_all.to_csv(out_dir / "all_drug_summary.tsv", sep="\t", index=False)
    ppi_all.to_csv(out_dir / "all_string_ppi_edges.tsv", sep="\t", index=False)
    metric_all.to_csv(out_dir / "all_node_metrics.tsv", sep="\t", index=False)
    hub_all.to_csv(out_dir / "all_hub_genes.tsv", sep="\t", index=False)
    enrich_all.to_csv(out_dir / "all_string_enrichment.tsv", sep="\t", index=False)

    summary = audit_df.merge(manifest_df[["pair_id", "analysis_type", "phenotype_domain", "gene_list"]], on=["pair_id", "analysis_type"], how="left")
    summary.to_csv(out_dir / "pair_analysis_summary.tsv", sep="\t", index=False)

    if not drug_all.empty:
        rec_rows = []
        for (analysis_type, drug), sub in drug_all.groupby(["analysis_type", "drug_name"]):
            rec_rows.append(
                {
                    "analysis_type": analysis_type,
                    "drug_name": drug,
                    "n_pairs": sub["pair_id"].nunique(),
                    "pair_list": ";".join(sorted(sub["pair_id"].unique())),
                    "n_target_genes_total": len(set(";".join(sub["target_gene_list"].dropna()).split(";")) - {""}),
                    "target_gene_list": ";".join(sorted(set(";".join(sub["target_gene_list"].dropna()).split(";")) - {""})),
                    "approved_any": bool(sub["approved"].fillna(False).any()),
                    "max_target_genes_per_pair": sub["n_target_genes"].max(),
                }
            )
        drug_rec = pd.DataFrame(rec_rows).sort_values(["analysis_type", "n_pairs", "max_target_genes_per_pair"], ascending=[True, False, False])
    else:
        drug_rec = pd.DataFrame()
    drug_rec.to_csv(out_dir / "drug_recurrence_summary.tsv", sep="\t", index=False)

    if not dgidb_all.empty:
        gene_drug = (
            dgidb_all.groupby(["analysis_type", "gene_symbol"])
            .agg(
                n_drugs=("drug_name", "nunique"),
                n_approved_drugs=("approved", lambda s: int(pd.Series(s).fillna(False).sum())),
                drug_list=("drug_name", lambda s: ";".join(sorted(set(s.dropna())))),
                pair_list=("pair_id", lambda s: ";".join(sorted(set(s.dropna())))),
            )
            .reset_index()
            .sort_values(["analysis_type", "n_drugs"], ascending=[True, False])
        )
    else:
        gene_drug = pd.DataFrame()
    gene_drug.to_csv(out_dir / "gene_druggability_summary.tsv", sep="\t", index=False)

    if not summary.empty:
        plot_overview_bar(
            summary[summary["analysis_type"].eq("GeneA_or_B")].assign(short_pair=lambda x: x["pair_id"].str.replace("ALLERGIC_RHINITIS_GCST90038664__", "AR__", regex=False).str.replace("NASAL_POLYPS_GCST90018883__", "NP__", regex=False)),
            dirs["fig_summary"] / "GeneA_or_B_string_edge_counts_by_pair.pdf",
            "GeneA+B STRING PPI edge counts by pair",
            "short_pair",
            "n_string_edges",
        )
        plot_overview_bar(
            summary[summary["analysis_type"].eq("GeneA_or_B")].assign(short_pair=lambda x: x["pair_id"].str.replace("ALLERGIC_RHINITIS_GCST90038664__", "AR__", regex=False).str.replace("NASAL_POLYPS_GCST90018883__", "NP__", regex=False)),
            dirs["fig_summary"] / "GeneA_or_B_drug_counts_by_pair.pdf",
            "GeneA+B DGIdb drug counts by pair",
            "short_pair",
            "n_drugs",
        )

    errors = int(
        audit_df[["dgidb_error", "string_ppi_error", "string_enrichment_error"]]
        .fillna("")
        .apply(lambda col: col.astype(str).str.len() > 0)
        .any(axis=1)
        .sum()
    )
    warnings = int((manifest_df["n_genes"] < 2).sum())
    audit_summary = [
        "network_pharmacology_audit_summary",
        f"analysis_timestamp={pd.Timestamp.now().isoformat()}",
        f"n_pairs={len(pair_ids)}",
        f"n_pair_analysis_runs={len(manifest_df)}",
        f"n_GeneA_runs={int(manifest_df['analysis_type'].eq('GeneA').sum())}",
        f"n_GeneA_or_B_runs={int(manifest_df['analysis_type'].eq('GeneA_or_B').sum())}",
        f"n_runs_with_zero_genes={int(manifest_df['n_genes'].eq(0).sum())}",
        f"n_runs_with_lt2_genes={int((manifest_df['n_genes'] < 2).sum())}",
        f"n_total_dgidb_interactions={len(dgidb_all)}",
        f"n_total_unique_drugs={dgidb_all['drug_name'].nunique() if not dgidb_all.empty else 0}",
        f"n_total_string_edges={len(ppi_all)}",
        f"n_total_string_enrichment_terms={len(enrich_all)}",
        f"errors={errors}",
        f"warnings={warnings}",
        "usable_for_supplement=yes" if errors == 0 else "usable_for_supplement=check_failed_api_rows",
    ]
    (dirs["audit"] / "audit_summary.txt").write_text("\n".join(audit_summary) + "\n")

    readme = f"""# score_final pairwise network pharmacology

## 分析口径

本目录基于 `score_final` 最终 TWAS 评分体系做网络药理学整理，不重跑 TWAS/MAGMA/PLACO/CPASSOC/coloc/SuSiE。

- Gene-A = Locus-A + MAGMA Bonferroni positive + FUSION GTEx v8 49-tissue candidate TWAS supported。
- Gene-B = Locus-B + MAGMA Bonferroni positive + FUSION GTEx v8 49-tissue candidate TWAS supported。
- TWAS-supported = 同一 disease pair 内，任意 GTEx v8 tissue 的 pair x tissue BH-FDR <= 0.05。
- 本分析每个疾病对做两套输入：`GeneA` 和 `GeneA_or_B`。

## 数据库和方法

- DGIdb GraphQL API：提取 drug-gene interactions，并汇总每个 drug 命中的 target genes。
- STRING API human species 9606：提取 PPI network，`required_score=400`。
- STRING enrichment API：补充每个输入基因集的 STRING 功能富集注释。
- 图形使用 matplotlib/networkx 生成 PDF。

## 主要输出

- `network_pharmacology_input_manifest.tsv`：每个 pair 每套输入的基因数和基因列表。
- `pair_analysis_summary.tsv`：每个 pair 的 DGIdb 药物数、STRING 边数、富集 term 数。
- `all_dgidb_interactions.tsv`：所有 drug-gene interaction 明细。
- `all_drug_summary.tsv`：每个 pair 的 drug 层面汇总。
- `drug_recurrence_summary.tsv`：跨 pair 复现出现的 drug。
- `gene_druggability_summary.tsv`：每个基因对应的药物可靶向性汇总。
- `all_string_ppi_edges.tsv`、`all_node_metrics.tsv`、`all_hub_genes.tsv`：PPI 边、节点指标和 hub genes。
- `all_string_enrichment.tsv`：STRING 功能富集注释。
- `figures/`：每个 pair 的 PPI、drug-gene network、top drug 图，以及总体汇总图。

## 运行汇总

- disease pairs: {len(pair_ids)}
- pair-analysis runs: {len(manifest_df)}，其中 GeneA {int(manifest_df['analysis_type'].eq('GeneA').sum())} 个，GeneA+B {int(manifest_df['analysis_type'].eq('GeneA_or_B').sum())} 个。
- GeneA 输入为空的 pair 数: {int((manifest_df['analysis_type'].eq('GeneA') & manifest_df['n_genes'].eq(0)).sum())}
- DGIdb interactions: {len(dgidb_all)}
- unique drugs: {dgidb_all['drug_name'].nunique() if not dgidb_all.empty else 0}
- STRING PPI edges: {len(ppi_all)}
- STRING enrichment terms: {len(enrich_all)}
- API/audit errors: {errors}
- warnings: {warnings}

## 解释限制

网络药理学结果是功能和药物可靶向性注释，不重新定义 Gene-A/Gene-B，也不构成药物疗效或因果证据。GeneA-only 结果更严格；GeneA+B 结果覆盖主 TWAS-supported prioritized genes，更适合作为补充候选药物/通路线索。
"""
    (out_dir / "README_network_pharmacology_final.md").write_text(readme)
    print("\n".join(audit_summary))
    return 0


if __name__ == "__main__":
    sys.exit(main())
