#!/usr/bin/env python3
from __future__ import annotations

from pathlib import Path
import os
import re

import pandas as pd
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests

project_root = os.environ.get("NASAL_PROJECT_ROOT", "").strip()
if not project_root:
    raise SystemExit("Set NASAL_PROJECT_ROOT to the directory containing results_end.")
ROOT = Path(project_root)
FINAL = ROOT / "results_end/final_gtexv8_twas_total_domain_fuma_enrichment"
OUT = ROOT / "results_end/integrated_fuma_tissue_cell_ora"
SUPP = OUT / "cell_marker_ora/gene_grade_background_ora"
GMT_DIR = OUT / "cell_marker_ora/marker_gmt"
DOMAINS = ["asthma_lower_airway", "atopic_allergic", "autoimmune_IBD", "ENT_infection"]
SOURCES = {
    "PanglaoDB": GMT_DIR / "panglaodb_human_cell_markers.gmt",
    "CellMarker": GMT_DIR / "cellmarker_human_cell_markers.gmt",
    "TabulaSapiens": GMT_DIR / "tabula_sapiens_cell_markers.gmt",
    "airway_nasal_literature": GMT_DIR / "airway_nasal_literature_cell_markers.gmt",
    "combined": GMT_DIR / "combined_cell_markers_nonredundant.gmt",
    "airway_immune_focused": GMT_DIR / "airway_immune_focused_cell_markers.gmt",
}


def clean_gene(x) -> str:
    x = str(x).strip().upper()
    return "" if x in {"", "NA", "NAN", "NONE"} else x


def uniq(vals) -> list[str]:
    return sorted({g for g in (clean_gene(x) for x in vals) if g})


def read_list(path: Path) -> list[str]:
    return uniq(path.read_text().splitlines())


def write_list(path: Path, genes) -> int:
    vals = uniq(genes)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(vals) + ("\n" if vals else ""))
    return len(vals)


def read_tsv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t")


def parse_gmt(path: Path, source: str) -> list[dict]:
    terms = []
    for line in path.read_text().splitlines():
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 3:
            continue
        genes = uniq(parts[2:])
        if len(genes) >= 3:
            terms.append({"marker_source": source, "term_name": parts[0], "source_description": parts[1], "genes": genes})
    return terms


def term_meta(term: str) -> tuple[str, str]:
    lo = str(term).lower().replace("_", " ")
    cell = str(term)
    tissue = ""
    for p in ["nasal", "airway", "bronchial", "lung", "ciliated", "basal", "secretory", "goblet", "club", "ionocyte", "epithelial", "cd4", "th2", "activated t", "t cell", "b cell", "plasma", "mast", "eosinophil", "neutrophil", "monocyte", "macrophage", "dendritic", "fibroblast", "endothelial", "smooth muscle", "keratinocyte", "intestinal", "colon"]:
        if p in lo:
            cell = p
            break
    for t in ["nasal", "airway", "bronchial", "lung", "blood", "colon", "intestin", "skin", "spleen"]:
        if t in lo:
            tissue = t
            break
    return tissue, cell


def run_ora(analysis_level: str, analysis_id: str, comparison: str, input_genes: list[str], background_genes: list[str], terms: list[dict], source: str) -> pd.DataFrame:
    inp = set(input_genes)
    bg = set(background_genes)
    rows = []
    for term in terms:
        marker = set(term["genes"]) & bg
        term_size = len(marker)
        if term_size < 5 or term_size > 500:
            continue
        a = len(inp & marker)
        b = len(inp) - a
        c = term_size - a
        d = len(bg) - len(inp) - c
        if min(a, b, c, d) < 0:
            continue
        odds, p = fisher_exact([[a, b], [c, d]], alternative="greater")
        tissue, cell = term_meta(term["term_name"])
        rows.append({
            "analysis_level": analysis_level,
            "analysis_id": analysis_id,
            "comparison": comparison,
            "marker_source": source,
            "tissue": tissue,
            "cell_type": cell,
            "term_name": term["term_name"],
            "input_gene_count": len(inp),
            "background_gene_count": len(bg),
            "term_size_in_background": term_size,
            "overlap_count": a,
            "gene_ratio": f"{a}/{len(inp)}",
            "background_ratio": f"{term_size}/{len(bg)}",
            "odds_ratio": odds,
            "p_value": p,
            "overlap_genes": ";".join(sorted(inp & marker)),
            "marker_genes_in_background": ";".join(sorted(marker)),
        })
    df = pd.DataFrame(rows)
    if not df.empty:
        df["FDR"] = multipletests(df["p_value"].fillna(1), method="fdr_bh")[1]
        df = df.sort_values(["FDR", "p_value", "term_name"])
    else:
        df["FDR"] = []
    return df


def main() -> None:
    SUPP.mkdir(parents=True, exist_ok=True)
    (SUPP / "input_gene_sets").mkdir(exist_ok=True)
    (SUPP / "ora_results").mkdir(exist_ok=True)

    gene_a = read_tsv(FINAL / "final_gene_tables/final_geneA_pair_gene_records.tsv")
    gene_b = read_tsv(FINAL / "final_gene_tables/final_geneB_pair_gene_records.tsv")
    pair_records = read_tsv(FINAL / "final_gene_tables/final_twas_supported_pair_gene_records.tsv")
    ab_bg = read_list(FINAL / "input_gene_sets/background_AB_locus_3782_unique_genes.txt")
    magma_bg = read_list(FINAL / "input_gene_sets/background_MAGMA_candidate_205_unique_genes.txt")

    gene_a_set = uniq(gene_a["gene_symbol"])
    gene_b_set = uniq(gene_b["gene_symbol"])
    gene_ab_set = uniq(gene_a_set + gene_b_set)
    primary_set = uniq(pair_records["gene_symbol"])
    remaining_ab = sorted(set(ab_bg) - set(primary_set))
    remaining_magma = sorted(set(magma_bg) - set(primary_set))

    tasks = []
    def add(level, aid, comp, genes, bg, background_name):
        ip = SUPP / "input_gene_sets" / f"{level}_{aid}_{comp}_input.txt"
        bp = SUPP / "input_gene_sets" / f"{level}_{aid}_{comp}_{background_name}_background.txt"
        write_list(ip, genes)
        write_list(bp, bg)
        tasks.append({
            "analysis_level": level,
            "analysis_id": aid,
            "comparison": comp,
            "background_name": background_name,
            "input_gene_file": str(ip),
            "n_input_genes": len(uniq(genes)),
            "background_gene_file": str(bp),
            "n_background_genes": len(uniq(bg)),
            "n_input_not_in_background": len(set(uniq(genes)) - set(uniq(bg))),
            "input_not_in_background": ";".join(sorted(set(uniq(genes)) - set(uniq(bg)))),
            "run_ORA": len(uniq(genes)) >= 5 and len(set(uniq(genes)) - set(uniq(bg))) == 0,
            "interpretation": "primary_gene_grade_support" if comp in {"GeneA", "GeneB", "GeneA_or_B"} else "background_remainder_negative_control",
        })

    for comp, genes in [("GeneA", gene_a_set), ("GeneB", gene_b_set), ("GeneA_or_B", gene_ab_set), ("remaining_AB_background_not_TWAS_primary", remaining_ab), ("remaining_MAGMA_background_not_TWAS_primary", remaining_magma)]:
        add("global", "all_pairs", comp, genes, ab_bg, "AB1456")
        if set(genes).issubset(set(magma_bg)):
            add("global", "all_pairs", comp, genes, magma_bg, "MAGMA205")

    for d in DOMAINS:
        d_ab = read_list(FINAL / f"input_gene_sets/domain_backgrounds/background_AB_{d}_unique_genes.txt")
        d_mag = read_list(FINAL / f"input_gene_sets/domain_backgrounds_conditional_MAGMA/background_MAGMA_{d}_unique_genes.txt")
        da = uniq(gene_a.loc[gene_a["phenotype_domain"].eq(d), "gene_symbol"])
        db = uniq(gene_b.loc[gene_b["phenotype_domain"].eq(d), "gene_symbol"])
        dab = uniq(da + db)
        dprimary = uniq(pair_records.loc[pair_records["phenotype_domain"].eq(d), "gene_symbol"])
        d_rem_ab = sorted(set(d_ab) - set(dprimary))
        d_rem_mag = sorted(set(d_mag) - set(dprimary))
        for comp, genes in [("GeneA", da), ("GeneB", db), ("GeneA_or_B", dab), ("remaining_domain_AB_background_not_TWAS_primary", d_rem_ab), ("remaining_domain_MAGMA_background_not_TWAS_primary", d_rem_mag)]:
            add("domain", d, comp, genes, d_ab, f"{d}_AB")
            if set(genes).issubset(set(d_mag)):
                add("domain", d, comp, genes, d_mag, f"{d}_MAGMA")

    manifest = pd.DataFrame(tasks)
    manifest.to_csv(SUPP / "gene_grade_background_ORA_manifest.tsv", sep="\t", index=False)

    all_terms = {source: parse_gmt(path, source) for source, path in SOURCES.items()}
    all_rows = []
    audit = []
    for _, task in manifest.iterrows():
        if not bool(task["run_ORA"]):
            audit.append({**task.to_dict(), "marker_source": "all", "n_tested_terms": 0, "n_FDR05": 0, "n_nominal_p05": 0})
            continue
        inp = read_list(Path(task["input_gene_file"]))
        bg = read_list(Path(task["background_gene_file"]))
        for source, terms in all_terms.items():
            df = run_ora(task.analysis_level, task.analysis_id, task.comparison, inp, bg, terms, source)
            out = SUPP / "ora_results" / f"{task.analysis_level}_{task.analysis_id}_{task.comparison}_vs_{task.background_name}_{source}.tsv"
            df.to_csv(out, sep="\t", index=False)
            if not df.empty:
                all_rows.append(df)
            audit.append({**task.to_dict(), "marker_source": source, "n_tested_terms": len(df), "n_FDR05": int((df["FDR"] <= 0.05).sum()) if not df.empty else 0, "n_nominal_p05": int((df["p_value"] < 0.05).sum()) if not df.empty else 0})

    all_df = pd.concat(all_rows, ignore_index=True) if all_rows else pd.DataFrame()
    all_df.to_csv(SUPP / "gene_grade_background_cell_marker_ORA_all.tsv", sep="\t", index=False)
    if not all_df.empty:
        all_df[all_df["FDR"] <= 0.05].to_csv(SUPP / "gene_grade_background_cell_marker_ORA_FDR05.tsv", sep="\t", index=False)
        all_df.sort_values(["analysis_level", "analysis_id", "comparison", "background_gene_count", "FDR", "p_value"]).groupby(["analysis_level", "analysis_id", "comparison", "background_gene_count"], as_index=False).head(10).to_csv(SUPP / "gene_grade_background_cell_marker_ORA_top10.tsv", sep="\t", index=False)
    else:
        pd.DataFrame().to_csv(SUPP / "gene_grade_background_cell_marker_ORA_FDR05.tsv", sep="\t", index=False)
        pd.DataFrame().to_csv(SUPP / "gene_grade_background_cell_marker_ORA_top10.tsv", sep="\t", index=False)
    audit_df = pd.DataFrame(audit)
    audit_df.to_csv(SUPP / "gene_grade_background_cell_marker_ORA_audit.tsv", sep="\t", index=False)

    summary = []
    summary.append(f"GeneA unique genes: {len(gene_a_set)}")
    summary.append(f"GeneB unique genes: {len(gene_b_set)}")
    summary.append(f"GeneA_or_B unique genes: {len(gene_ab_set)}")
    summary.append(f"AB background unique genes: {len(ab_bg)}")
    summary.append(f"MAGMA candidate background unique genes: {len(magma_bg)}")
    summary.append(f"remaining AB background not TWAS primary: {len(remaining_ab)}")
    summary.append(f"remaining MAGMA background not TWAS primary: {len(remaining_magma)}")
    if not all_df.empty:
        summary.append(f"total ORA rows: {len(all_df)}")
        summary.append(f"FDR<=0.05 rows: {int((all_df['FDR'] <= 0.05).sum())}")
        summary.append("FDR<=0.05 by comparison/source:")
        sig = all_df[all_df["FDR"] <= 0.05]
        if sig.empty:
            summary.append("none")
        else:
            tab = sig.groupby(["analysis_level", "analysis_id", "comparison", "marker_source"]).size().reset_index(name="n_FDR05").sort_values("n_FDR05", ascending=False)
            summary.append(tab.to_string(index=False))
        summary.append("Top global GeneA/GeneB combined nominal terms:")
        top = all_df[(all_df["analysis_level"] == "global") & (all_df["comparison"].isin(["GeneA", "GeneB", "GeneA_or_B"])) & (all_df["marker_source"] == "combined")].sort_values(["FDR", "p_value"]).head(20)
        summary.append(top[["comparison", "background_gene_count", "term_name", "overlap_count", "p_value", "FDR", "overlap_genes"]].to_string(index=False))
    (SUPP / "README_gene_grade_background_cell_marker_ORA.md").write_text(
        "# Gene-grade and Background Cell Marker ORA\n\n"
        "This is a supplementary sensitivity analysis requested after the integrated branch. It reuses the existing marker GMT files and does not rerun TWAS/MAGMA/FUSION/FUMA Cell Type.\n\n"
        "Analyses include Gene-A, Gene-B, GeneA-or-B, AB-background remainder, and MAGMA-background remainder. Backgrounds include global/domain A/B unique gene backgrounds and, where the input is a subset, MAGMA candidate backgrounds.\n\n"
        "The background remainder analyses are negative-control/descriptive checks, not core-gene enrichment claims.\n\n"
        "```text\n" + "\n".join(summary) + "\n```\n"
    )
    (SUPP / "summary_gene_grade_background_cell_marker_ORA.txt").write_text("\n".join(summary) + "\n")
    print("\n".join(summary))


if __name__ == "__main__":
    main()
