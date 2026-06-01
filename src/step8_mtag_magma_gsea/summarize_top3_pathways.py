#!/usr/bin/env python3
from pathlib import Path

import pandas as pd


ROOT = Path("/platform_data/p_user/p010/phase0")
FILES = {
    "GSEA_GO_BP": "results/step8_pathway/pair_specific_gsea_go_bp.tsv",
    "GSEA_KEGG": "results/step8_pathway/pair_specific_gsea_kegg.tsv",
    "GSEA_Reactome": "results/step8_pathway/pair_specific_gsea_reactome.tsv",
    "ORA_GO_BP": "results/step8_pathway/pair_specific_ora_go_bp.tsv",
    "ORA_KEGG": "results/step8_pathway/pair_specific_ora_kegg.tsv",
    "ORA_Reactome": "results/step8_pathway/pair_specific_ora_reactome.tsv",
}


def main():
    rows = []
    for method, rel in FILES.items():
        path = ROOT / rel
        if not path.exists() or path.stat().st_size == 0:
            continue
        df = pd.read_csv(path, sep="\t", dtype=str, low_memory=False)
        if df.empty:
            continue
        pcol = "p.adjust" if "p.adjust" in df.columns else "pvalue"
        df[pcol] = pd.to_numeric(df[pcol], errors="coerce")
        df["_pvalue_num"] = pd.to_numeric(df.get("pvalue", df[pcol]), errors="coerce")
        desc_col = "pathway_description" if "pathway_description" in df.columns else "pathway_id"
        id_col = "pathway_id" if "pathway_id" in df.columns else desc_col
        df = df.sort_values(["pair_id", pcol, "_pvalue_num"], na_position="last")
        for pair_id, group in df.groupby("pair_id", sort=True):
            for rank, (_, r) in enumerate(group.head(3).iterrows(), start=1):
                rows.append(
                    {
                        "method": method,
                        "pair_id": pair_id,
                        "trait1": r.get("trait1", ""),
                        "trait2": r.get("trait2", ""),
                        "rank": rank,
                        "pathway_id": r.get(id_col, ""),
                        "pathway_description": r.get(desc_col, ""),
                        "pvalue": r.get("pvalue", ""),
                        "p_adjust": r.get("p.adjust", ""),
                        "NES": r.get("NES", ""),
                        "setSize": r.get("setSize", ""),
                        "overlap": r.get("overlap", ""),
                        "overlap_genes": r.get("overlap_genes", ""),
                    }
                )
    out = pd.DataFrame(rows)
    outpath = ROOT / "results/step8_pathway/top3_pathways_by_pair_method.tsv"
    outpath.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(outpath, sep="\t", index=False)
    print(outpath)
    print("rows", len(out))
    if not out.empty:
        print("methods", ",".join(sorted(out["method"].unique())))
        print("pairs", out["pair_id"].nunique())
        print(out.head(30).to_csv(sep="\t", index=False))


if __name__ == "__main__":
    main()
