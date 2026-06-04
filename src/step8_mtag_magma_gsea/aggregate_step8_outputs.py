#!/usr/bin/env python3
from __future__ import annotations

import os
from pathlib import Path

import pandas as pd


ROOT = Path(os.environ.get("PROJECT_ROOT", "/platform_data/p_user/p010/phase0"))
MAGMA = ROOT / "results/step8_magma"
PATHWAY = ROOT / "results/step8_pathway"


def concat_files(files, out):
    frames = []
    for path in files:
        if path.exists() and path.stat().st_size > 0:
            frames.append(pd.read_csv(path, sep="\t", dtype=str))
    if frames:
        df = pd.concat(frames, ignore_index=True)
    else:
        df = pd.DataFrame()
    out.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out, sep="\t", index=False)
    print(out, len(df))


def main():
    concat_files(sorted((MAGMA / "per_pair").glob("*.gene_results.tsv")), MAGMA / "pair_specific_gene_results.tsv")
    per = PATHWAY / "per_pair"
    concat_files(sorted(per.glob("*.pair_specific_gsea_go_bp.tsv")), PATHWAY / "pair_specific_gsea_go_bp.tsv")
    concat_files(sorted(per.glob("*.pair_specific_gsea_kegg.tsv")), PATHWAY / "pair_specific_gsea_kegg.tsv")
    concat_files(sorted(per.glob("*.pair_specific_gsea_reactome.tsv")), PATHWAY / "pair_specific_gsea_reactome.tsv")


if __name__ == "__main__":
    main()
