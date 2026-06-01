#!/usr/bin/env python3
from pathlib import Path

import pandas as pd


ROOT = Path("/platform_data/p_user/p010/phase0")
inp = ROOT / "results/step8_pathway/top3_pathways_by_pair_method.tsv"
out = ROOT / "results/step8_pathway/top3_pathways_by_pair_method.xlsx"
df = pd.read_csv(inp, sep="\t", dtype=str)
with pd.ExcelWriter(out, engine="openpyxl") as writer:
    df.to_excel(writer, sheet_name="all_top3", index=False)
    for method, group in df.groupby("method", sort=True):
        group.to_excel(writer, sheet_name=method[:31], index=False)
print(out)
