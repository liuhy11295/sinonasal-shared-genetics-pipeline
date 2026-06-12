from pathlib import Path
import os

import pandas as pd


project_root = os.environ.get("PROJECT_ROOT", "").strip()
if not project_root:
    raise SystemExit("Set PROJECT_ROOT to the analysis project root.")
root = Path(project_root) / "results"
gqc_path = root / "qc/global_qc_report.tsv"

placo = pd.read_csv(root / "placo_genomewide/placo_genomewide_significant_snps.tsv", sep="\t")
cpassoc = pd.read_csv(root / "cpassoc_genomewide/cpassoc_genomewide_significant_snps.tsv", sep="\t")
gqc = pd.read_csv(gqc_path, sep="\t")
metrics = dict(zip(gqc["metric"].astype(str), gqc["value"].astype(str)))

metrics.update(
    {
        "PLACO_genomewide_significant_rows_total_P5e8": str(len(placo)),
        "PLACO_genomewide_P_lt_5e8_SNPs": str(int((placo["PLACO_p"].astype(float) < 5e-8).sum())),
        "CPASSOC_genomewide_significant_rows_total_P5e8": str(len(cpassoc)),
        "CPASSOC_genomewide_any_P_lt_5e8_SNPs": str(
            int(((cpassoc["SHet_p"].astype(float) < 5e-8) | (cpassoc["SHom_p"].astype(float) < 5e-8)).sum())
        ),
    }
)

ordered = [
    "pair_total",
    "LDSC_positive_pairs",
    "LAVA_positive_pairs",
    "PLACO_SNPs",
    "CPASSOC_SNPs",
    "coloc_loci",
    "SuSiE_loci",
    "MAGMA_significant_genes",
    "GO_significant_pathways",
    "KEGG_significant_pathways",
    "Reactome_significant_pathways",
    "PLACO_genomewide_pairs",
    "CPASSOC_genomewide_pairs",
    "PLACO_genomewide_significant_rows_total_P5e8",
    "PLACO_genomewide_P_lt_5e8_SNPs",
    "CPASSOC_genomewide_significant_rows_total_P5e8",
    "CPASSOC_genomewide_any_P_lt_5e8_SNPs",
]

pd.DataFrame([{"metric": k, "value": metrics[k]} for k in ordered if k in metrics]).to_csv(
    gqc_path, sep="\t", index=False
)

with (root / "update_log.md").open("a") as fh:
    fh.write("\n## Threshold QC metrics\n\n")
    fh.write("- Restored positive-result QC for `P<5e-8`; BH-FDR is not used as the PLACO/CPASSOC positive-result threshold.\n")

print("RESTORED_P5E8_QC")
