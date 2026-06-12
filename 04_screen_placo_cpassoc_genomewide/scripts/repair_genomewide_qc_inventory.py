from pathlib import Path
from datetime import datetime
import os

import pandas as pd


project_root = os.environ.get("PROJECT_ROOT", "").strip()
if not project_root:
    raise SystemExit("Set PROJECT_ROOT to the analysis project root.")
root = Path(project_root) / "results"

placo = pd.read_csv(root / "placo_genomewide/placo_genomewide_significant_snps.tsv", sep="\t")
cpassoc = pd.read_csv(root / "cpassoc_genomewide/cpassoc_genomewide_significant_snps.tsv", sep="\t")
placo_qc = pd.read_csv(root / "placo_genomewide/placo_qc.tsv", sep="\t")
cpassoc_qc = pd.read_csv(root / "cpassoc_genomewide/cpassoc_qc.tsv", sep="\t")

gqc_path = root / "qc/global_qc_report.tsv"
if gqc_path.exists():
    gqc = pd.read_csv(gqc_path, sep="\t")
    metrics = dict(zip(gqc["metric"].astype(str), gqc["value"].astype(str)))
else:
    metrics = {}

metrics.update(
    {
        "PLACO_genomewide_significant_rows_total_P5e8": str(len(placo)),
        "CPASSOC_genomewide_significant_rows_total_P5e8": str(len(cpassoc)),
        "PLACO_genomewide_P_lt_5e8_SNPs": str(int((placo["PLACO_p"].astype(float) < 5e-8).sum())),
        "CPASSOC_genomewide_SHet_P_lt_5e8_SNPs": str(int((cpassoc["SHet_p"].astype(float) < 5e-8).sum())),
        "CPASSOC_genomewide_SHom_P_lt_5e8_SNPs": str(int((cpassoc["SHom_p"].astype(float) < 5e-8).sum())),
        "CPASSOC_genomewide_any_P_lt_5e8_SNPs": str(
            int(((cpassoc["SHet_p"].astype(float) < 5e-8) | (cpassoc["SHom_p"].astype(float) < 5e-8)).sum())
        ),
        "PLACO_genomewide_pairs": str(int(placo_qc.shape[0])),
        "CPASSOC_genomewide_pairs": str(int(cpassoc_qc.shape[0])),
    }
)
pd.DataFrame([{"metric": k, "value": v} for k, v in metrics.items()]).to_csv(gqc_path, sep="\t", index=False)

desc = {
    "placo_genomewide": "Genome-wide PLACO significant results, QC, logs, and benchmark files",
    "cpassoc_genomewide": "Genome-wide CPASSOC significant results, QC, logs, and benchmark files",
    "integrated": "Integrated pair/SNP/locus summaries",
    "qc": "Global and method-level QC reports",
}

rows = []
for p in sorted(root.rglob("*")):
    if not p.is_file():
        continue
    n_rows = ""
    n_cols = ""
    if p.suffix == ".tsv":
        try:
            with p.open("rb") as fh:
                line_count = sum(1 for _ in fh)
            n_rows = max(line_count - 1, 0)
            with p.open("r", errors="replace") as fh:
                header = fh.readline().rstrip("\n")
            n_cols = len(header.split("\t")) if header else 0
        except Exception:
            n_rows = "NA"
            n_cols = "NA"
    rel_parts = p.relative_to(root).parts
    part = rel_parts[0] if rel_parts else ""
    rows.append(
        {
            "file_path": p.as_posix(),
            "file_type": p.suffix.lstrip(".") or "no_extension",
            "n_rows": n_rows,
            "n_cols": n_cols,
            "file_size": p.stat().st_size,
            "description": desc.get(part, part or "results file"),
            "last_updated_step": (
                "genomewide_placo_cpassoc_sensitivity"
                if part in {"placo_genomewide", "cpassoc_genomewide", "integrated", "qc"}
                else "pre_existing_or_previous_step"
            ),
        }
    )
pd.DataFrame(rows).to_csv(root / "final_project_inventory.tsv", sep="\t", index=False)

with (root / "update_log.md").open("a") as fh:
    fh.write("\n## QC metric clarification\n\n")
    fh.write(
        f"- Updated global QC at {datetime.now().isoformat(timespec='seconds')} to report total retained "
        "genome-wide PLACO/CPASSOC rows using P<5e-8 as the positive-result threshold.\n"
    )
    fh.write("- Expanded `final_project_inventory.tsv` to list all files under `results/`.\n")

print("QC_AND_INVENTORY_REPAIRED")
