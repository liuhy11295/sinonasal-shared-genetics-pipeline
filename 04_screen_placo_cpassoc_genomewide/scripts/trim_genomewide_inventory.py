from pathlib import Path

import pandas as pd


root = Path("/platform_data/p_user/p010/phase0/results")

core_files = [
    ("placo_genomewide/placo_genomewide_significant_snps.tsv", "positive PLACO SNPs retained by P<5e-8"),
    ("placo_genomewide/placo_genomewide_summary.tsv", "per-pair genome-wide PLACO positive SNP counts"),
    ("placo_genomewide/placo_qc.tsv", "per-pair genome-wide PLACO QC"),
    ("placo_genomewide/README.md", "genome-wide PLACO result notes"),
    ("cpassoc_genomewide/cpassoc_genomewide_significant_snps.tsv", "positive CPASSOC SNPs retained by P<5e-8"),
    ("cpassoc_genomewide/cpassoc_genomewide_summary.tsv", "per-pair genome-wide CPASSOC positive SNP counts"),
    ("cpassoc_genomewide/cpassoc_qc.tsv", "per-pair genome-wide CPASSOC QC"),
    ("cpassoc_genomewide/README.md", "genome-wide CPASSOC result notes"),
    ("integrated/snp_master_table.tsv", "integrated positive genome-wide PLACO/CPASSOC SNP table"),
    ("integrated/locus_master_table.tsv", "integrated positive SNP-position locus table"),
    ("integrated/pair_snp_locus_summary.tsv", "pair-level genome-wide PLACO/CPASSOC support summary"),
    ("integrated/pair_master_summary.tsv", "updated pair master summary with genome-wide support fields"),
    ("qc/global_qc_report.tsv", "global positive-result QC summary"),
    ("README.md", "project-level result notes"),
    ("method_notes.md", "method notes distinguishing candidate-based and genome-wide analyses"),
    ("update_log.md", "update log for genome-wide PLACO/CPASSOC sensitivity"),
]


def table_shape(path: Path):
    if path.suffix != ".tsv" or not path.exists():
        return "", ""
    with path.open("rb") as fh:
        n_rows = max(sum(1 for _ in fh) - 1, 0)
    with path.open("r", errors="replace") as fh:
        header = fh.readline().rstrip("\n")
    n_cols = len(header.split("\t")) if header else 0
    return n_rows, n_cols


rows = []
for rel, description in core_files:
    p = root / rel
    n_rows, n_cols = table_shape(p)
    rows.append(
        {
            "file_path": p.as_posix(),
            "file_type": p.suffix.lstrip(".") or "no_extension",
            "n_rows": n_rows,
            "n_cols": n_cols,
            "file_size": p.stat().st_size if p.exists() else "",
            "description": description,
            "last_updated_step": "genomewide_placo_cpassoc_sensitivity",
        }
    )

pd.DataFrame(rows).to_csv(root / "final_project_inventory.tsv", sep="\t", index=False)

gqc_path = root / "qc/global_qc_report.tsv"
gqc = pd.read_csv(gqc_path, sep="\t")
keep = [
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
    "PLACO_genomewide_significant_rows_total_P5e8",
    "CPASSOC_genomewide_significant_rows_total_P5e8",
    "PLACO_genomewide_pairs",
    "CPASSOC_genomewide_pairs",
]
gqc = gqc[gqc["metric"].isin(keep)]
gqc.to_csv(gqc_path, sep="\t", index=False)

with (root / "update_log.md").open("a") as fh:
    fh.write("\n## Inventory trim\n\n")
    fh.write("- Trimmed `final_project_inventory.tsv` back to positive-result core files only.\n")
    fh.write("- Kept global QC focused on retained positive rows and pair counts, without expanded threshold-specific inventory.\n")

print("TRIMMED_TO_POSITIVE_CORE_FILES")
