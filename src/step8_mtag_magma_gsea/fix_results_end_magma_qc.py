from pathlib import Path
import csv
import math
from collections import Counter, defaultdict


ROOT = Path.cwd() / "result_5_31" / "results-end"
MAGMA = ROOT / "magma" / "gene_results.tsv"
OUTDIR = ROOT / "magma"
QC = ROOT / "qc" / "global_qc_report.tsv"
INV = ROOT / "final_project_inventory.tsv"


def to_float(value):
    try:
        if value is None or value == "":
            return math.nan
        return float(value)
    except Exception:
        return math.nan


def is_num(value):
    return not math.isnan(value)


def write_dicts(path, columns, data):
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=columns, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(data)


rows = []
pvals_by_pair = defaultdict(list)
with MAGMA.open("r", encoding="utf-8", errors="ignore", newline="") as handle:
    reader = csv.DictReader(handle, delimiter="\t")
    fields = reader.fieldnames or []
    pcol = "p" if "p" in fields else "P"
    gene_key = "gene_id" if "gene_id" in fields else "gene_symbol"
    for index, row in enumerate(reader):
        pair = row.get("pair_id", "")
        pvalue = to_float(row.get(pcol, ""))
        bonf = to_float(row.get("gene_bonf_p", ""))
        if not is_num(bonf):
            tested = to_float(row.get("n_tested_genes", ""))
            bonf = pvalue * tested if is_num(pvalue) and is_num(tested) else math.nan
        record = {
            "pair_id": pair,
            "trait1": row.get("trait1", ""),
            "trait2": row.get("trait2", ""),
            "gene": row.get(gene_key, ""),
            "p": pvalue,
            "bonf": bonf,
            "fdr": math.nan,
        }
        rows.append(record)
        if is_num(pvalue):
            pvals_by_pair[pair].append((pvalue, index))

for pair, values in pvals_by_pair.items():
    values.sort(key=lambda item: item[0])
    count = len(values)
    qvalues = [
        min(pvalue * count / rank, 1.0)
        for rank, (pvalue, _index) in enumerate(values, start=1)
    ]
    for idx in range(count - 2, -1, -1):
        qvalues[idx] = min(qvalues[idx], qvalues[idx + 1])
    for (_pvalue, index), qvalue in zip(values, qvalues):
        rows[index]["fdr"] = qvalue

pair_meta = {}
pair_genes = defaultdict(set)
pair_bonf = Counter()
pair_fdr = Counter()
pair_nominal = Counter()
pair_min_p = defaultdict(lambda: math.inf)
pair_min_bonf = defaultdict(lambda: math.inf)
pair_min_fdr = defaultdict(lambda: math.inf)
bonf_genes = set()
fdr_genes = set()
nominal_genes = set()
bonf_rows = 0
fdr_rows = 0
nominal_rows = 0

for row in rows:
    pair = row["pair_id"]
    pair_meta[pair] = (row["trait1"], row["trait2"])
    if row["gene"]:
        pair_genes[pair].add(row["gene"])
    if is_num(row["p"]):
        pair_min_p[pair] = min(pair_min_p[pair], row["p"])
    if is_num(row["bonf"]):
        pair_min_bonf[pair] = min(pair_min_bonf[pair], row["bonf"])
    if is_num(row["fdr"]):
        pair_min_fdr[pair] = min(pair_min_fdr[pair], row["fdr"])
    if is_num(row["bonf"]) and row["bonf"] <= 0.05:
        pair_bonf[pair] += 1
        bonf_rows += 1
        if row["gene"]:
            bonf_genes.add(row["gene"])
    if is_num(row["fdr"]) and row["fdr"] < 0.05:
        pair_fdr[pair] += 1
        fdr_rows += 1
        if row["gene"]:
            fdr_genes.add(row["gene"])
    if is_num(row["p"]) and row["p"] < 0.05:
        pair_nominal[pair] += 1
        nominal_rows += 1
        if row["gene"]:
            nominal_genes.add(row["gene"])

summary_rows = []
for pair in sorted(pair_meta):
    trait1, trait2 = pair_meta[pair]
    summary_rows.append(
        {
            "pair_id": pair,
            "trait1": trait1,
            "trait2": trait2,
            "tested_genes": len(pair_genes[pair]),
            "bonferroni_sig_genes": pair_bonf[pair],
            "fdr05_genes": pair_fdr[pair],
            "nominal_p05_genes": pair_nominal[pair],
            "min_p": pair_min_p[pair] if pair_min_p[pair] != math.inf else "",
            "min_bonf_p": (
                pair_min_bonf[pair] if pair_min_bonf[pair] != math.inf else ""
            ),
            "min_fdr": pair_min_fdr[pair] if pair_min_fdr[pair] != math.inf else "",
        }
    )
summary_rows.sort(
    key=lambda row: (
        -row["bonferroni_sig_genes"],
        -row["fdr05_genes"],
        -row["nominal_p05_genes"],
        row["pair_id"],
    )
)

write_dicts(
    OUTDIR / "pair_gene_count_summary.tsv",
    [
        "pair_id",
        "trait1",
        "trait2",
        "tested_genes",
        "bonferroni_sig_genes",
        "fdr05_genes",
        "nominal_p05_genes",
        "min_p",
        "min_bonf_p",
        "min_fdr",
    ],
    summary_rows,
)


def top_pair(counter):
    return counter.most_common(1)[0] if counter else ("", 0)


top_bonf, n_top_bonf = top_pair(pair_bonf)
top_fdr, n_top_fdr = top_pair(pair_fdr)
top_nominal, n_top_nominal = top_pair(pair_nominal)
write_dicts(
    OUTDIR / "magma_threshold_comparison.tsv",
    [
        "threshold",
        "n_gene_rows",
        "n_unique_genes",
        "n_pairs_with_any_gene",
        "top_pair_by_count",
        "top_pair_gene_count",
    ],
    [
        {
            "threshold": "Bonferroni<=0.05",
            "n_gene_rows": bonf_rows,
            "n_unique_genes": len(bonf_genes),
            "n_pairs_with_any_gene": sum(1 for value in pair_bonf.values() if value > 0),
            "top_pair_by_count": top_bonf,
            "top_pair_gene_count": n_top_bonf,
        },
        {
            "threshold": "FDR<0.05_pairwise_BH",
            "n_gene_rows": fdr_rows,
            "n_unique_genes": len(fdr_genes),
            "n_pairs_with_any_gene": sum(1 for value in pair_fdr.values() if value > 0),
            "top_pair_by_count": top_fdr,
            "top_pair_gene_count": n_top_fdr,
        },
        {
            "threshold": "Nominal_P<0.05",
            "n_gene_rows": nominal_rows,
            "n_unique_genes": len(nominal_genes),
            "n_pairs_with_any_gene": sum(
                1 for value in pair_nominal.values() if value > 0
            ),
            "top_pair_by_count": top_nominal,
            "top_pair_gene_count": n_top_nominal,
        },
    ],
)

with QC.open("r", encoding="utf-8", errors="ignore", newline="") as handle:
    qc_rows = list(csv.DictReader(handle, delimiter="\t"))
remove_metrics = {
    "MAGMA_significant_genes",
    "MAGMA_significant_gene_rows_bonferroni",
    "MAGMA_significant_unique_genes_bonferroni",
    "MAGMA_pairs_with_bonferroni_genes",
    "MAGMA_significant_gene_rows_fdr05",
    "MAGMA_significant_unique_genes_fdr05",
    "MAGMA_nominal_p05_gene_rows",
    "MAGMA_nominal_p05_unique_genes",
}
qc_rows = [row for row in qc_rows if row.get("metric") not in remove_metrics]
qc_rows.extend(
    [
        {
            "metric": "MAGMA_significant_gene_rows_bonferroni",
            "value": str(bonf_rows),
        },
        {
            "metric": "MAGMA_significant_unique_genes_bonferroni",
            "value": str(len(bonf_genes)),
        },
        {
            "metric": "MAGMA_pairs_with_bonferroni_genes",
            "value": str(sum(1 for value in pair_bonf.values() if value > 0)),
        },
        {"metric": "MAGMA_significant_gene_rows_fdr05", "value": str(fdr_rows)},
        {
            "metric": "MAGMA_significant_unique_genes_fdr05",
            "value": str(len(fdr_genes)),
        },
        {"metric": "MAGMA_nominal_p05_gene_rows", "value": str(nominal_rows)},
        {
            "metric": "MAGMA_nominal_p05_unique_genes",
            "value": str(len(nominal_genes)),
        },
    ]
)
write_dicts(QC, ["metric", "value"], qc_rows)


def count_data_lines(path):
    with path.open("r", encoding="utf-8", errors="ignore") as handle:
        return max(sum(1 for _line in handle) - 1, 0)


with INV.open("r", encoding="utf-8", errors="ignore", newline="") as handle:
    inv_rows = list(csv.DictReader(handle, delimiter="\t"))
updated_paths = {
    "magma/pair_gene_count_summary.tsv",
    "magma/magma_threshold_comparison.tsv",
    "qc/global_qc_report.tsv",
}
inv_rows = [row for row in inv_rows if row.get("file_path") not in updated_paths]
for rel_path, description in [
    (
        "magma/pair_gene_count_summary.tsv",
        "Per-pair MAGMA tested genes and threshold counts.",
    ),
    (
        "magma/magma_threshold_comparison.tsv",
        "MAGMA Bonferroni, FDR, and nominal threshold comparison.",
    ),
    (
        "qc/global_qc_report.tsv",
        "Global QC report with corrected MAGMA gene-count metrics.",
    ),
]:
    path = ROOT / rel_path
    inv_rows.append(
        {
            "file_path": rel_path,
            "file_type": path.suffix.lstrip("."),
            "rows": str(count_data_lines(path)),
            "size_bytes": str(path.stat().st_size),
            "description": description,
        }
    )
inv_rows.sort(key=lambda row: row.get("file_path", ""))
write_dicts(
    INV,
    ["file_path", "file_type", "rows", "size_bytes", "description"],
    inv_rows,
)

print("updated_qc", QC)
print("summary_rows", len(summary_rows))
print("bonf_rows", bonf_rows, "bonf_unique_genes", len(bonf_genes))
print("fdr_rows", fdr_rows, "fdr_unique_genes", len(fdr_genes))
print("nominal_rows", nominal_rows, "nominal_unique_genes", len(nominal_genes))
