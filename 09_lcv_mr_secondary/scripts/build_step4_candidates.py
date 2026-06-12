#!/usr/bin/env python3
import csv
import os
from collections import defaultdict
from pathlib import Path

project_root = os.environ.get("PROJECT_ROOT", "").strip()
if not project_root:
    raise SystemExit("Set PROJECT_ROOT to the analysis project root.")
PROJECT_ROOT = Path(project_root)
STEP2 = PROJECT_ROOT / "results/phase0_extension/step2_input_tables"
STEP3 = PROJECT_ROOT / "results/phase0_extension/step3_susie_coloc"
STEP4 = PROJECT_ROOT / "results/phase0_extension/step4_lcv_mr"


def read_tsv(path):
    with path.open(newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def write_tsv(path, rows, fields):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, delimiter="\t", fieldnames=fields, extrasaction="ignore")
        w.writeheader()
        for row in rows:
            w.writerow({field: row.get(field, "") for field in fields})


def fnum(x):
    try:
        return float(x)
    except Exception:
        return None


priority = read_tsv(STEP2 / "priority_pairs.tsv")
lcv = {r["pair_id"]: r for r in read_tsv(STEP4 / "lcv_results.tsv")}
coloc = read_tsv(STEP2 / "coloc_positive_loci.tsv")
coloc_susie = read_tsv(STEP3 / "coloc_susie_results.tsv")

lava_counts = {r["pair_id"]: int(float(r.get("n_LAVA_loci") or 0)) for r in priority}
ldsc_sig = {
    r["pair_id"]
    for r in priority
    if "LDSC_significant" in (r.get("priority_reason") or "").split(";")
}
coloc_h4 = {r["pair_id"] for r in coloc if fnum(r.get("PP.H4")) is not None and fnum(r["PP.H4"]) >= 0.5}
susie_shared = {
    r["pair_id"]
    for r in coloc_susie
    if fnum(r.get("PP.H4")) is not None and fnum(r["PP.H4"]) >= 0.8
}
lcv_sig = {
    pid for pid, r in lcv.items()
    if fnum(r.get("p_value")) is not None and fnum(r["p_value"]) < 0.05
}

rows = []
for r in priority:
    pid = r["pair_id"]
    reasons = []
    if pid in lcv_sig:
        reasons.append("LCV_significant")
    if pid in ldsc_sig and lava_counts.get(pid, 0) > 0:
        reasons.append("LDSC_significant_and_LAVA_positive")
    if pid in coloc_h4:
        reasons.append("coloc_H4_positive")
    if pid in susie_shared:
        reasons.append("coloc_susie_shared_signal")
    if not reasons:
        continue
    rows.append({
        "pair_id": pid,
        "trait1": r["trait1"],
        "trait2": r["trait2"],
        "selection_reason": ";".join(reasons),
    })

write_tsv(STEP4 / "mr_candidate_pairs.tsv", rows, ["pair_id", "trait1", "trait2", "selection_reason"])
write_tsv(STEP4 / "mr_candidate_summary.tsv", [{"n_candidate_pairs": len(rows)}], ["n_candidate_pairs"])
print(f"n_candidate_pairs\t{len(rows)}")
