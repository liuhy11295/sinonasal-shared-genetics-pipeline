#!/usr/bin/env python3
import csv
import os
from collections import Counter, defaultdict
from pathlib import Path

project_root = os.environ.get("PROJECT_ROOT", "")
if not project_root:
    raise SystemExit("Set PROJECT_ROOT")
PROJECT_ROOT = Path(project_root).expanduser().resolve()
STEP2 = PROJECT_ROOT / "results/phase0_extension/step2_input_tables"
STEP3 = PROJECT_ROOT / "results/phase0_extension/step3_susie_coloc"
STEP4 = PROJECT_ROOT / "results/phase0_extension/step4_lcv_mr"
PER = STEP4 / "per_pair"


def read_tsv(path):
    if not path.exists() or path.stat().st_size == 0:
        return []
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


feas = []
mr = []
sens = []
for p in sorted(PER.glob("*.mr_feasibility.tsv")):
    feas.extend(read_tsv(p))
for p in sorted(PER.glob("*.mr_results.tsv")):
    mr.extend(read_tsv(p))
for p in sorted(PER.glob("*.mr_sensitivity.tsv")):
    sens.extend(read_tsv(p))

write_tsv(STEP4 / "mr_feasibility_report.tsv", feas, ["pair_id", "trait1", "trait2", "n_IV", "status", "reason"])
write_tsv(STEP4 / "mr_results.tsv", mr, ["pair_id", "exposure", "outcome", "method", "beta", "se", "p_value", "OR", "CI_lower", "CI_upper"])
write_tsv(STEP4 / "mr_sensitivity.tsv", sens, ["pair_id", "exposure", "outcome", "egger_intercept", "egger_p", "heterogeneity_q", "heterogeneity_p", "mr_presso_global_p", "n_outlier"])

priority = {r["pair_id"]: r for r in read_tsv(STEP2 / "priority_pairs.tsv")}
lcv = {r["pair_id"]: r for r in read_tsv(STEP4 / "lcv_results.tsv")}
candidates = {r["pair_id"]: r for r in read_tsv(STEP4 / "mr_candidate_pairs.tsv")}
high_conf_pairs = {r["pair_id"] for r in read_tsv(STEP3 / "high_confidence_shared_variants.tsv")}
mr_sig = defaultdict(list)
for r in mr:
    p = fnum(r.get("p_value"))
    if p is not None and p < 0.05:
        mr_sig[r["pair_id"]].append(r)

summary = []
for pid, cand in candidates.items():
    lcv_r = lcv.get(pid, {})
    lcv_p = fnum(lcv_r.get("p_value"))
    lcv_support = "yes" if lcv_p is not None and lcv_p < 0.05 else "no"
    mr_support = "yes" if pid in mr_sig else "no"
    has_high_conf = pid in high_conf_pairs
    if lcv_support == "yes" and mr_support == "yes" and has_high_conf:
        cls = "strong_causal_evidence"
    elif (lcv_support == "yes" and mr_support == "yes") or ((lcv_support == "yes" or mr_support == "yes") and has_high_conf):
        cls = "moderate_causal_evidence"
    elif lcv_support == "yes" or mr_support == "yes":
        cls = "weak_causal_evidence"
    else:
        cls = "shared_genetics_only"
    notes = []
    if lcv_r.get("interpretation"):
        notes.append("LCV=" + lcv_r.get("interpretation", ""))
    if pid in mr_sig:
        notes.append("MR_significant_methods=" + str(len(mr_sig[pid])))
    if has_high_conf:
        notes.append("coloc_susie_high_confidence")
    summary.append({
        "pair_id": pid,
        "trait1": cand["trait1"],
        "trait2": cand["trait2"],
        "LCV_support": lcv_support,
        "MR_support": mr_support,
        "final_classification": cls,
        "note": ";".join(notes),
    })

write_tsv(STEP4 / "causal_evidence_summary.tsv", summary, ["pair_id", "trait1", "trait2", "LCV_support", "MR_support", "final_classification", "note"])

counts = Counter(r["final_classification"] for r in summary)
lcv_sig_n = sum(1 for r in lcv.values() if fnum(r.get("p_value")) is not None and fnum(r["p_value"]) < 0.05)
eligible_n = sum(1 for r in feas if r.get("status") == "eligible")
mr_sig_n = len(mr_sig)
out = [{
    "lcv_analyzed_pairs": len(lcv),
    "lcv_significant_pairs": lcv_sig_n,
    "mr_candidate_pairs": len(candidates),
    "mr_eligible_pairs": eligible_n,
    "mr_significant_pairs": mr_sig_n,
    "strong_causal_evidence": counts["strong_causal_evidence"],
    "moderate_causal_evidence": counts["moderate_causal_evidence"],
    "shared_genetics_only": counts["shared_genetics_only"],
}]
fields = ["lcv_analyzed_pairs", "lcv_significant_pairs", "mr_candidate_pairs", "mr_eligible_pairs", "mr_significant_pairs", "strong_causal_evidence", "moderate_causal_evidence", "shared_genetics_only"]
write_tsv(STEP4 / "step4_summary.tsv", out, fields)
print("\t".join(fields))
print("\t".join(str(out[0][field]) for field in fields))
