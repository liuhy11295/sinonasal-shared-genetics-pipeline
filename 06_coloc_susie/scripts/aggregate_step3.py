#!/usr/bin/env python3
import csv
import os
from collections import OrderedDict
from pathlib import Path

project_root = os.environ.get("PROJECT_ROOT", "")
if not project_root:
    raise SystemExit("Set PROJECT_ROOT")
PROJECT_ROOT = Path(project_root).expanduser().resolve()
STEP2 = PROJECT_ROOT / "results/phase0_extension/step2_input_tables"
STEP3 = PROJECT_ROOT / "results/phase0_extension/step3_susie_coloc"
PER = STEP3 / "per_locus"


def read_tsv(path):
    if not path.exists() or path.stat().st_size == 0:
        return []
    with path.open(newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def write_tsv(path, rows, fields):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, delimiter="\t", fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fields})


def fnum(x):
    try:
        return float(x)
    except Exception:
        return None


susie_fields = [
    "pair_id", "trait1", "trait2", "locus_id", "CHR", "START", "END", "trait",
    "credible_set_id", "SNP", "PIP", "is_cs_lead",
]
coloc_susie_fields = [
    "pair_id", "trait1", "trait2", "locus_id", "shared_signal_id", "SNP",
    "PIP_trait1", "PIP_trait2", "PP.H4",
]
high_conf_fields = [
    "pair_id", "trait1", "trait2", "locus_id", "SNP", "CHR", "BP",
    "PIP_trait1", "PIP_trait2", "PP.H4",
]
qc_fields = [
    "pair_id", "trait1", "trait2", "locus_id", "status", "reason",
    "n_snps", "n_ref_snps", "n_shared_snps", "n_cs_trait1", "n_cs_trait2",
]
candidate_fields = [
    "pair_id", "trait1", "trait2", "locus_id", "SNP", "CHR", "BP", "EA", "OA",
    "BETA_trait1", "SE_trait1", "P_trait1", "BETA_trait2", "SE_trait2", "P_trait2",
    "source",
]

susie_rows = []
coloc_rows = []
qc_rows = []
for path in sorted(PER.glob("*.susie_finemap.tsv")):
    susie_rows.extend(read_tsv(path))
for path in sorted(PER.glob("*.coloc_susie.tsv")):
    coloc_rows.extend(read_tsv(path))
for path in sorted(PER.glob("*.qc.tsv")):
    qc_rows.extend(read_tsv(path))

write_tsv(STEP3 / "susie_finemap_results.tsv", susie_rows, susie_fields)
write_tsv(STEP3 / "coloc_susie_results.tsv", coloc_rows, coloc_susie_fields)
write_tsv(STEP3 / "susie_qc_report.tsv", qc_rows, qc_fields)

candidate_rows = read_tsv(STEP2 / "candidate_snps.tsv")
candidate_by_key = {
    (r.get("pair_id"), r.get("locus_id"), r.get("SNP")): r for r in candidate_rows
}
candidate_out = list(candidate_rows)
seen_full = {
    (r.get("pair_id"), r.get("locus_id"), r.get("SNP"), r.get("source")) for r in candidate_out
}
new_candidate_count = 0

susie_lookup = OrderedDict()
for r in susie_rows:
    pip = fnum(r.get("PIP"))
    if pip is not None and pip >= 0.1:
        key = (r.get("pair_id"), r.get("locus_id"), r.get("SNP"))
        if key not in susie_lookup:
            susie_lookup[key] = r

for key, r in susie_lookup.items():
    pair_id, locus_id, snp = key
    base = candidate_by_key.get(key, {})
    out = {
        "pair_id": pair_id,
        "trait1": r.get("trait1"),
        "trait2": r.get("trait2"),
        "locus_id": locus_id,
        "SNP": snp,
        "CHR": r.get("CHR"),
        "BP": base.get("BP", ""),
        "EA": base.get("EA", ""),
        "OA": base.get("OA", ""),
        "BETA_trait1": base.get("BETA_trait1", ""),
        "SE_trait1": base.get("SE_trait1", ""),
        "P_trait1": base.get("P_trait1", ""),
        "BETA_trait2": base.get("BETA_trait2", ""),
        "SE_trait2": base.get("SE_trait2", ""),
        "P_trait2": base.get("P_trait2", ""),
        "source": "SuSiE_high_PIP",
    }
    full = (pair_id, locus_id, snp, out["source"])
    if full not in seen_full:
        seen_full.add(full)
        candidate_out.append(out)
        new_candidate_count += 1

high_conf = []
for r in coloc_rows:
    pp4 = fnum(r.get("PP.H4"))
    pp3 = fnum(r.get("PP.H3"))
    pip1 = fnum(r.get("PIP_trait1"))
    pip2 = fnum(r.get("PIP_trait2"))
    if pp4 is None or pip1 is None or pip2 is None:
        continue
    if pp4 >= 0.8 and (pp3 is None or pp4 > pp3) and pip1 >= 0.1 and pip2 >= 0.1:
        key = (r.get("pair_id"), r.get("locus_id"), r.get("SNP"))
        base = candidate_by_key.get(key, {})
        susie_base = susie_lookup.get(key, {})
        high_conf.append({
            "pair_id": r.get("pair_id"),
            "trait1": r.get("trait1"),
            "trait2": r.get("trait2"),
            "locus_id": r.get("locus_id"),
            "SNP": r.get("SNP"),
            "CHR": susie_base.get("CHR", base.get("CHR", "")),
            "BP": base.get("BP", ""),
            "PIP_trait1": r.get("PIP_trait1"),
            "PIP_trait2": r.get("PIP_trait2"),
            "PP.H4": r.get("PP.H4"),
        })
        out = {
            "pair_id": r.get("pair_id"),
            "trait1": r.get("trait1"),
            "trait2": r.get("trait2"),
            "locus_id": r.get("locus_id"),
            "SNP": r.get("SNP"),
            "CHR": susie_base.get("CHR", base.get("CHR", "")),
            "BP": base.get("BP", ""),
            "EA": base.get("EA", ""),
            "OA": base.get("OA", ""),
            "BETA_trait1": base.get("BETA_trait1", ""),
            "SE_trait1": base.get("SE_trait1", ""),
            "P_trait1": base.get("P_trait1", ""),
            "BETA_trait2": base.get("BETA_trait2", ""),
            "SE_trait2": base.get("SE_trait2", ""),
            "P_trait2": base.get("P_trait2", ""),
            "source": "coloc_susie_shared_signal",
        }
        full = (out["pair_id"], out["locus_id"], out["SNP"], out["source"])
        if full not in seen_full:
            seen_full.add(full)
            candidate_out.append(out)
            new_candidate_count += 1

write_tsv(STEP3 / "high_confidence_shared_variants.tsv", high_conf, high_conf_fields)
write_tsv(STEP3 / "candidate_snps.tsv", candidate_out, candidate_fields)

summary_fields = [
    "analyzed_loci", "success_loci", "failed_loci", "credible_set_total",
    "coloc_susie_shared_signal_count", "high_confidence_shared_variants",
    "new_candidate_snps",
]
success_loci = sum(1 for r in qc_rows if r.get("status") == "success")
failed_loci = sum(1 for r in qc_rows if r.get("status") != "success")
credible_set_total = 0
for r in qc_rows:
    for key in ("n_cs_trait1", "n_cs_trait2"):
        try:
            credible_set_total += int(float(r.get(key) or 0))
        except Exception:
            pass
shared_signal_count = len({(r.get("pair_id"), r.get("locus_id"), r.get("shared_signal_id")) for r in coloc_rows})
summary = [{
    "analyzed_loci": len(qc_rows),
    "success_loci": success_loci,
    "failed_loci": failed_loci,
    "credible_set_total": credible_set_total,
    "coloc_susie_shared_signal_count": shared_signal_count,
    "high_confidence_shared_variants": len(high_conf),
    "new_candidate_snps": new_candidate_count,
}]
write_tsv(STEP3 / "step3_summary.tsv", summary, summary_fields)
print("\t".join(summary_fields))
print("\t".join(str(summary[0][k]) for k in summary_fields))
