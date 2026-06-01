#!/usr/bin/env python3
import csv
from pathlib import Path

PROJECT_ROOT = Path("/platform_data/p_user/p010/phase0")
STEP2 = PROJECT_ROOT / "results/phase0_extension/step2_input_tables"
STEP5 = PROJECT_ROOT / "results/phase0_extension/step5_placo_cpassoc"
PER = STEP5 / "per_pair"

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
        if x in ("", "NA", "nan", "None"):
            return None
        return float(x)
    except Exception:
        return None

placo_fields = ["pair_id","trait1","trait2","SNP","CHR","BP","EA","OA","Z_trait1","Z_trait2","P_trait1","P_trait2","PLACO_stat","PLACO_p"]
cpassoc_fields = ["pair_id","trait1","trait2","SNP","CHR","BP","EA","OA","P_trait1","P_trait2","SHet","SHet_p","SHom","SHom_p"]
qc_fields = ["pair_id","trait1","trait2","method","status","reason"]
sig_fields = ["pair_id","trait1","trait2","SNP","CHR","BP","EA","OA","method","p_value","source"]

placo = []
cpassoc = []
qc = []
for path in sorted(PER.glob("*.placo.tsv")):
    placo.extend(read_tsv(path))
for path in sorted(PER.glob("*.cpassoc.tsv")):
    cpassoc.extend(read_tsv(path))
for path in sorted(PER.glob("*.qc.tsv")):
    qc.extend(read_tsv(path))

write_tsv(STEP5 / "placo_results.tsv", placo, placo_fields)
write_tsv(STEP5 / "cpassoc_results.tsv", cpassoc, cpassoc_fields)
write_tsv(STEP5 / "cross_trait_qc_report.tsv", qc, qc_fields)

sig = []
seen = set()
for row in placo:
    p = fnum(row.get("PLACO_p"))
    if p is not None and p < 5e-8:
        key = (row["pair_id"], row["SNP"], "PLACO")
        if key not in seen:
            seen.add(key)
            sig.append({
                "pair_id": row["pair_id"], "trait1": row["trait1"], "trait2": row["trait2"],
                "SNP": row["SNP"], "CHR": row["CHR"], "BP": row["BP"], "EA": row["EA"], "OA": row["OA"],
                "method": "PLACO", "p_value": row["PLACO_p"], "source": "placo_results.tsv"
            })
for row in cpassoc:
    for method, col in (("CPASSOC_SHet", "SHet_p"), ("CPASSOC_SHom", "SHom_p")):
        p = fnum(row.get(col))
        if p is not None and p < 5e-8:
            key = (row["pair_id"], row["SNP"], method)
            if key not in seen:
                seen.add(key)
                sig.append({
                    "pair_id": row["pair_id"], "trait1": row["trait1"], "trait2": row["trait2"],
                    "SNP": row["SNP"], "CHR": row["CHR"], "BP": row["BP"], "EA": row["EA"], "OA": row["OA"],
                    "method": method, "p_value": row[col], "source": "cpassoc_results.tsv"
                })
write_tsv(STEP5 / "cross_trait_significant_snps.tsv", sig, sig_fields)

old_candidate = read_tsv(STEP2 / "candidate_snps.tsv")
candidate_fields = ["pair_id","trait1","trait2","locus_id","SNP","CHR","BP","EA","OA","BETA_trait1","SE_trait1","P_trait1","BETA_trait2","SE_trait2","P_trait2","source"]
existing = {
    (r.get("pair_id",""), r.get("SNP",""), r.get("source",""))
    for r in old_candidate
}
new_candidates = []
for row in sig:
    source = "PLACO_significant_SNP" if row["method"] == "PLACO" else "CPASSOC_significant_SNP"
    key = (row["pair_id"], row["SNP"], source)
    if key in existing:
        continue
    existing.add(key)
    new_candidates.append({
        "pair_id": row["pair_id"], "trait1": row["trait1"], "trait2": row["trait2"],
        "locus_id": "", "SNP": row["SNP"], "CHR": row["CHR"], "BP": row["BP"],
        "EA": row["EA"], "OA": row["OA"],
        "BETA_trait1": "", "SE_trait1": "", "P_trait1": "",
        "BETA_trait2": "", "SE_trait2": "", "P_trait2": "",
        "source": source,
    })
write_tsv(STEP5 / "candidate_snps.tsv", old_candidate + new_candidates, candidate_fields)

placo_pairs = len({r["pair_id"] for r in placo if r.get("SNP")})
cpassoc_pairs = len({r["pair_id"] for r in cpassoc if r.get("SNP")})
placo_sig = sum(1 for r in sig if r["method"] == "PLACO")
cpassoc_sig = sum(1 for r in sig if r["method"].startswith("CPASSOC"))
failed = [r for r in qc if r.get("status") != "PASS"]
summary = [{
    "placo_analyzed_pairs": placo_pairs,
    "cpassoc_analyzed_pairs": cpassoc_pairs,
    "placo_significant_snps": placo_sig,
    "cpassoc_significant_snps": cpassoc_sig,
    "cross_trait_significant_snps": len(sig),
    "new_candidate_snps": len(new_candidates),
    "failed_records": len(failed),
}]
write_tsv(STEP5 / "step5_summary.tsv", summary, list(summary[0].keys()))
print("step5_summary", summary[0])
