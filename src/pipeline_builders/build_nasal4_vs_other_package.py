#!/usr/bin/env python3
"""Build the formal nasal4-vs-other Slurm package.

The previous nasal4 package is an input/QC package for the four replacement
nasal inflammatory GWAS. This script builds the formal Phase 0 package where
each replacement nasal GWAS is paired against the selected non-nasal comparator
traits from the project-level nasal-vs-other design.
"""

from __future__ import annotations

import csv
import shutil
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]
PKG = PROJECT_ROOT / "results" / "phase0_server" / "nasal4_vs_other_package"
NASAL4_PKG = PROJECT_ROOT / "results" / "phase0_server" / "nasal4_package"
PAIR_DESIGN = PROJECT_ROOT / "results" / "phase0_shared_genetics" / "disease_tables" / "phase0_nasal_vs_other_pair_manifest.tsv"
OLD_TRAIT_MANIFEST = PROJECT_ROOT / "results" / "phase0_server" / "ldsc_package" / "phase0_v3_ldsc_trait_manifest.tsv"
OLD_LAVA_R = PROJECT_ROOT / "results" / "phase0_server" / "local_rg_package" / "run_phase0_v3_lava.R"
OLD_LAVA_VERIFY = PROJECT_ROOT / "results" / "phase0_server" / "local_rg_package" / "verify_local_rg_outputs.py"
OLD_COLOC_EXTRACT = PROJECT_ROOT / "results" / "phase0_server" / "coloc_package" / "extract_coloc_regions.py"
OLD_COLOC_R = PROJECT_ROOT / "results" / "phase0_server" / "coloc_package" / "run_phase0_v3_coloc.R"
OLD_COLOC_VERIFY = PROJECT_ROOT / "results" / "phase0_server" / "coloc_package" / "verify_coloc_outputs.py"
OLD_MTAG_FILTER = PROJECT_ROOT / "results" / "phase0_server" / "mtag_package" / "filter_mtag_candidates.py"


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, fields: list[str], rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def write_text(path: Path, text: str, executable: bool = False) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text.strip() + "\n", encoding="utf-8", newline="\n")
    if executable:
        path.chmod(0o755)


def copy_file(src: Path, dst: Path, executable: bool = False) -> None:
    dst.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(src, dst)
    if executable:
        dst.chmod(0o755)


def slug_from_trait_manifest(row: dict[str, str]) -> str:
    return f"{row['trait']}.{row['accession']}"


def build_manifests() -> tuple[int, int]:
    nasal = read_tsv(NASAL4_PKG / "nasal4_manifest.tsv")
    old_pairs = read_tsv(PAIR_DESIGN)
    old_traits = {row["accession"]: row for row in read_tsv(OLD_TRAIT_MANIFEST)}

    comparators: dict[str, dict[str, str]] = {}
    for row in old_pairs:
        trait_b = row["trait_b"]
        if trait_b not in comparators:
            comparators[trait_b] = {
                "trait_id": trait_b,
                "trait_label": row["trait_name_b"],
                "module": row["module_b"],
            }

    trait_rows: list[dict[str, object]] = []
    for row in nasal:
        trait_id = row["trait_id"]
        n = int(row["cases"]) + int(row["controls"])
        trait_rows.append(
            {
                "trait_id": trait_id,
                "trait_label": row["trait_label"],
                "role": "replacement_nasal_anchor",
                "module": "M0_sinonasal_anchor",
                "source": row["source"],
                "accession": row["accession"],
                "cases": row["cases"],
                "controls": row["controls"],
                "sample_size": n,
                "common_file": f"results/phase0_server/nasal4_package/common/{trait_id}.common.tsv.gz",
                "ldsc_sumstats_file": f"results/phase0_server/nasal4_package/sumstats/{trait_id}.sumstats.gz",
                "mtag_input_file": f"results/phase0_server/nasal4_vs_other_package/inputs/{trait_id}.mtag.tsv.gz",
            }
        )

    for trait_id, meta in comparators.items():
        old = old_traits[trait_id]
        trait_rows.append(
            {
                "trait_id": trait_id,
                "trait_label": meta["trait_label"],
                "role": "other_comparator",
                "module": meta["module"],
                "source": "FinnGen_R12_existing_phase0",
                "accession": trait_id,
                "cases": old["cases"],
                "controls": old["controls"],
                "sample_size": old["sample_size"],
                "common_file": old["munged_file"],
                "ldsc_sumstats_file": old["ldsc_sumstats_file"],
                "mtag_input_file": f"results/phase0_server/nasal4_vs_other_package/inputs/{slug_from_trait_manifest(old)}.mtag.tsv.gz",
            }
        )

    trait_by_id = {row["trait_id"]: row for row in trait_rows}
    pair_rows: list[dict[str, object]] = []
    for nrow in nasal:
        nasal_id = nrow["trait_id"]
        for trait_b, brow in comparators.items():
            a = trait_by_id[nasal_id]
            b = trait_by_id[trait_b]
            pair_id = f"{nasal_id}__{trait_b}"
            pair_rows.append(
                {
                    "pair_id": pair_id,
                    "trait_a": nasal_id,
                    "trait_b": trait_b,
                    "module_a": a["module"],
                    "module_b": b["module"],
                    "trait_name_a": a["trait_label"],
                    "trait_name_b": b["trait_label"],
                    "pair_mode": "nasal4_vs_other",
                    "source_a": a["source"],
                    "source_b": b["source"],
                    "sumstats_a": a["common_file"],
                    "sumstats_b": b["common_file"],
                    "ldsc_sumstats_a": a["ldsc_sumstats_file"],
                    "ldsc_sumstats_b": b["ldsc_sumstats_file"],
                    "mtag_input_a": a["mtag_input_file"],
                    "mtag_input_b": b["mtag_input_file"],
                    "mtag_out_prefix": f"results/phase0_server/nasal4_vs_other_package/mtag/{pair_id}/{pair_id}",
                    "mtag_out_dir": f"results/phase0_server/nasal4_vs_other_package/mtag/{pair_id}",
                    "trait_a_cases": a["cases"],
                    "trait_a_controls": a["controls"],
                    "trait_b_cases": b["cases"],
                    "trait_b_controls": b["controls"],
                }
            )

    trait_fields = [
        "trait_id",
        "trait_label",
        "role",
        "module",
        "source",
        "accession",
        "cases",
        "controls",
        "sample_size",
        "common_file",
        "ldsc_sumstats_file",
        "mtag_input_file",
    ]
    pair_fields = [
        "pair_id",
        "trait_a",
        "trait_b",
        "module_a",
        "module_b",
        "trait_name_a",
        "trait_name_b",
        "pair_mode",
        "source_a",
        "source_b",
        "sumstats_a",
        "sumstats_b",
        "ldsc_sumstats_a",
        "ldsc_sumstats_b",
        "mtag_input_a",
        "mtag_input_b",
        "mtag_out_prefix",
        "mtag_out_dir",
        "trait_a_cases",
        "trait_a_controls",
        "trait_b_cases",
        "trait_b_controls",
    ]
    write_tsv(PKG / "nasal4_vs_other_trait_manifest.tsv", trait_fields, trait_rows)
    write_tsv(PKG / "nasal4_vs_other_pair_manifest.tsv", pair_fields, pair_rows)
    write_tsv(PKG / "phase0_v3_local_rg_pair_manifest.tsv", pair_fields, pair_rows)
    write_tsv(
        PKG / "phase0_v3_lava_gwas_info.tsv",
        ["phenotype", "cases", "controls", "filename"],
        [
            {
                "phenotype": row["trait_id"],
                "cases": row["cases"],
                "controls": row["controls"],
                "filename": row["ldsc_sumstats_file"],
            }
            for row in trait_rows
        ],
    )
    return len(trait_rows), len(pair_rows)


def write_ldsc_scripts() -> None:
    write_text(
        PKG / "summarize_nasal4_vs_other_ldsc.py",
        r'''
#!/usr/bin/env python3
import csv
import math
import re
from pathlib import Path

PKG = Path("results/phase0_server/nasal4_vs_other_package")

def read_tsv(path):
    with open(path, encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))

def write_tsv(path, fields, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)

def read_text(path):
    return path.read_text(encoding="utf-8", errors="replace") if path.exists() else ""

def first(pattern, text):
    match = re.search(pattern, text, flags=re.IGNORECASE)
    return match.group(1) if match else "NA"

def finite(value):
    try:
        return math.isfinite(float(value))
    except (TypeError, ValueError):
        return False

def parse_h2(path):
    text = read_text(path)
    return {
        "status": "OK" if path.exists() and "Analysis finished" in text and finite(first(r"Total Observed scale h2:\s*([-+0-9.eE]+)", text)) else "CHECK_LOG",
        "h2": first(r"Total Observed scale h2:\s*([-+0-9.eE]+)", text),
        "h2_se": first(r"Total Observed scale h2:\s*[-+0-9.eE]+\s*\(([-+0-9.eE]+)\)", text),
        "lambda_gc": first(r"Lambda GC:\s*([-+0-9.eE]+)", text),
        "intercept": first(r"Intercept:\s*([-+0-9.eE]+)", text),
    }

def parse_rg(path):
    text = read_text(path)
    rg = first(r"Genetic Correlation:\s*([-+0-9.eE]+)", text)
    se = first(r"Genetic Correlation:\s*[-+0-9.eE]+\s*\(([-+0-9.eE]+)\)", text)
    p = first(r"P:\s*([-+0-9.eE]+)", text)
    note = ""
    if path.exists() and finite(rg) and finite(se) and finite(p):
        status = "OK"
    elif path.exists() and "Analysis finished" in text and ("FloatingPointError" in text or "low" in text or "nan" in text.lower()):
        status = "UNSTABLE_RG_NA"
        note = "LDSC completed but rg is NA/unstable; likely low h2 or weak signal."
    else:
        status = "CHECK_LOG"
    return {
        "status": status,
        "rg": rg,
        "se": se,
        "p": p,
        "h2_obs": first(r"Total Observed scale h2:\s*([-+0-9.eE]+)", text),
        "h2_obs_se": first(r"Total Observed scale h2:\s*[-+0-9.eE]+\s*\(([-+0-9.eE]+)\)", text),
        "h2_int": first(r"Intercept:\s*([-+0-9.eE]+)", text),
        "gcov_int": first(r"Total Observed scale gencov:\s*([-+0-9.eE]+)", text),
        "note": note,
    }

def main():
    traits = read_tsv(PKG / "nasal4_vs_other_trait_manifest.tsv")
    pairs = read_tsv(PKG / "nasal4_vs_other_pair_manifest.tsv")
    h2_rows = []
    for trait in traits:
        path = PKG / "h2" / f"{trait['trait_id']}.log"
        parsed = parse_h2(path)
        parsed.update({"trait_id": trait["trait_id"], "trait_label": trait["trait_label"], "role": trait["role"], "module": trait["module"], "log": str(path)})
        h2_rows.append(parsed)
    rg_rows = []
    for pair in pairs:
        path = PKG / "rg" / f"{pair['pair_id']}.log"
        parsed = parse_rg(path)
        parsed.update({"pair_id": pair["pair_id"], "trait_a": pair["trait_a"], "trait_b": pair["trait_b"], "trait_name_a": pair["trait_name_a"], "trait_name_b": pair["trait_name_b"], "log": str(path)})
        rg_rows.append(parsed)
    write_tsv(PKG / "nasal4_vs_other_ldsc_h2_summary.tsv", ["trait_id","trait_label","role","module","status","h2","h2_se","lambda_gc","intercept","log"], h2_rows)
    write_tsv(PKG / "nasal4_vs_other_ldsc_rg_summary.tsv", ["pair_id","trait_a","trait_b","trait_name_a","trait_name_b","status","rg","se","p","h2_obs","h2_obs_se","h2_int","gcov_int","log","note"], rg_rows)
    verify = [
        {"check": "h2", "total": len(h2_rows), "ok": sum(r["status"] == "OK" for r in h2_rows), "unstable_or_skipped": 0, "failures": sum(r["status"] == "CHECK_LOG" for r in h2_rows)},
        {"check": "rg", "total": len(rg_rows), "ok": sum(r["status"] == "OK" for r in rg_rows), "unstable_or_skipped": sum(r["status"] == "UNSTABLE_RG_NA" for r in rg_rows), "failures": sum(r["status"] == "CHECK_LOG" for r in rg_rows)},
    ]
    write_tsv(PKG / "nasal4_vs_other_ldsc_verify.tsv", ["check","total","ok","unstable_or_skipped","failures"], verify)
    failures = sum(int(row["failures"]) for row in verify)
    print(f"nasal4_vs_other_ldsc_summary h2={len(h2_rows)} rg={len(rg_rows)} failures={failures}")
    return 1 if failures else 0

if __name__ == "__main__":
    raise SystemExit(main())
''',
        executable=True,
    )
    write_text(
        PKG / "run_nasal4_vs_other_ldsc.sh",
        r'''
#!/usr/bin/env bash
set -euo pipefail

PROJECT_ROOT="${1:-${PROJECT_ROOT:-/platform_data/p_user/p010/phase0}}"
cd "${PROJECT_ROOT}"
source results/phase0_server/server_transfer_package/activate_conda_runtime.sh "${PROJECT_ROOT}"
source results/phase0_server/server_transfer_package/no_conda_server_paths.example.env

: "${LDSC_DIR:?Set LDSC_DIR}"
: "${LD_REF_PREFIX:?Set LD_REF_PREFIX}"
: "${WLD_REF_PREFIX:?Set WLD_REF_PREFIX}"

PKG_DIR="results/phase0_server/nasal4_vs_other_package"
LOG_DIR="${PKG_DIR}/logs"
mkdir -p "${LOG_DIR}" "${PKG_DIR}/h2" "${PKG_DIR}/rg"

JOBS="${PHASE0_LDSC_JOBS:-${SLURM_NTASKS:-40}}"
echo "nasal4_vs_other LDSC jobs=${JOBS}"

python3 - <<'PY'
import csv, pathlib, sys
pkg = pathlib.Path("results/phase0_server/nasal4_vs_other_package")
missing = []
for row in csv.DictReader((pkg / "nasal4_vs_other_trait_manifest.tsv").open(), delimiter="\t"):
    for field in ["ldsc_sumstats_file"]:
        p = pathlib.Path(row[field])
        if not p.exists() or p.stat().st_size == 0:
            missing.append(str(p))
if missing:
    print("Missing LDSC sumstats:", *missing, sep="\n", file=sys.stderr)
    raise SystemExit(1)
PY

run_h2() {
  local trait="$1"
  local sumstats="$2"
  local out="${PKG_DIR}/h2/${trait}"
  if [[ -s "${out}.log" ]] && grep -q "Analysis finished" "${out}.log"; then
    return 0
  fi
  python "${LDSC_DIR}/ldsc.py" --h2 "${sumstats}" --ref-ld-chr "${LD_REF_PREFIX}" --w-ld-chr "${WLD_REF_PREFIX}" --out "${out}"
}

run_rg() {
  local pair_id="$1"
  local sumstats_a="$2"
  local sumstats_b="$3"
  local out="${PKG_DIR}/rg/${pair_id}"
  if [[ -s "${out}.log" ]] && grep -q "Analysis finished" "${out}.log"; then
    return 0
  fi
  python "${LDSC_DIR}/ldsc.py" --rg "${sumstats_a},${sumstats_b}" --ref-ld-chr "${LD_REF_PREFIX}" --w-ld-chr "${WLD_REF_PREFIX}" --out "${out}"
}

limit_jobs() {
  while (( $(jobs -rp | wc -l) >= JOBS )); do wait -n || true; done
}

while IFS=$'\t' read -r trait label role module source accession cases controls sample common sumstats mtag; do
  limit_jobs
  run_h2 "${trait}" "${sumstats}" >"${LOG_DIR}/h2_${trait}.stdout" 2>"${LOG_DIR}/h2_${trait}.stderr" &
done < <(tail -n +2 "${PKG_DIR}/nasal4_vs_other_trait_manifest.tsv")
wait

while IFS=$'\t' read -r pair_id trait_a trait_b module_a module_b name_a name_b pair_mode source_a source_b common_a common_b ldsc_a ldsc_b rest; do
  limit_jobs
  run_rg "${pair_id}" "${ldsc_a}" "${ldsc_b}" >"${LOG_DIR}/rg_${pair_id}.stdout" 2>"${LOG_DIR}/rg_${pair_id}.stderr" &
done < <(tail -n +2 "${PKG_DIR}/nasal4_vs_other_pair_manifest.tsv")
wait

python3 "${PKG_DIR}/summarize_nasal4_vs_other_ldsc.py"
''',
        executable=True,
    )


def write_mtag_scripts() -> None:
    write_text(
        PKG / "prepare_nasal4_vs_other_mtag_inputs.py",
        r'''
#!/usr/bin/env python3
import csv
import gzip
import math
from pathlib import Path

PKG = Path("results/phase0_server/nasal4_vs_other_package")
OUT_FIELDS = ["SNP", "CHR", "BP", "A1", "A2", "Z", "P", "EAF", "N"]

def open_text(path):
    return gzip.open(path, "rt", encoding="utf-8", errors="replace", newline="") if str(path).endswith(".gz") else open(path, encoding="utf-8", errors="replace", newline="")

def to_float(value):
    try:
        return float(value)
    except (TypeError, ValueError):
        return math.nan

def convert(source, dest):
    dest = Path(dest)
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_name(dest.name + ".tmp")
    if tmp.exists():
        tmp.unlink()
    read_rows = written = skipped = 0
    try:
        with open_text(source) as src, gzip.open(tmp, "wt", encoding="utf-8", newline="") as out:
            reader = csv.DictReader(src, delimiter="\t")
            writer = csv.DictWriter(out, fieldnames=OUT_FIELDS, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            for row in reader:
                read_rows += 1
                z = to_float(row.get("Z"))
                if not math.isfinite(z):
                    beta = to_float(row.get("BETA"))
                    se = to_float(row.get("SE"))
                    z = beta / se if math.isfinite(beta) and math.isfinite(se) and se > 0 else math.nan
                p = to_float(row.get("P"))
                eaf = to_float(row.get("EAF", row.get("FRQ")))
                n = to_float(row.get("N"))
                if not all(math.isfinite(x) for x in [z, p, eaf, n]) or not (0 <= p <= 1) or not (0 < eaf < 1) or n <= 0:
                    skipped += 1
                    continue
                chrom = row.get("CHR") or row.get("chr")
                bp = row.get("BP") or row.get("pos") or row.get("position")
                if not row.get("SNP") or not chrom or not bp or not row.get("A1") or not row.get("A2"):
                    skipped += 1
                    continue
                writer.writerow({"SNP": row["SNP"], "CHR": chrom, "BP": bp, "A1": row["A1"], "A2": row["A2"], "Z": f"{z:.12g}", "P": row["P"], "EAF": f"{eaf:.12g}", "N": f"{n:.12g}"})
                written += 1
    except Exception:
        if tmp.exists():
            tmp.unlink()
        raise
    if written == 0:
        if tmp.exists():
            tmp.unlink()
        raise RuntimeError(f"{source} produced zero MTAG rows")
    tmp.replace(dest)
    return read_rows, written, skipped

def main():
    selected_pairs = []
    with (PKG / "nasal4_vs_other_mtag_pair_manifest.tsv").open(encoding="utf-8", newline="") as handle:
        selected_pairs = list(csv.DictReader(handle, delimiter="\t"))
    selected_traits = {row["trait_a"] for row in selected_pairs} | {row["trait_b"] for row in selected_pairs}
    if not selected_traits:
        raise RuntimeError("No selected MTAG traits found; run select_nasal4_vs_other_mtag_pairs.py first")

    rows = []
    seen = set()
    with (PKG / "nasal4_vs_other_trait_manifest.tsv").open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["trait_id"] not in selected_traits:
                continue
            source = row["common_file"]
            dest = row["mtag_input_file"]
            if dest in seen:
                continue
            seen.add(dest)
            dest_path = Path(dest)
            if dest_path.exists() and dest_path.stat().st_size > 0:
                rows.append({"trait_id": row["trait_id"], "source_file": source, "mtag_input_file": dest, "read_rows": "NA", "written_rows": "NA", "skipped_rows": "NA", "status": "reused"})
                continue
            read_rows, written, skipped = convert(source, dest)
            rows.append({"trait_id": row["trait_id"], "source_file": source, "mtag_input_file": dest, "read_rows": read_rows, "written_rows": written, "skipped_rows": skipped, "status": "created"})
    out = PKG / "nasal4_vs_other_mtag_input_qc.tsv"
    with out.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["trait_id","source_file","mtag_input_file","read_rows","written_rows","skipped_rows","status"], delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    print(f"nasal4_vs_other_mtag_inputs selected_pairs={len(selected_pairs)} selected_traits={len(selected_traits)} prepared_traits={len(rows)}")

if __name__ == "__main__":
    main()
''',
        executable=True,
    )
    write_text(
        PKG / "select_nasal4_vs_other_mtag_pairs.py",
        r'''
#!/usr/bin/env python3
import argparse
import csv
import math
from pathlib import Path

PKG = Path("results/phase0_server/nasal4_vs_other_package")

def read_tsv(path):
    with open(path, encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--fdr", type=float, default=0.05)
    args = parser.parse_args()
    pairs = {row["pair_id"]: row for row in read_tsv(PKG / "nasal4_vs_other_pair_manifest.tsv")}
    rg_rows = read_tsv(PKG / "nasal4_vs_other_ldsc_rg_summary.tsv")

    valid = []
    for row in rg_rows:
        if row.get("status") != "OK":
            continue
        try:
            p = float(row.get("p", "nan"))
        except ValueError:
            continue
        if math.isfinite(p) and row["pair_id"] in pairs:
            valid.append((p, row))

    valid.sort(key=lambda item: item[0])
    q_by_pair = {}
    running_q = 1.0
    for idx in range(len(valid) - 1, -1, -1):
        rank = idx + 1
        p, row = valid[idx]
        running_q = min(running_q, p * len(valid) / rank)
        q_by_pair[row["pair_id"]] = min(running_q, 1.0)

    selected = []
    for p, row in valid:
        q = q_by_pair[row["pair_id"]]
        if q <= args.fdr:
            out = dict(pairs[row["pair_id"]])
            out["rg"] = row.get("rg", "")
            out["rg_p"] = row.get("p", "")
            out["rg_bh_q"] = f"{q:.12g}"
            out["rg_screen_method"] = "BH_FDR_on_LDSC_OK_pairs"
            out["rg_screen_threshold"] = f"{args.fdr:.12g}"
            selected.append(out)
    fields = list(next(iter(pairs.values())).keys()) + ["rg", "rg_p", "rg_bh_q", "rg_screen_method", "rg_screen_threshold"]
    out_path = PKG / "nasal4_vs_other_mtag_pair_manifest.tsv"
    with out_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(selected)
    print(f"nasal4_vs_other_mtag_pairs selected={len(selected)} total={len(pairs)} valid_ldsc_p={len(valid)} method=BH fdr={args.fdr}")
    if not selected:
        raise SystemExit(1)

if __name__ == "__main__":
    main()
''',
        executable=True,
    )
    write_text(
        PKG / "run_nasal4_vs_other_mtag.sh",
        r'''
#!/usr/bin/env bash
set -euo pipefail

PROJECT_ROOT="${1:-${PROJECT_ROOT:-/platform_data/p_user/p010/phase0}}"
cd "${PROJECT_ROOT}"
source results/phase0_server/server_transfer_package/activate_conda_runtime.sh "${PROJECT_ROOT}"
source results/phase0_server/server_transfer_package/no_conda_server_paths.example.env
: "${MTAG_DIR:?Set MTAG_DIR}"

PKG_DIR="results/phase0_server/nasal4_vs_other_package"
LOG_DIR="${PKG_DIR}/logs"
mkdir -p "${LOG_DIR}" "${PKG_DIR}/inputs" "${PKG_DIR}/mtag"

JOBS="${PHASE0_MTAG_JOBS:-${SLURM_NTASKS:-40}}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
echo "nasal4_vs_other MTAG concurrent pair jobs=${JOBS}"

python3 "${PKG_DIR}/select_nasal4_vs_other_mtag_pairs.py" | tee "${LOG_DIR}/mtag_select_pairs.log"
python3 "${PKG_DIR}/prepare_nasal4_vs_other_mtag_inputs.py" | tee "${LOG_DIR}/mtag_prepare_inputs.log"

run_pair() {
  local pair_id="$1" input_a="$2" input_b="$3" out_prefix="$4"
  mkdir -p "$(dirname "${out_prefix}")"
  if [[ -s "${out_prefix}_trait_1.txt" && -s "${out_prefix}_trait_2.txt" && -s "${out_prefix}.log" ]] && grep -q "MTAG complete" "${out_prefix}.log"; then
    return 0
  fi
  rm -f "${out_prefix}_trait_1.txt" "${out_prefix}_trait_2.txt" "${out_prefix}_omega_hat.txt" "${out_prefix}_sigma_hat.txt"
  /usr/bin/time -v -o "${LOG_DIR}/mtag_${pair_id}.time.log" \
    python "${MTAG_DIR}/mtag.py" \
      --sumstats "${input_a},${input_b}" \
      --out "${out_prefix}" \
      --snp_name SNP --chr_name CHR --bpos_name BP --a1_name A1 --a2_name A2 --eaf_name EAF --z_name Z --p_name P --n_name N \
      --output_pval_max "${PHASE0_MTAG_DISPLAY_P:-1e-5}" \
      --output_anchor_trait 1 \
      --output_manhattan_bin_bp "${PHASE0_MTAG_PLOT_BIN_BP:-1000000}" \
      --output_lead_pval "${PHASE0_MTAG_CANDIDATE_P:-5e-8}" \
      --output_lead_window_bp "${PHASE0_MTAG_LEAD_WINDOW_BP:-500000}" \
      --stream_stdout \
      >"${LOG_DIR}/mtag_${pair_id}.log" 2>&1
}

JOB_STATUS=0
limit_jobs() {
  while (( $(jobs -rp | wc -l) >= JOBS )); do
    wait -n || JOB_STATUS=1
  done
}

while IFS=$'\t' read -r pair_id trait_a trait_b module_a module_b name_a name_b pair_mode source_a source_b common_a common_b ldsc_a ldsc_b input_a input_b out_prefix out_dir rest; do
  limit_jobs
  run_pair "${pair_id}" "${input_a}" "${input_b}" "${out_prefix}" &
done < <(tail -n +2 "${PKG_DIR}/nasal4_vs_other_mtag_pair_manifest.tsv")
wait || JOB_STATUS=1
if (( JOB_STATUS != 0 )); then
  echo "At least one MTAG pair failed; inspect ${LOG_DIR}/mtag_*.log" >&2
  exit 1
fi

python3 results/phase0_server/mtag_package/filter_mtag_candidates.py \
  --pair-manifest "${PKG_DIR}/nasal4_vs_other_mtag_pair_manifest.tsv" \
  --out "${PKG_DIR}/figure_mtag_manhattan.tsv" \
  --summary "${PKG_DIR}/nasal4_vs_other_mtag_filter_summary.tsv" | tee "${LOG_DIR}/mtag_filter_candidates.log"
python3 "${PKG_DIR}/verify_nasal4_vs_other_mtag_outputs.py"
''',
        executable=True,
    )
    write_text(
        PKG / "verify_nasal4_vs_other_mtag_outputs.py",
        r'''
#!/usr/bin/env python3
import csv
from pathlib import Path

PKG = Path("results/phase0_server/nasal4_vs_other_package")
FATAL = ("Traceback", "Analysis terminated from error", "Could not find SIGNED_SUMSTAT")

def read_tsv(path):
    with open(path, encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))

rows = []
for pair in read_tsv(PKG / "nasal4_vs_other_mtag_pair_manifest.tsv"):
    prefix = Path(pair["mtag_out_prefix"])
    log = Path(str(prefix) + ".log")
    t1 = Path(str(prefix) + "_trait_1.txt")
    t2 = Path(str(prefix) + "_trait_2.txt")
    detail = ""
    if not log.exists() or log.stat().st_size == 0:
        detail = "missing_log"
    else:
        text = log.read_text(encoding="utf-8", errors="replace")
        fatal = next((x for x in FATAL if x in text), "")
        if fatal:
            detail = "fatal_log:" + fatal
        elif "MTAG complete" not in text:
            detail = "missing_completion_marker"
        elif "Compact output retained" not in text:
            detail = "missing_compact_output_marker"
    if not detail and (not t1.exists() or t1.stat().st_size == 0 or not t2.exists() or t2.stat().st_size == 0):
        detail = "missing_trait_output"
    rows.append({"output_file": f"mtag_pair::{pair['pair_id']}", "status": "present" if not detail else "failed_mtag", "detail": detail, "bytes": (t1.stat().st_size if t1.exists() else 0) + (t2.stat().st_size if t2.exists() else 0)})
figure = PKG / "figure_mtag_manhattan.tsv"
rows.append({"output_file": "figure_mtag_manhattan.tsv", "status": "present" if figure.exists() and figure.stat().st_size > 0 else "missing_or_empty", "detail": "", "bytes": figure.stat().st_size if figure.exists() else 0})
with (PKG / "nasal4_vs_other_mtag_verify.tsv").open("w", encoding="utf-8", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=["output_file","status","detail","bytes"], delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
failures = [r for r in rows if r["status"] != "present"]
print(f"nasal4_vs_other_mtag_verify checked={len(rows)} failures={len(failures)}")
raise SystemExit(1 if failures else 0)
''',
        executable=True,
    )


def write_lava_scripts() -> None:
    text = OLD_LAVA_R.read_text(encoding="utf-8")
    text = text.replace('"results/phase0_server/local_rg_package"', '"results/phase0_server/nasal4_vs_other_package"')
    text = text.replace("phase0_v3_local_rg_pair_manifest.tsv", "phase0_v3_local_rg_pair_manifest.tsv")
    write_text(PKG / "run_nasal4_vs_other_lava.R", text, executable=True)
    copy_file(OLD_LAVA_VERIFY, PKG / "verify_nasal4_vs_other_lava_outputs.py", executable=True)
    write_text(
        PKG / "run_nasal4_vs_other_lava.sh",
        r'''
#!/usr/bin/env bash
set -euo pipefail

PROJECT_ROOT="${1:-${PROJECT_ROOT:-/platform_data/p_user/p010/phase0}}"
cd "${PROJECT_ROOT}"
source results/phase0_server/server_transfer_package/activate_conda_runtime.sh "${PROJECT_ROOT}"
source results/phase0_server/server_transfer_package/no_conda_server_paths.example.env
: "${LAVA_REF_PREFIX:?Set LAVA_REF_PREFIX}"
: "${LAVA_LOCUS_FILE:?Set LAVA_LOCUS_FILE}"

PKG_DIR="results/phase0_server/nasal4_vs_other_package"
LOG_DIR="${PKG_DIR}/logs"
mkdir -p "${LOG_DIR}" "${PKG_DIR}/lava_parts" "${PKG_DIR}/lava_raw"

JOBS="${PHASE0_LAVA_JOBS:-${SLURM_NTASKS:-22}}"
if (( JOBS > 22 )); then JOBS=22; fi
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
echo "nasal4_vs_other LAVA chromosome jobs=${JOBS}"

: > "${LOG_DIR}/01_lava.log"
running=0
status=0
for chr in {1..22}; do
  Rscript "${PKG_DIR}/run_nasal4_vs_other_lava.R" "${PKG_DIR}" worker "${chr}" >"${LOG_DIR}/lava_chr_${chr}.log" 2>&1 &
  running=$((running + 1))
  if (( running >= JOBS )); then
    wait -n || status=1
    running=$((running - 1))
  fi
done
while (( running > 0 )); do
  wait -n || status=1
  running=$((running - 1))
done
for chr in {1..22}; do
  printf '\n## chr_%s\n' "${chr}" >> "${LOG_DIR}/01_lava.log"
  cat "${LOG_DIR}/lava_chr_${chr}.log" >> "${LOG_DIR}/01_lava.log"
done
if (( status != 0 )); then exit 1; fi
Rscript "${PKG_DIR}/run_nasal4_vs_other_lava.R" "${PKG_DIR}" merge | tee -a "${LOG_DIR}/01_lava.log"
python3 "${PKG_DIR}/verify_nasal4_vs_other_lava_outputs.py" --package-dir "${PKG_DIR}" --out "${PKG_DIR}/nasal4_vs_other_lava_verify.tsv"
''',
        executable=True,
    )


def write_coloc_scripts() -> None:
    write_text(
        PKG / "extract_coloc_regions.py",
        r'''
#!/usr/bin/env python3
import argparse
import csv
import gzip
from collections import defaultdict
from pathlib import Path


def open_text(path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace", newline="")
    return open(path, encoding="utf-8", errors="replace", newline="")


def extract_many(tasks):
    by_source = defaultdict(list)
    counts = {}
    headers = {}
    for task in tasks:
        task["dest"] = Path(task["dest"])
        task["dest"].parent.mkdir(parents=True, exist_ok=True)
        task["start"] = int(task["start"])
        task["end"] = int(task["end"])
        task["chr"] = str(task["chr"]).removeprefix("chr")
        by_source[task["source"]].append(task)
        counts[str(task["dest"])] = 0

    for source, source_tasks in by_source.items():
        tasks_by_chr = defaultdict(list)
        for task in source_tasks:
            tasks_by_chr[task["chr"]].append(task)
        for chrom in tasks_by_chr:
            tasks_by_chr[chrom].sort(key=lambda row: row["start"])

        writers = {}
        handles = {}
        with open_text(source) as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if reader.fieldnames is None:
                continue
            for task in source_tasks:
                headers[str(task["dest"])] = reader.fieldnames
            for row in reader:
                row_chr = (row.get("CHR") or row.get("chr") or row.get("#chrom") or "").removeprefix("chr")
                row_pos = row.get("BP") or row.get("pos") or row.get("position") or ""
                try:
                    pos = int(float(row_pos))
                except ValueError:
                    continue
                for task in tasks_by_chr.get(row_chr, []):
                    if task["start"] > pos:
                        break
                    if pos > task["end"]:
                        continue
                    dest = str(task["dest"])
                    if dest not in writers:
                        out = gzip.open(task["dest"], "wt", encoding="utf-8", newline="")
                        handles[dest] = out
                        writer = csv.DictWriter(out, fieldnames=reader.fieldnames, delimiter="\t", lineterminator="\n")
                        writer.writeheader()
                        writers[dest] = writer
                    writers[dest].writerow(row)
                    counts[dest] += 1
        for out in handles.values():
            out.close()

    for dest in counts:
        path = Path(dest)
        if path.exists():
            continue
        fieldnames = headers.get(dest, ["SNP", "CHR", "BP", "A1", "A2", "BETA", "SE", "P", "N"])
        with gzip.open(path, "wt", encoding="utf-8", newline="") as out:
            writer = csv.DictWriter(out, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
            writer.writeheader()
    return counts


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--locus-manifest", default="results/phase0_server/coloc_package/phase0_v3_coloc_locus_manifest.tsv")
    parser.add_argument("--out", default="results/phase0_server/coloc_package/phase0_v3_region_extract_qc.tsv")
    args = parser.parse_args()
    loci = list(csv.DictReader(open(args.locus_manifest, encoding="utf-8"), delimiter="\t"))
    extract_tasks = []
    for locus in loci:
        for side in ["a", "b"]:
            extract_tasks.append({
                "locus_id": locus["locus_id"],
                "pair_id": locus["pair_id"],
                "side": side,
                "source": locus[f"trait_{side}_file"],
                "dest": locus[f"region_{side}_file"],
                "chr": locus["chr"],
                "start": locus["start"],
                "end": locus["end"],
            })
    counts = extract_many(extract_tasks)

    rows = []
    for locus in loci:
        chrom = locus["chr"]
        start = int(locus["start"])
        end = int(locus["end"])
        n_a = counts.get(locus["region_a_file"], 0)
        n_b = counts.get(locus["region_b_file"], 0)
        rows.append({
            "locus_id": locus["locus_id"],
            "pair_id": locus["pair_id"],
            "window": f"{chrom}:{start}-{end}",
            "trait_a_rows": n_a,
            "trait_b_rows": n_b,
            "status": "present" if n_a > 0 and n_b > 0 else "missing_or_empty",
        })
    fields = ["locus_id", "pair_id", "window", "trait_a_rows", "trait_b_rows", "status"]
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    failures = [row for row in rows if row["status"] != "present"]
    print(f"region_extract loci={len(rows)} failures={len(failures)}")
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
''',
        executable=True,
    )
    copy_file(OLD_COLOC_VERIFY, PKG / "verify_nasal4_vs_other_coloc_outputs.py", executable=True)
    coloc_r_path = PKG / "run_nasal4_vs_other_coloc.R"
    if coloc_r_path.exists():
        coloc_r = coloc_r_path.read_text(encoding="utf-8")
    else:
        coloc_r = OLD_COLOC_R.read_text(encoding="utf-8")
        coloc_r = coloc_r.replace('"results/phase0_server/coloc_package"', '"results/phase0_server/nasal4_vs_other_package"')
        coloc_r = coloc_r.replace('"results/phase0_server/mtag_package/figure_mtag_manhattan.tsv"', '"results/phase0_server/nasal4_vs_other_package/figure_mtag_manhattan.tsv"')
    write_text(PKG / "run_nasal4_vs_other_coloc.R", coloc_r, executable=True)
    build_loci_path = PKG / "build_nasal4_vs_other_coloc_loci.py"
    if not build_loci_path.exists():
        raise FileNotFoundError(build_loci_path)
    write_text(PKG / "build_nasal4_vs_other_coloc_loci.py", build_loci_path.read_text(encoding="utf-8"), executable=True)
    write_text(
        PKG / "run_nasal4_vs_other_coloc.sh",
        r'''
#!/usr/bin/env bash
set -euo pipefail

PROJECT_ROOT="${1:-${PROJECT_ROOT:-/platform_data/p_user/p010/phase0}}"
cd "${PROJECT_ROOT}"
source results/phase0_server/server_transfer_package/activate_conda_runtime.sh "${PROJECT_ROOT}"
source results/phase0_server/server_transfer_package/no_conda_server_paths.example.env
PKG_DIR="results/phase0_server/nasal4_vs_other_package"
LOG_DIR="${PKG_DIR}/logs"
mkdir -p "${LOG_DIR}" "${PKG_DIR}/regions" "${PKG_DIR}/coloc" "${PKG_DIR}/regions_merged" "${PKG_DIR}/coloc_merged"
export PHASE0_LAVA_COLOC_FDR="${PHASE0_LAVA_COLOC_FDR:-0.05}"
echo "PHASE0_LAVA_COLOC_FDR=${PHASE0_LAVA_COLOC_FDR}"

python3 "${PKG_DIR}/build_nasal4_vs_other_coloc_loci.py" | tee "${LOG_DIR}/coloc_build_loci.log"
/usr/bin/time -v -o "${LOG_DIR}/coloc.time.log" \
  bash -c "python3 '${PKG_DIR}/extract_coloc_regions.py' --locus-manifest '${PKG_DIR}/phase0_v3_coloc_locus_manifest.tsv' --out '${PKG_DIR}/nasal4_vs_other_region_extract_qc.tsv' && Rscript '${PKG_DIR}/run_nasal4_vs_other_coloc.R' '${PKG_DIR}'" \
  >"${LOG_DIR}/coloc.log" 2>&1
python3 "${PKG_DIR}/verify_nasal4_vs_other_coloc_outputs.py" --package-dir "${PKG_DIR}" --out "${PKG_DIR}/nasal4_vs_other_coloc_verify.tsv"
python3 "${PKG_DIR}/summarize_nasal4_memory.py" "${LOG_DIR}/coloc.time.log" > "${LOG_DIR}/coloc.memory.tsv"
''',
        executable=True,
    )


def write_slurm_and_docs(trait_count: int, pair_count: int) -> None:
    copy_file(NASAL4_PKG / "summarize_nasal4_memory.py", PKG / "summarize_nasal4_memory.py", executable=True)
    base = r'''
#!/bin/bash
#SBATCH -J {job}
#SBATCH -o %j.{job}.out
#SBATCH -e %j.{job}.err
#SBATCH -D /platform_data/p_user/p010/phase0
#SBATCH -p compute_nodes
#SBATCH -N 1
#SBATCH -n {cpus}
#SBATCH --mem={mem}
#SBATCH -t {time}

set -euo pipefail
export SLURM_CONF="${{SLURM_CONF:-/usr/local/slurm/etc/slurm.conf}}"
export PROJECT_ROOT="${{PROJECT_ROOT:-/platform_data/p_user/p010/phase0}}"
export CONDA_ROOT="${{CONDA_ROOT:-/platform_data/p_user/p010/miniconda3}}"
export PHASE0_CONDA_ENV="${{PHASE0_CONDA_ENV:-/platform_data/p_user/p010/conda_envs/crswnp_phase0}}"
export PHASE0_PY2_ENV="${{PHASE0_PY2_ENV:-/platform_data/p_user/p010/conda_envs/ldsc_mtag_py2}}"
export THREADS="${{SLURM_NTASKS:-{cpus}}}"
export PHASE0_RUN_TMP="${{PROJECT_ROOT}}/.tmp_slurm/${{SLURM_JOB_ID:-{job}}}"
mkdir -p "${{PHASE0_RUN_TMP}}"

echo "START $(date)"
hostname
echo "THREADS=${{THREADS}}"
echo "MEM={mem}"
{exports}
{command} "${{PROJECT_ROOT}}"
echo "DONE $(date)"
'''
    jobs = [
        ("nasal4vo_ldsc", 40, "256G", "2-00:00:00", "export PHASE0_LDSC_JOBS=${THREADS}", "results/phase0_server/nasal4_vs_other_package/run_nasal4_vs_other_ldsc.sh"),
        ("nasal4vo_lava", 22, "256G", "2-00:00:00", "export PHASE0_LAVA_JOBS=${THREADS}", "results/phase0_server/nasal4_vs_other_package/run_nasal4_vs_other_lava.sh"),
        ("nasal4vo_mtag", 40, "400G", "3-00:00:00", "export PHASE0_MTAG_JOBS=${THREADS}", "results/phase0_server/nasal4_vs_other_package/run_nasal4_vs_other_mtag.sh"),
        ("nasal4vo_coloc", 8, "128G", "1-00:00:00", "export PHASE0_COLOC_JOBS=${THREADS}", "results/phase0_server/nasal4_vs_other_package/run_nasal4_vs_other_coloc.sh"),
    ]
    for job, cpus, mem, time, exports, command in jobs:
        write_text(PKG / f"{job}_{cpus}c{mem.lower()}.slurm", base.format(job=job, cpus=cpus, mem=mem, time=time, exports=exports, command=command))
    write_text(
        PKG / "submit_nasal4_vs_other_ldsc_lava.sh",
        r'''
#!/usr/bin/env bash
set -euo pipefail
export PATH=/usr/local/slurm/bin:/usr/local/slurm/sbin:/usr/bin:/bin:$PATH
export SLURM_CONF="${SLURM_CONF:-/usr/local/slurm/etc/slurm.conf}"
cd /platform_data/p_user/p010/phase0/results/phase0_server/nasal4_vs_other_package
COMMON_EXPORTS="ALL,PROJECT_ROOT=/platform_data/p_user/p010/phase0,CONDA_ROOT=/platform_data/p_user/p010/miniconda3,PHASE0_CONDA_ENV=/platform_data/p_user/p010/conda_envs/crswnp_phase0,PHASE0_PY2_ENV=/platform_data/p_user/p010/conda_envs/ldsc_mtag_py2,SLURM_CONF=${SLURM_CONF}"
LDSC_JOBID="$(sbatch --parsable --export="${COMMON_EXPORTS}" nasal4vo_ldsc_40c256g.slurm)"
LAVA_JOBID="$(sbatch --parsable --export="${COMMON_EXPORTS}" nasal4vo_lava_22c256g.slurm)"
echo "ldsc=${LDSC_JOBID} lava=${LAVA_JOBID}"
''',
        executable=True,
    )
    write_text(
        PKG / "submit_nasal4_vs_other_after_ldsc.sh",
        r'''
#!/usr/bin/env bash
set -euo pipefail
if [[ $# -lt 1 ]]; then
  echo "Usage: $0 <ldsc_jobid>" >&2
  exit 2
fi
LDSC_JOBID="$1"
export PATH=/usr/local/slurm/bin:/usr/local/slurm/sbin:/usr/bin:/bin:$PATH
export SLURM_CONF="${SLURM_CONF:-/usr/local/slurm/etc/slurm.conf}"
cd /platform_data/p_user/p010/phase0/results/phase0_server/nasal4_vs_other_package
COMMON_EXPORTS="ALL,PROJECT_ROOT=/platform_data/p_user/p010/phase0,CONDA_ROOT=/platform_data/p_user/p010/miniconda3,PHASE0_CONDA_ENV=/platform_data/p_user/p010/conda_envs/crswnp_phase0,PHASE0_PY2_ENV=/platform_data/p_user/p010/conda_envs/ldsc_mtag_py2,SLURM_CONF=${SLURM_CONF}"
MTAG_JOBID="$(sbatch --parsable --dependency="afterok:${LDSC_JOBID}" --export="${COMMON_EXPORTS}" nasal4vo_mtag_40c400g.slurm)"
COLOC_JOBID="$(sbatch --parsable --dependency="afterok:${MTAG_JOBID}" --export="${COMMON_EXPORTS}" nasal4vo_coloc_8c128g.slurm)"
echo "mtag=${MTAG_JOBID} coloc=${COLOC_JOBID}"
''',
        executable=True,
    )
    write_text(
        PKG / "README.md",
        f'''
# Nasal4-vs-Other Formal Phase 0 Package

This is the formal replacement chain for the project design:

- nasal side: four non-FinnGen replacement chronic nasal inflammatory GWAS.
- comparator side: the selected non-nasal system disease panel from the project
  `nasal_vs_other` manifest.
- pair mode: `nasal4_vs_other`; nasal-nasal and other-other pairs are excluded.

Counts:

- traits: {trait_count}
- formal pairs: {pair_count}

The older `nasal4_package` LDSC/LAVA outputs are input/QC only because they
compare the four nasal traits against each other. They are not formal T4/T5/T6
results.

Slurm entrypoints:

- LDSC: `nasal4vo_ldsc_40c256g.slurm`
- LAVA: `nasal4vo_lava_22c256g.slurm`
- MTAG: `nasal4vo_mtag_40c400g.slurm`
- COLOC: `nasal4vo_coloc_8c128g.slurm`
''',
    )


def main() -> int:
    PKG.mkdir(parents=True, exist_ok=True)
    trait_count, pair_count = build_manifests()
    write_ldsc_scripts()
    write_lava_scripts()
    write_mtag_scripts()
    write_coloc_scripts()
    write_slurm_and_docs(trait_count, pair_count)
    print(f"Built {PKG} traits={trait_count} pairs={pair_count}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
