#!/usr/bin/env python3
import csv
import gzip
import hashlib
import os
import shutil
import subprocess
import sys
from datetime import datetime
from pathlib import Path


if not os.environ.get("BASE") or not os.environ.get("ARCHIVE"):
    raise SystemExit("Set BASE and ARCHIVE for the final 31-pair LCV/MR analysis.")
BASE = Path(os.environ["BASE"])
OUT = Path(os.environ.get("OUT", str(BASE / "results_final31")))
ARCHIVE = Path(os.environ["ARCHIVE"])
MANIFEST_DIR = BASE / "manifests"
LCV_DIR = BASE / "lcv"
MR_DIR = BASE / "mr"


def read_tsv(path):
    if not path.exists() or path.stat().st_size == 0:
        return []
    with path.open(encoding="utf-8-sig", newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def write_tsv(path, rows, fields):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, delimiter="\t", fieldnames=fields, extrasaction="ignore")
        w.writeheader()
        for row in rows:
            w.writerow({field: row.get(field, "") for field in fields})


def open_text(path):
    return gzip.open(path, "rt", encoding="utf-8", errors="replace") if str(path).endswith(".gz") else path.open(encoding="utf-8", errors="replace")


def count_rows(path):
    if not path.exists() or not path.is_file():
        return ""
    try:
        opener = gzip.open if str(path).endswith(".gz") else open
        with opener(path, "rt", encoding="utf-8", errors="replace") as fh:
            n = sum(1 for _ in fh)
        return max(0, n - 1)
    except Exception as exc:
        return f"unreadable:{exc}"


def header_cols(path):
    if not path.exists() or not path.is_file():
        return ""
    try:
        with open_text(path) as fh:
            header = fh.readline().strip()
        sep = "\t" if "\t" in header else ","
        return ";".join(header.split(sep)[:30])
    except Exception as exc:
        return f"unreadable:{exc}"


def shell_bool(cmd):
    try:
        return subprocess.run(cmd, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL).returncode == 0
    except Exception:
        return False


def backup_existing_outputs():
    existing = [
        p for p in OUT.iterdir()
        if p.name not in {"scripts"} and not p.name.startswith("backup_") and p.exists()
    ] if OUT.exists() else []
    if not existing:
        return ""
    ts = datetime.now().strftime("%Y%m%d_%H%M%S")
    bdir = OUT / f"backup_{ts}"
    bdir.mkdir(parents=True, exist_ok=True)
    for p in existing:
        dst = bdir / p.name
        if p.is_dir():
            shutil.copytree(p, dst, dirs_exist_ok=True)
        else:
            shutil.copy2(p, dst)
    return str(bdir)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    ARCHIVE.mkdir(parents=True, exist_ok=True)
    backup_dir = backup_existing_outputs()

    final_pairs_path = MANIFEST_DIR / "final_31_pairs.tsv"
    trait_manifest_path = MANIFEST_DIR / "final_24_traits_input_manifest.tsv"
    raw_manifest_path = MR_DIR / "raw_gwas" / "MR_raw_GWAS_download_manifest.tsv"
    pairs = read_tsv(final_pairs_path)
    traits = {r["trait_id"]: r for r in read_tsv(trait_manifest_path)}
    raw = {r["trait_id"]: r for r in read_tsv(raw_manifest_path)}

    rows = []
    for r in pairs:
        t1, t2 = r["trait1"], r["trait2"]
        t1_lcv = BASE / traits.get(t1, {}).get("lcv_input_file", "")
        t2_lcv = BASE / traits.get(t2, {}).get("lcv_input_file", "")
        t1_raw = MR_DIR / "raw_gwas" / "downloads" / Path(raw.get(t1, {}).get("mr_download_url", "")).name
        t2_raw = MR_DIR / "raw_gwas" / "downloads" / Path(raw.get(t2, {}).get("mr_download_url", "")).name
        t1_std = OUT / "MR" / "standardized" / f"{t1}.standardized.tsv.gz"
        t2_std = OUT / "MR" / "standardized" / f"{t2}.standardized.tsv.gz"
        rows.append({
            **r,
            "trait1_lcv_file": str(t1_lcv),
            "trait2_lcv_file": str(t2_lcv),
            "trait1_file": str(t1_lcv),
            "trait2_file": str(t2_lcv),
            "trait1_raw_gwas_url": raw.get(t1, {}).get("mr_download_url", ""),
            "trait2_raw_gwas_url": raw.get(t2, {}).get("mr_download_url", ""),
            "trait1_raw_gwas_file": str(t1_raw),
            "trait2_raw_gwas_file": str(t2_raw),
            "trait1_standard_gwas": str(t1_std),
            "trait2_standard_gwas": str(t2_std),
            "trait1_lcv_exists": str(t1_lcv.exists()),
            "trait2_lcv_exists": str(t2_lcv.exists()),
            "trait1_mr_source": raw.get(t1, {}).get("source", ""),
            "trait2_mr_source": raw.get(t2, {}).get("source", ""),
        })
    fields = list(rows[0].keys()) if rows else ["pair_id", "trait1", "trait2"]
    write_tsv(OUT / "00_pair_manifest_final31.tsv", rows, fields)

    inv_targets = []
    for p in sorted(MANIFEST_DIR.glob("*.tsv")):
        inv_targets.append((p, "manifest", "pair/trait/source manifest"))
    for p in sorted((LCV_DIR / "input_sumstats").glob("*.sumstats.gz")):
        inv_targets.append((p, "lcv_sumstats", "LCV input"))
    for p in sorted((LCV_DIR / "reference_ldscores").glob("LDscore.*.l2.ldscore.gz")):
        inv_targets.append((p, "lcv_ldscore", "LCV LD-score reference"))
    inv_targets.append((raw_manifest_path, "mr_manifest", "MR GWAS download manifest"))
    for p in sorted((MR_DIR / "reference_plink").glob("*.bed")):
        inv_targets.append((p, "mr_ld_bed", "MR clumping PLINK reference"))
    for p in sorted((MR_DIR / "scripts").glob("*")):
        if p.is_file():
            inv_targets.append((p, "script", "historical project script"))
    inv_rows = []
    for p, ftype, usage in inv_targets:
        inv_rows.append({
            "file_path": str(p),
            "file_type": ftype,
            "exists": str(p.exists()),
            "line_count": count_rows(p),
            "key_columns": header_cols(p),
            "usage": usage,
        })
    write_tsv(OUT / "00_input_inventory.tsv", inv_rows, ["file_path", "file_type", "exists", "line_count", "key_columns", "usage"])

    with (OUT / "00_column_mapping_notes.md").open("w", encoding="utf-8") as fh:
        fh.write("# LCVMR final31 column mapping notes\n\n")
        fh.write(f"- Selected pair manifest: `{final_pairs_path}` because it explicitly contains final 31 pairs.\n")
        fh.write("- `trait1` and `trait2` are retained from the selected pair manifest.\n")
        fh.write("- `trait1_file` and `trait2_file` point to local LCV munged sumstats for traceability.\n")
        fh.write("- `trait*_standard_gwas` points to the standardized MR GWAS files to be generated under `results_final31/MR/standardized/`.\n")
        fh.write("- MR raw URLs are read from `mr/raw_gwas/MR_raw_GWAS_download_manifest.tsv`; no alternative data source is selected here.\n")
        if backup_dir:
            fh.write(f"- Existing outputs were copied to backup directory `{backup_dir}` before rewriting standardized final31 outputs.\n")
        fh.write("- LCV package availability was checked locally; absence is recorded in LCV QC rather than silently rerunning with a different method.\n")

    lcv_tool = shell_bool("Rscript -e 'quit(status=!requireNamespace(\"LCV\", quietly=TRUE))'")
    run_lcv_file = next(BASE.glob("**/RunLCV.R"), None)
    lcv_qc = [{
        "check": "selected_pair_manifest_rows",
        "value": str(len(pairs)),
        "status": "ok" if len(pairs) == 31 else "error",
        "notes": str(final_pairs_path),
    }, {
        "check": "lcv_input_sumstats_count",
        "value": str(len(list((LCV_DIR / "input_sumstats").glob("*.sumstats.gz")))),
        "status": "ok",
        "notes": "Expected 24 unique traits.",
    }, {
        "check": "R_LCV_package_available",
        "value": str(lcv_tool),
        "status": "ok" if lcv_tool else "failed",
        "notes": "LCV cannot be rerun locally unless the LCV R package or RunLCV.R is installed.",
    }, {
        "check": "RunLCV_R_found",
        "value": str(run_lcv_file or ""),
        "status": "ok" if run_lcv_file else "failed",
        "notes": "No local RunLCV.R found under LCVMR_data.",
    }]
    write_tsv(OUT / "LCV" / "lcv_qc_summary.tsv", lcv_qc, ["check", "value", "status", "notes"])
    with (OUT / "LCV" / "lcv_run_log.txt").open("w", encoding="utf-8") as fh:
        fh.write(f"{datetime.now().isoformat()} Prepared LCV inventory for final31.\n")
        if lcv_tool or run_lcv_file:
            fh.write("LCV execution tool appears available; run script can be added if needed.\n")
        else:
            fh.write("LCV execution skipped: local LCV package and RunLCV.R were not found. Existing historical results will be normalized with missing final31 pairs marked failed.\n")

    old_lcv_path = OUT / "lcv_results.tsv"
    old_lcv = {r["pair_id"]: r for r in read_tsv(old_lcv_path)}
    lcv_rows = []
    for r in rows:
        old = old_lcv.get(r["pair_id"], {})
        if old:
            status, notes = "historical_result_available", "Historical LCV result copied from pre-existing results_final31/lcv_results.tsv."
        else:
            status, notes = "failed", "No historical LCV result for this final31 pair; local LCV execution package/script unavailable."
        lcv_rows.append({
            "pair_id": r["pair_id"],
            "trait1": r["trait1"],
            "trait2": r["trait2"],
            "gcp": old.get("gcp", ""),
            "gcp_se": old.get("gcp_se", ""),
            "gcp_p": old.get("p_value", old.get("gcp_p", "")),
            "rho_estimate": "",
            "genetic_correlation": "",
            "h2_trait1": "",
            "h2_trait2": "",
            "status": status,
            "notes": notes,
        })
    write_tsv(OUT / "LCV" / "lcv_pair_results.tsv", lcv_rows, ["pair_id", "trait1", "trait2", "gcp", "gcp_se", "gcp_p", "rho_estimate", "genetic_correlation", "h2_trait1", "h2_trait2", "status", "notes"])

    source_notes = OUT / "MR" / "mr_data_source_decisions.md"
    source_notes.parent.mkdir(parents=True, exist_ok=True)
    with source_notes.open("w", encoding="utf-8") as fh:
        fh.write("# MR data source decisions\n\n")
        fh.write("- Primary source: existing project manifest `mr/raw_gwas/MR_raw_GWAS_download_manifest.tsv`.\n")
        fh.write("- No alternate GWAS source was substituted.\n")
        fh.write("- LCV munged sumstats are not used for MR because they contain `Z` but not raw `BETA/SE/P`.\n")

    print(f"Prepared final31 inputs at {OUT}")
    print(f"Pair manifest rows: {len(rows)}")
    print(f"Backup directory: {backup_dir or 'none'}")


if __name__ == "__main__":
    main()
