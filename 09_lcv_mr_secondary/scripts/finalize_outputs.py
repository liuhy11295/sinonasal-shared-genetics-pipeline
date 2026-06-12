#!/usr/bin/env python3

import csv
import hashlib
import os
import shutil
from collections import Counter
from datetime import datetime
from pathlib import Path


if not os.environ.get("BASE") or not os.environ.get("ARCHIVE"):
    raise SystemExit("Set BASE and ARCHIVE for the final 31-pair LCV/MR analysis.")
BASE = Path(os.environ["BASE"])
OUT = Path(os.environ.get("OUT", BASE / "results_final31"))
ARCHIVE = Path(os.environ["ARCHIVE"])


def read_tsv(path):
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path, fieldnames, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)


def md5sum(path, chunk_size=8 * 1024 * 1024):
    digest = hashlib.md5()
    with path.open("rb") as handle:
        while chunk := handle.read(chunk_size):
            digest.update(chunk)
    return digest.hexdigest()


def display_number(value):
    if value in (None, ""):
        return "NA"
    try:
        return f"{float(value):.3g}"
    except ValueError:
        return value


def build_summary():
    pairs = read_tsv(OUT / "00_pair_manifest_final31.tsv")
    lcv = {row["pair_id"]: row for row in read_tsv(OUT / "LCV/lcv_pair_results.tsv")}
    mr_rows = read_tsv(OUT / "MR/mr_results_all.tsv")
    ivw = {
        (row["pair_id"], row["exposure"], row["outcome"]): row
        for row in mr_rows
        if row.get("method") == "IVW"
    }
    instruments = {
        (row["pair_id"], row["exposure"], row["outcome"]): row
        for row in read_tsv(OUT / "MR/mr_instrument_summary.tsv")
    }

    fields = [
        "pair_id",
        "trait1",
        "trait2",
        "LCV_gcp",
        "LCV_gcp_se",
        "LCV_gcp_p",
        "LCV_status",
        "MR_trait1_to_trait2_IVW_beta",
        "MR_trait1_to_trait2_IVW_p",
        "MR_trait1_to_trait2_IVW_OR",
        "MR_trait1_to_trait2_status",
        "MR_trait2_to_trait1_IVW_beta",
        "MR_trait2_to_trait1_IVW_p",
        "MR_trait2_to_trait1_IVW_OR",
        "MR_trait2_to_trait1_status",
        "interpretation_note",
    ]
    output = []
    direction_qc = []
    for pair in pairs:
        pair_id, trait1, trait2 = pair["pair_id"], pair["trait1"], pair["trait2"]
        lcv_row = lcv.get(pair_id, {})
        direction_rows = []
        for exposure, outcome in ((trait1, trait2), (trait2, trait1)):
            key = (pair_id, exposure, outcome)
            inst = instruments.get(key, {})
            mr = ivw.get(key, {})
            status = inst.get("status", "failed")
            direction_rows.append((mr, status))
            direction_qc.append(
                {
                    "pair_id": pair_id,
                    "exposure": exposure,
                    "outcome": outcome,
                    "direction": f"{exposure}_to_{outcome}",
                    "n_instruments_raw": inst.get("n_instruments_raw", ""),
                    "n_instruments_after_strength_filter": inst.get(
                        "n_instruments_after_strength_filter", ""
                    ),
                    "n_instruments_after_clump": inst.get(
                        "n_instruments_after_clump", ""
                    ),
                    "n_instruments_after_harmonise": inst.get(
                        "n_instruments_after_harmonise", ""
                    ),
                    "n_palindromic_removed": inst.get("n_palindromic_removed", ""),
                    "status": status,
                    "notes": mr.get("notes", ""),
                }
            )

        nominal = []
        for label, (mr, status) in zip(("trait1_to_trait2", "trait2_to_trait1"), direction_rows):
            if status == "eligible" and mr.get("pval"):
                try:
                    if float(mr["pval"]) < 0.05:
                        nominal.append(label)
                except ValueError:
                    pass
        if nominal:
            mr_note = "nominal IVW P<0.05 in " + ", ".join(nominal)
        else:
            mr_note = "no eligible direction with nominal IVW P<0.05"
        lcv_note = lcv_row.get("notes", "LCV unavailable")
        note = (
            f"LCV: {lcv_note}; MR: {mr_note}. "
            "Exploratory directional evidence only; evidence grades unchanged."
        )
        (mr12, status12), (mr21, status21) = direction_rows
        output.append(
            {
                "pair_id": pair_id,
                "trait1": trait1,
                "trait2": trait2,
                "LCV_gcp": lcv_row.get("gcp", ""),
                "LCV_gcp_se": lcv_row.get("gcp_se", ""),
                "LCV_gcp_p": lcv_row.get("gcp_p", ""),
                "LCV_status": lcv_row.get("status", "failed"),
                "MR_trait1_to_trait2_IVW_beta": mr12.get("beta", ""),
                "MR_trait1_to_trait2_IVW_p": mr12.get("pval", ""),
                "MR_trait1_to_trait2_IVW_OR": mr12.get("or", ""),
                "MR_trait1_to_trait2_status": status12,
                "MR_trait2_to_trait1_IVW_beta": mr21.get("beta", ""),
                "MR_trait2_to_trait1_IVW_p": mr21.get("pval", ""),
                "MR_trait2_to_trait1_IVW_OR": mr21.get("or", ""),
                "MR_trait2_to_trait1_status": status21,
                "interpretation_note": note,
            }
        )

    write_tsv(OUT / "LCVMR_final31_summary.tsv", fields, output)
    write_tsv(
        OUT / "MR/mr_direction_qc_summary.tsv",
        list(direction_qc[0]),
        direction_qc,
    )
    return pairs, output, direction_qc


def build_raw_inventory():
    manifest = read_tsv(BASE / "mr/raw_gwas/MR_raw_GWAS_download_manifest.tsv")
    rows = []
    for row in manifest:
        url = row["mr_download_url"].strip('"')
        path = BASE / "mr/raw_gwas/downloads" / Path(url).name
        exists = path.is_file() and path.stat().st_size > 0
        rows.append(
            {
                "trait_id": row["trait_id"].strip('"'),
                "source": row["source"].strip('"'),
                "url": url,
                "local_file": str(path),
                "exists": str(exists).lower(),
                "size_bytes": path.stat().st_size if exists else 0,
                "md5": md5sum(path) if exists else "",
            }
        )
    write_tsv(OUT / "MR/mr_raw_gwas_file_inventory.tsv", list(rows[0]), rows)
    return rows


def update_input_inventory(raw_rows):
    inventory_path = OUT / "00_input_inventory.tsv"
    fields = [
        "file_path",
        "file_type",
        "exists",
        "line_count",
        "key_columns",
        "usage",
    ]
    existing = read_tsv(inventory_path)
    retained = [
        row
        for row in existing
        if row.get("file_type") not in {"mr_raw_gwas", "mr_standardized_gwas"}
    ]
    additions = []
    for row in raw_rows:
        additions.append(
            {
                "file_path": row["local_file"],
                "file_type": "mr_raw_gwas",
                "exists": row["exists"],
                "line_count": "",
                "key_columns": (
                    "SNP;CHR;BP;A1;A2;BETA;SE;P;FRQ;N "
                    "(source columns mapped during standardization)"
                ),
                "usage": f'MR raw GWAS; source={row["source"]}; md5={row["md5"]}',
            }
        )
        trait_id = row["trait_id"]
        standardized = OUT / "MR/standardized" / f"{trait_id}.standardized.tsv.gz"
        additions.append(
            {
                "file_path": str(standardized),
                "file_type": "mr_standardized_gwas",
                "exists": str(standardized.is_file() and standardized.stat().st_size > 0),
                "line_count": "",
                "key_columns": "SNP;CHR;BP;A1;A2;BETA;SE;P;FRQ;N",
                "usage": "standardized local MR input",
            }
        )
    write_tsv(inventory_path, fields, retained + additions)


def archive_outputs():
    excluded_parts = {
        "backup_20260607_235212",
        "standardized",
        "tmp",
        "downloads",
        "instruments",
        "outcome_lookup",
        "presso_checkpoints",
    }
    if ARCHIVE.exists() and any(ARCHIVE.iterdir()):
        stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        backup = ARCHIVE.parent / f"LCVMR_backup_{stamp}"
        shutil.move(str(ARCHIVE), str(backup))
    ARCHIVE.mkdir(parents=True, exist_ok=True)

    copied = []
    for source in sorted(OUT.rglob("*")):
        if not source.is_file():
            continue
        relative = source.relative_to(OUT)
        if any(part in excluded_parts or part.startswith("backup_") for part in relative.parts):
            continue
        destination = ARCHIVE / relative
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, destination)
        copied.append((source, destination))
    return copied


def write_archive_index(copied):
    rows = []
    for source, destination in copied:
        rows.append(
            {
                "archive_path": str(destination),
                "size_bytes": destination.stat().st_size,
                "md5": md5sum(destination),
                "source_path": str(source),
            }
        )
    write_tsv(ARCHIVE / "archive_index.tsv", list(rows[0]), rows)


def write_report(pairs, summary, direction_qc, raw_rows, archive_verified):
    lcv_counts = Counter(row["LCV_status"] for row in summary)
    mr_counts = Counter(row["status"] for row in direction_qc)
    lcv_direction_counts = Counter()
    for row in read_tsv(OUT / "LCV/lcv_pair_results.tsv"):
        note = row.get("notes", "")
        if "trait1_to_trait2" in note:
            lcv_direction_counts["trait1_to_trait2"] += 1
        elif "trait2_to_trait1" in note:
            lcv_direction_counts["trait2_to_trait1"] += 1
        else:
            lcv_direction_counts["undetermined"] += 1
    sensitivity = read_tsv(OUT / "MR/mr_sensitivity_all.tsv")
    heterogeneity_sig = sum(
        bool(row.get("heterogeneity_p"))
        and float(row["heterogeneity_p"]) < 0.05
        for row in sensitivity
    )
    egger_sig = sum(
        bool(row.get("egger_p")) and float(row["egger_p"]) < 0.05
        for row in sensitivity
    )
    presso_global_sig = sum(
        bool(row.get("mr_presso_global_p"))
        and float(row["mr_presso_global_p"]) < 0.05
        for row in sensitivity
    )
    eligible_ivw = []
    for row in summary:
        for suffix, direction in (
            ("trait1_to_trait2", f'{row["trait1"]} -> {row["trait2"]}'),
            ("trait2_to_trait1", f'{row["trait2"]} -> {row["trait1"]}'),
        ):
            p = row[f"MR_{suffix}_IVW_p"]
            if p:
                eligible_ivw.append((float(p), row["pair_id"], direction))
    nominal = [row for row in eligible_ivw if row[0] < 0.05]

    main_files = [
        "LCV/lcv_pair_results.tsv",
        "LCV/lcv_qc_summary.tsv",
        "MR/mr_results_all.tsv",
        "MR/mr_sensitivity_all.tsv",
        "MR/mr_instrument_summary.tsv",
        "MR/mr_direction_qc_summary.tsv",
        "MR/mr_qc_summary.tsv",
        "MR/mr_presso_qc.tsv",
        "MR/mr_raw_gwas_file_inventory.tsv",
        "LCVMR_final31_summary.tsv",
    ]
    downloaded = [
        f'- `{row["trait_id"]}`: `{row["local_file"]}` ({row["size_bytes"]} bytes)'
        for row in raw_rows
    ]
    nominal_lines = [
        f"- `{pair_id}` ({direction}), IVW P={display_number(p)}"
        for p, pair_id, direction in sorted(nominal)
    ] or ["- None at nominal IVW P<0.05."]
    report = f"""# LCVMR final31 completion report

## Pair manifest

- Selected manifest: `{OUT / "00_pair_manifest_final31.tsv"}`
- Source manifest: `{BASE / "manifests/final_31_pairs.tsv"}`
- Confirmed disease-pair count: {len(pairs)}
- Selection basis: the manifest is explicitly named `final_31_pairs.tsv` and matches the integrated final pair summary.

## LCV

- Successful pairs: {lcv_counts.get("success", 0)}
- Failed pairs: {len(pairs) - lcv_counts.get("success", 0)}
- Trait1-to-trait2 partial genetic causality: {lcv_direction_counts.get("trait1_to_trait2", 0)}
- Trait2-to-trait1 partial genetic causality: {lcv_direction_counts.get("trait2_to_trait1", 0)}
- Undetermined/shared-genetics-only pairs: {lcv_direction_counts.get("undetermined", 0)}
- Result: `{OUT / "LCV/lcv_pair_results.tsv"}`

## Bidirectional MR

- Expected directions: {len(pairs) * 2}
- Eligible/successful directions: {mr_counts.get("eligible", 0)}
- Insufficient-instrument directions: {mr_counts.get("insufficient_instruments", 0)}
- Failed directions: {mr_counts.get("failed", 0)}
- Instrument definition: P <= 5e-8, F-statistic >= 10, EUR LD clumping r2 < 0.001 within 10 Mb, palindromic variants removed, minimum four harmonized instruments.
- Primary estimator: IVW. MR-Egger, weighted median and MR-PRESSO were sensitivity analyses.
- Sensitivity flags among eligible directions: heterogeneity P<0.05 in {heterogeneity_sig}; Egger-intercept P<0.05 in {egger_sig}; MR-PRESSO global P<0.05 in {presso_global_sig}.
- Directional analyses are exploratory and do not modify SNP, locus or gene evidence grades.

Nominal IVW findings:

{chr(10).join(nominal_lines)}

## Downloaded MR data

All 24 trait-level raw GWAS files were downloaded from the URLs in the MR source manifest. Large raw files were not duplicated in the archive; their paths, sizes and MD5 values are recorded in `MR/mr_raw_gwas_file_inventory.tsv`.

{chr(10).join(downloaded)}

## Main results

{chr(10).join(f"- `{OUT / path}`" for path in main_files)}

## Archive

- Archive directory: `{ARCHIVE}`
- Archive completed: yes
- Key-result MD5 comparison between primary output and archive: {"passed" if archive_verified else "failed"}
- Archive index: `{ARCHIVE / "archive_index.tsv"}`

## Manual review

- Treat nominal MR P values as exploratory; no multiple-testing correction was used to redefine the primary evidence hierarchy.
- Review directions with significant heterogeneity, Egger intercept or MR-PRESSO global tests in `MR/mr_sensitivity_all.tsv`.
- MR-PRESSO uses 2,000 simulation draws; boundary-form empirical P values are retained in the sensitivity table.
- LCV and MR support directional interpretation only and do not establish a causal variant, gene, mechanism or therapeutic effect.
"""
    (OUT / "TASK_LCVMR_final31_completion_report.md").write_text(
        report, encoding="utf-8"
    )


def verify_key_files():
    keys = [
        "00_pair_manifest_final31.tsv",
        "LCV/lcv_pair_results.tsv",
        "LCV/lcv_qc_summary.tsv",
        "MR/mr_results_all.tsv",
        "MR/mr_sensitivity_all.tsv",
        "MR/mr_instrument_summary.tsv",
        "MR/mr_direction_qc_summary.tsv",
        "MR/mr_qc_summary.tsv",
        "MR/mr_raw_gwas_file_inventory.tsv",
        "LCVMR_final31_summary.tsv",
        "TASK_LCVMR_final31_completion_report.md",
    ]
    return all(
        (OUT / relative).is_file()
        and (ARCHIVE / relative).is_file()
        and md5sum(OUT / relative) == md5sum(ARCHIVE / relative)
        for relative in keys
    )


def main():
    pairs, summary, direction_qc = build_summary()
    if len(pairs) != 31 or len(summary) != 31 or len(direction_qc) != 62:
        raise RuntimeError("final31 row-count invariant failed")
    raw_rows = build_raw_inventory()
    if len(raw_rows) != 24 or not all(row["exists"] == "true" for row in raw_rows):
        raise RuntimeError("MR raw GWAS inventory is incomplete")
    update_input_inventory(raw_rows)

    write_report(pairs, summary, direction_qc, raw_rows, archive_verified=False)
    copied = archive_outputs()
    write_archive_index(copied)
    initial_verified = verify_key_files()
    write_report(pairs, summary, direction_qc, raw_rows, archive_verified=initial_verified)
    shutil.copy2(
        OUT / "TASK_LCVMR_final31_completion_report.md",
        ARCHIVE / "TASK_LCVMR_final31_completion_report.md",
    )
    write_archive_index(
        [
            (source, ARCHIVE / source.relative_to(OUT))
            for source, _ in copied
        ]
    )
    if not verify_key_files():
        raise RuntimeError("primary/archive MD5 verification failed")
    print(f"summary_rows={len(summary)}")
    print(f"mr_direction_rows={len(direction_qc)}")
    print(f"archive_files={len(copied)}")
    print("archive_md5_verification=passed")


if __name__ == "__main__":
    main()
