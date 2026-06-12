#!/usr/bin/env python3
"""Submission-readiness audit for the five final supplementary tables."""

from __future__ import annotations

import csv
import os
import importlib.util
import zipfile
from collections import Counter
from pathlib import Path
from xml.etree import ElementTree as ET


paper_root = os.environ.get("PAPER_ROOT", "").strip()
results_root = os.environ.get("RESULTS_ROOT", "").strip()
if not paper_root or not results_root:
    raise SystemExit("Set PAPER_ROOT and RESULTS_ROOT before auditing the supplementary tables.")
ROOT = Path(paper_root)
TABLE_ROOT = Path(os.environ.get("TABLE_OUTPUT_ROOT", ROOT / "Table"))
FINAL = TABLE_ROOT / "source_data" / "final"
XLSX_DIR = TABLE_ROOT / "supplementary"
REPORT = TABLE_ROOT / "qa" / "submission_readiness_audit.tsv"
MISSING = {"", "NA", "NaN", "N/A", "null", "None"}


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def load_builder():
    path = TABLE_ROOT / "scripts" / "build_final_tables.py"
    spec = importlib.util.spec_from_file_location("table_builder", path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader
    spec.loader.exec_module(module)
    return module


def add(results, check, passed, detail):
    results.append(
        {
            "check": check,
            "status": "PASS" if passed else "FAIL",
            "detail": detail,
        }
    )


def main() -> None:
    results = []
    xlsx_files = sorted(XLSX_DIR.glob("Supplementary_Table_*.xlsx"))
    tsv_files = sorted(FINAL.glob("ST*.tsv"))
    add(results, "Final workbook count", len(xlsx_files) == 5, str(len(xlsx_files)))
    add(results, "Final worksheet TSV count", len(tsv_files) == 9, str(len(tsv_files)))

    expected_sheets = {
        "Supplementary_Table_1_GWAS_datasets_and_QC.xlsx": ["GWAS_and_power_QC"],
        "Supplementary_Table_2_Analysis_parameters_and_definitions.xlsx": [
            "Analysis_parameters",
            "Harmonization_QC",
        ],
        "Supplementary_Table_3_Final_31_disease_pairs.xlsx": ["Final_31_pairs"],
        "Supplementary_Table_4_Prioritized_loci_and_genes.xlsx": [
            "Prioritized_loci",
            "Prioritized_genes",
        ],
        "Supplementary_Table_5_Downstream_supporting_evidence.xlsx": [
            "Significant_pathways",
            "Approved_drug_targets",
            "Directional_evidence",
        ],
    }
    namespace = {"m": "http://schemas.openxmlformats.org/spreadsheetml/2006/main"}
    for path in xlsx_files:
        try:
            with zipfile.ZipFile(path) as archive:
                workbook = ET.fromstring(archive.read("xl/workbook.xml"))
                names = [
                    node.attrib["name"]
                    for node in workbook.find("m:sheets", namespace)
                ]
            add(
                results,
                f"Workbook structure: {path.name}",
                names == expected_sheets[path.name],
                ", ".join(names),
            )
        except Exception as error:
            add(results, f"Workbook structure: {path.name}", False, repr(error))

    total_missing = 0
    for path in tsv_files:
        rows = read_tsv(path)
        missing = sum(
            value.strip() in MISSING
            for row in rows
            for value in row.values()
        )
        total_missing += missing
        widths = {len(row) for row in rows}
        add(
            results,
            f"No blank/NA cells: {path.name}",
            missing == 0 and widths == {len(rows[0])},
            f"rows={len(rows)}; columns={len(rows[0])}; missing={missing}",
        )
    add(results, "Total blank/NA cells", total_missing == 0, str(total_missing))

    builder = load_builder()
    expected_builds = {
        "ST1_GWAS_and_power_QC.tsv": builder.build_st1(),
        "ST2_Analysis_parameters.tsv": builder.build_st2_parameters(),
        "ST2_Harmonization_QC.tsv": builder.build_st2_qc(),
        "ST3_Final_31_pairs.tsv": builder.build_st3(),
        "ST4_Prioritized_loci.tsv": builder.build_st4_loci(),
        "ST4_Prioritized_genes.tsv": builder.build_st4_genes(),
        "ST5_Significant_pathways.tsv": builder.build_pathways(),
        "ST5_Approved_drug_targets.tsv": builder.build_drugs(),
        "ST5_Directional_evidence.tsv": builder.build_directional(),
    }
    for name, (expected_rows, expected_fields) in expected_builds.items():
        observed = read_tsv(FINAL / name)
        normalized_expected = [
            {field: row.get(field, "") for field in expected_fields}
            for row in expected_rows
        ]
        add(
            results,
            f"Source-derived content match: {name}",
            observed == normalized_expected,
            f"observed={len(observed)}; expected={len(normalized_expected)}",
        )

    pairs = read_tsv(FINAL / "ST3_Final_31_pairs.tsv")
    pair_ids = {row["pair_id"] for row in pairs}
    p_labels = {row["P"] for row in pairs}
    add(
        results,
        "Final pair set",
        len(pairs) == 31
        and len(pair_ids) == 31
        and p_labels == {f"P{i}" for i in range(1, 32)},
        "31 unique pairs with P1-P31",
    )

    loci = read_tsv(FINAL / "ST4_Prioritized_loci.tsv")
    genes = read_tsv(FINAL / "ST4_Prioritized_genes.tsv")
    add(
        results,
        "Prioritized locus set",
        len(loci) == 181
        and len({(row["P"], row["locus_id"]) for row in loci}) == 181
        and {row["P"] for row in loci} <= p_labels,
        "181 unique pair-locus records",
    )
    add(
        results,
        "Prioritized gene set",
        len(genes) == 123 and len({row["gene_symbol"] for row in genes}) == 123,
        "123 unique genes",
    )

    pathways = read_tsv(FINAL / "ST5_Significant_pathways.tsv")
    add(
        results,
        "Pathway significance",
        len(pathways) == 305 and all(float(row["FDR"]) < 0.05 for row in pathways),
        "305 rows; all FDR < 0.05",
    )

    drugs = read_tsv(FINAL / "ST5_Approved_drug_targets.tsv")
    raw_drugs = read_tsv(
        Path(os.environ["RESULTS_ROOT"])
        / "network_pharmacology_res"
        / "all_dgidb_interactions.tsv"
    )
    approved_raw = {
        row["drug_name"] for row in raw_drugs if row["approved"] == "True"
    }
    add(
        results,
        "Approved drug annotations",
        len(drugs) == 127
        and {row["drug_name"] for row in drugs} == approved_raw
        and all(row["annotation_source"] == "DGIdb GraphQL API" for row in drugs),
        "127 unique DGIdb approved=True drugs",
    )

    directions = read_tsv(FINAL / "ST5_Directional_evidence.tsv")
    direction_counts = Counter(row["pair_id"] for row in directions)
    method_prefixes = ("MR_IVW", "MR_Egger", "weighted_median", "MR_PRESSO")
    eligible = [row for row in directions if row["MR_status"] == "eligible"]
    limited = [
        row for row in directions if row["MR_status"] == "insufficient_instruments"
    ]
    methods_complete = all(
        all(not row[f"{prefix}_beta"].startswith("Not estimable") for prefix in method_prefixes)
        for row in eligible
    ) and all(
        all(row[f"{prefix}_beta"].startswith("Not estimable") for prefix in method_prefixes)
        for row in limited
    )
    add(
        results,
        "Directional evidence coverage",
        len(directions) == 62
        and len(direction_counts) == 31
        and set(direction_counts.values()) == {2}
        and len(eligible) == 58
        and len(limited) == 4
        and methods_complete,
        "62 directions; 58 estimable; 4 explicitly instrument-limited; four MR methods",
    )

    prohibited = ("FUMA", "Cell marker", "WebCSEA", "MAGMA-only")
    prohibited_hits = []
    for path in tsv_files:
        text = path.read_text(encoding="utf-8")
        prohibited_hits.extend(term for term in prohibited if term.lower() in text.lower())
    add(
        results,
        "Deferred modules excluded",
        not prohibited_hits,
        "none" if not prohibited_hits else ", ".join(sorted(set(prohibited_hits))),
    )

    REPORT.parent.mkdir(parents=True, exist_ok=True)
    with REPORT.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["check", "status", "detail"],
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(results)

    failed = [row for row in results if row["status"] != "PASS"]
    print(f"Submission audit: {len(results) - len(failed)}/{len(results)} checks passed")
    for row in failed:
        print(f"FAIL: {row['check']}: {row['detail']}")
    if failed:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
