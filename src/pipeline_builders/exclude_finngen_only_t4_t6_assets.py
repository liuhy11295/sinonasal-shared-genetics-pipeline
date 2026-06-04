#!/usr/bin/env python3
"""Exclude legacy FinnGen-only T4-T6 outputs from manuscript-facing indexes."""

from __future__ import annotations

import csv
import shutil
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]
ASSET_ROOT = PROJECT_ROOT / "results" / "manuscript_assets_current_20260525"
EXCLUDED_ROOT = PROJECT_ROOT / "results" / "deprecated_finngen_only_20260528"
ARCHIVED_ASSET_ROOT = EXCLUDED_ROOT / "manuscript_assets_removed"
T5_DOWNLOAD_ROOT = PROJECT_ROOT / "results" / "downloads" / "phase0_t5_lava_20260527"


EXCLUDED_ASSET_PREFIXES = {
    "main_figures/Fig2_global_ldsc_available_layer/": "T4_LDSC",
    "supplementary_figures/SuppFig_S1_ldsc_h2_qc/": "T4_LDSC",
    "supplementary_tables/SuppTable_S3_ldsc_rg_qc/": "T4_LDSC",
}

EXCLUDED_ASSET_DIRS = {
    Path("main_figures/Fig2_global_ldsc_available_layer"): "T4_LDSC",
    Path("supplementary_figures/SuppFig_S1_ldsc_h2_qc"): "T4_LDSC",
    Path("supplementary_tables/SuppTable_S3_ldsc_rg_qc"): "T4_LDSC",
}

DEPRECATION_MARKER = "DEPRECATED_FINNGEN_ONLY_DO_NOT_USE.md"


def read_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return list(reader.fieldnames or []), list(reader)


def write_tsv(path: Path, fields: list[str], rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text.rstrip() + "\n", encoding="utf-8")


def unique_destination(path: Path) -> Path:
    if not path.exists():
        return path
    stem = path.stem
    suffix = path.suffix
    for idx in range(1, 1000):
        candidate = path.with_name(f"{stem}.archived{idx}{suffix}")
        if not candidate.exists():
            return candidate
    raise RuntimeError(f"Could not create unique archive destination for {path}")


def assert_inside(path: Path, root: Path) -> None:
    path.resolve().relative_to(root.resolve())


def archive_legacy_asset_dirs(exclusions: list[dict[str, object]]) -> None:
    for rel_dir, layer in EXCLUDED_ASSET_DIRS.items():
        src_dir = ASSET_ROOT / rel_dir
        if not src_dir.exists():
            continue
        assert_inside(src_dir, ASSET_ROOT)

        moved = 0
        for path in sorted(src_dir.rglob("*")):
            if not path.is_file():
                continue
            if path.name == DEPRECATION_MARKER:
                continue
            rel_path = path.relative_to(ASSET_ROOT)
            dst = unique_destination(ARCHIVED_ASSET_ROOT / rel_path)
            dst.parent.mkdir(parents=True, exist_ok=True)
            shutil.move(str(path), str(dst))
            moved += 1
            exclusions.append(
                {
                    "layer": layer,
                    "path": str(dst),
                    "scope": "current_manuscript_tree_archive",
                    "action": "moved_from_current_tree_to_deprecated_archive",
                    "reason": f"Legacy FinnGen-only result file removed from current manuscript assets; original relative path was {rel_path.as_posix()}.",
                }
            )

        write_text(
            src_dir / DEPRECATION_MARKER,
            f"""# Deprecated FinnGen-Only Result Directory

Legacy FinnGen-only result files from this current manuscript asset directory
were moved out of the formal current-results tree on 2026-05-28.

Archived location:

`results/deprecated_finngen_only_20260528/manuscript_assets_removed/{rel_dir.as_posix()}`

Do not cite the old T4/T5/T6 outputs from this directory. Use the replacement
nasal4-vs-other GWAS Catalog/Pan-UKB Slurm chain instead.
""",
        )


def add_archived_asset_exclusions(exclusions: list[dict[str, object]]) -> None:
    if not ARCHIVED_ASSET_ROOT.exists():
        return
    for path in sorted(ARCHIVED_ASSET_ROOT.rglob("*")):
        if not path.is_file():
            continue
        rel = path.relative_to(ARCHIVED_ASSET_ROOT)
        layer = "T4_LDSC"
        if "local_rg" in rel.as_posix().lower() or "lava" in rel.as_posix().lower():
            layer = "T5_LAVA"
        if "mtag" in rel.as_posix().lower():
            layer = "T6_MTAG"
        exclusions.append(
            {
                "layer": layer,
                "path": str(path),
                "scope": "current_manuscript_tree_archive",
                "action": "moved_from_current_tree_to_deprecated_archive",
                "reason": f"Legacy FinnGen-only result file removed from current manuscript assets; archived relative path is {rel.as_posix()}.",
            }
        )


def excluded_layer(asset_file: str) -> str | None:
    for prefix, layer in EXCLUDED_ASSET_PREFIXES.items():
        if asset_file.startswith(prefix):
            return layer
    if asset_file == "manuscript_current_result_counts.tsv":
        return "T4_LDSC"
    return None


def update_manifest(exclusions: list[dict[str, object]]) -> None:
    manifest = ASSET_ROOT / "manuscript_asset_manifest.tsv"
    fields, rows = read_tsv(manifest)
    retained: list[dict[str, object]] = []
    for row in rows:
        layer = excluded_layer(row.get("asset_file", ""))
        if layer is None:
            retained.append(row)
            continue
        exclusions.append(
            {
                "layer": layer,
                "path": str(ASSET_ROOT / row["asset_file"]),
                "scope": "manuscript_asset_manifest",
                "action": "removed_from_current_manifest",
                "reason": "Legacy FinnGen-only T4-T6 result layer excluded from formal manuscript-facing outputs.",
            }
        )
    write_tsv(manifest, fields, retained)


def update_availability(exclusions: list[dict[str, object]]) -> None:
    availability = ASSET_ROOT / "availability_status.tsv"
    fields, rows = read_tsv(availability)
    replacements = {
        ("Fig2", "global LDSC rg layer"): (
            "deprecated_excluded_finngen_only",
            "",
            "Legacy T4 used FinnGen-only same-source disease GWAS; formal global LDSC results await the nasal4 GWAS Catalog/Pan-UKB rerun.",
        ),
        ("SuppFig_S1", "LDSC h2 QC"): (
            "deprecated_excluded_finngen_only",
            "",
            "Legacy T4 heritability QC belongs to the FinnGen-only panel and is excluded from current formal results.",
        ),
        ("SuppFig_S2-S5", "local/partitioned rg"): (
            "deprecated_excluded_finngen_only",
            "",
            "Legacy T5 LAVA local-rg outputs/comparisons were generated from the FinnGen-only panel and are excluded.",
        ),
        ("SuppFig_S6", "MTAG Manhattan"): (
            "deprecated_excluded_finngen_only",
            "",
            "Legacy T6 MTAG result layer/recommendations are excluded; MTAG should use the replacement nasal4 chain.",
        ),
    }
    for row in rows:
        key = (row.get("manuscript_item", ""), row.get("component", ""))
        if key in replacements:
            old = dict(row)
            status, available_file, reason = replacements[key]
            row["status"] = status
            row["available_file"] = available_file
            row["reason_or_gap"] = reason
            exclusions.append(
                {
                    "layer": key[1],
                    "path": old.get("available_file", ""),
                    "scope": "availability_status",
                    "action": "status_changed_to_deprecated_excluded_finngen_only",
                    "reason": reason,
                }
            )
    write_tsv(availability, fields, rows)

    unavailable = ASSET_ROOT / "unavailable_planned_assets" / "planned_but_not_yet_available.tsv"
    unavailable_rows = [row for row in rows if row.get("status") != "available"]
    write_tsv(unavailable, fields, unavailable_rows)


def update_counts() -> None:
    counts = ASSET_ROOT / "manuscript_current_result_counts.tsv"
    rows = [
        {"metric": "formal_t4_t5_t6_result_layers_included", "value": "0"},
        {"metric": "excluded_legacy_result_layers", "value": "T4_LDSC;T5_LAVA;T6_MTAG"},
        {"metric": "exclusion_reason", "value": "legacy_finngen_only_same_source_panel"},
        {"metric": "replacement_chain", "value": "results/phase0_server/nasal4_vs_other_package"},
        {"metric": "replacement_status", "value": "pending_slurm_completion"},
    ]
    write_tsv(counts, ["metric", "value"], rows)


def add_t5_download_exclusions(exclusions: list[dict[str, object]]) -> None:
    if not T5_DOWNLOAD_ROOT.exists():
        return
    for path in sorted(T5_DOWNLOAD_ROOT.rglob("*")):
        if not path.is_file():
            continue
        if path.suffix.lower() not in {".tsv", ".csv", ".md", ".txt"}:
            continue
        layer = "T5_LAVA"
        if "mtag" in path.name.lower():
            layer = "T6_MTAG_RECOMMENDATION"
        exclusions.append(
            {
                "layer": layer,
                "path": str(path),
                "scope": "downloaded_legacy_result_bundle",
                "action": "deprecated_excluded_from_formal_results",
                "reason": "Downloaded legacy FinnGen-only result/comparison bundle retained only for audit.",
            }
        )


def write_deprecation_docs(exclusions: list[dict[str, object]]) -> None:
    manifest = EXCLUDED_ROOT / "excluded_t4_t5_t6_finngen_only_manifest.tsv"
    fields = ["layer", "path", "scope", "action", "reason"]
    write_tsv(manifest, fields, exclusions)
    write_text(
        EXCLUDED_ROOT / "README.md",
        """# Deprecated FinnGen-Only T4-T6 Outputs

This folder records result layers removed from current manuscript-facing use on
2026-05-28.

The excluded layers are:

- T4 LDSC global genetic-correlation and h2-QC result tables from the old
  FinnGen-only panel.
- T5 LAVA local-rg outputs and LDSC/LAVA comparison tables downloaded from the
  old FinnGen-only run.
- T6 MTAG result/recommendation layer from the same old panel, where present.

Original files were not deleted. They are retained only for audit and debugging.
Any old result files previously left in the current manuscript-asset tree were
moved into `manuscript_assets_removed/` under this deprecated folder.
Formal results should be taken from the replacement nasal4-vs-other chain under:

`results/phase0_server/nasal4_vs_other_package`
""",
    )
    write_text(
        T5_DOWNLOAD_ROOT / "DEPRECATED_FINNGEN_ONLY_DO_NOT_USE.md",
        """# Deprecated: FinnGen-Only T5/T6 Bundle

Do not use this downloaded bundle as a formal manuscript result. It came from
the old FinnGen-only T5 LAVA analysis and related MTAG-inclusion comparison
files. Keep it only for audit/debugging.

Use the replacement nasal4-vs-other GWAS Catalog/Pan-UKB Slurm chain instead.
""",
    )
    write_text(
        ASSET_ROOT / "README.md",
        """# Current Manuscript Assets

Legacy FinnGen-only T4/T5/T6 result layers have been excluded from current
manuscript-facing use on 2026-05-28.

What remains here is either workflow/provenance material or non-result planning
context. Do not cite the old T4 LDSC, T5 LAVA, or T6 MTAG outputs from this
folder as formal biological results.

The replacement formal analysis chain is:

`results/phase0_server/nasal4_vs_other_package`

See the exclusion audit table:

`results/deprecated_finngen_only_20260528/excluded_t4_t5_t6_finngen_only_manifest.tsv`
""",
    )


def main() -> None:
    if not ASSET_ROOT.exists():
        raise SystemExit(f"Missing asset root: {ASSET_ROOT}")
    exclusions: list[dict[str, object]] = []
    update_manifest(exclusions)
    update_availability(exclusions)
    update_counts()
    archive_legacy_asset_dirs(exclusions)
    add_archived_asset_exclusions(exclusions)
    add_t5_download_exclusions(exclusions)
    write_deprecation_docs(exclusions)
    print(f"Excluded legacy FinnGen-only entries: {len(exclusions)}")
    print(EXCLUDED_ROOT / "excluded_t4_t5_t6_finngen_only_manifest.tsv")


if __name__ == "__main__":
    main()
