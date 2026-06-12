#!/usr/bin/env python3
"""Build the five final supplementary table workbooks without external packages."""

from __future__ import annotations

import csv
import math
import os
import re
import zipfile
from collections import OrderedDict
from datetime import datetime, timezone
from pathlib import Path
from xml.sax.saxutils import escape


paper_root = os.environ.get("PAPER_ROOT", "").strip()
results_root = os.environ.get("RESULTS_ROOT", "").strip()
if not paper_root or not results_root:
    raise SystemExit("Set PAPER_ROOT and RESULTS_ROOT before building the supplementary tables.")
ROOT = Path(paper_root)
TABLE_ROOT = Path(os.environ.get("TABLE_OUTPUT_ROOT", ROOT / "Table"))
OUT_DIR = TABLE_ROOT / "supplementary"
SOURCE_OUT = TABLE_ROOT / "source_data" / "final"
QA_DIR = TABLE_ROOT / "qa"

RESULTS = Path(results_root)
ENRICHMENT = (
    RESULTS
    / "final_gtexv8_twas_total_domain_fuma_enrichment"
    / "enrichment_results"
)


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, str]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fields, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows({field: row.get(field, "") for field in fields} for row in rows)


def merge_left(
    left: list[dict[str, str]],
    right: list[dict[str, str]],
    key: str,
    right_fields: list[str],
) -> list[dict[str, str]]:
    index = {row[key]: row for row in right}
    output = []
    for row in left:
        merged = dict(row)
        matched = index.get(row.get(key, ""), {})
        for field in right_fields:
            merged[field] = matched.get(field, "")
        output.append(merged)
    return output


def unique_tokens(values: list[str]) -> list[str]:
    found = set()
    for value in values:
        if not value or value == "NA":
            continue
        for token in value.split(";"):
            token = token.strip()
            if token:
                found.add(token)
    return sorted(found)


def build_st1() -> tuple[list[dict[str, str]], list[str]]:
    metadata = read_tsv(
        ROOT
        / "method_table"
        / "Supplementary_Table_1_GWAS_source_and_phenotype_information.tsv"
    )
    power = read_tsv(
        ROOT / "method_table" / "Supplementary_Table_8_LDSC_h2_power_audit.tsv"
    )
    power_index = {row["trait_id"]: row for row in power}
    fields = [
        "trait_id",
        "trait_label",
        "phenotype_code_or_accession",
        "disease_group",
        "phenotype_module",
        "source_database",
        "consortium_or_biobank",
        "ancestry",
        "sample_size",
        "case_count",
        "control_count",
        "genome_build",
        "effect_scale",
        "ldsc_h2",
        "ldsc_h2_se",
        "ldsc_h2_z",
        "ldsc_mean_chi2",
        "ldsc_intercept",
        "ldsc_power_class",
        "trait_level_h2_qc_decision",
        "appears_in_final_31_pairs",
        "n_final_downstream_pairs",
        "external_metadata_source",
        "notes",
    ]
    rows = []
    for source in metadata:
        audit = power_index.get(source["trait_id"], {})
        row = {field: source.get(field, "") for field in fields}
        row["ldsc_mean_chi2"] = audit.get("mean_chi2", "")
        if not row["ldsc_mean_chi2"] or row["ldsc_mean_chi2"] == "NA":
            row["ldsc_mean_chi2"] = "Not available in the nasal-anchor audit"
        row["ldsc_power_class"] = (
            source.get("ldsc_h2_power_class", "") or audit.get("power_class", "")
        )
        rows.append(row)
    return rows, fields


def build_st2_parameters() -> tuple[list[dict[str, str]], list[str]]:
    software = read_tsv(
        ROOT
        / "method_table"
        / "Supplementary_Table_4_Software_parameters_and_reference_resources.tsv"
    )
    testing = read_tsv(
        ROOT
        / "method_table"
        / "Supplementary_Table_5_Multiple_testing_correction_universe.tsv"
    )
    grading = read_tsv(
        ROOT
        / "method_table"
        / "Supplementary_Table_7_Evidence_grading_definitions.tsv"
    )
    fields = [
        "section",
        "module_or_level",
        "item_or_grade",
        "software_or_package",
        "version",
        "analysis_stratum",
        "parameter_or_required_evidence",
        "reference_or_supporting_evidence",
        "testing_universe",
        "number_of_tests",
        "correction_method",
        "significance_threshold",
        "fdr_threshold",
        "interpretation_boundary",
        "notes",
    ]
    rows = []
    for source in software:
        if source["module"] in {
            "FUMA GENE2FUNC/SNP2GENE",
            "Cell marker ORA / WebCSEA",
        }:
            continue
        rows.append(
            {
                "section": "Software and reference resources",
                "module_or_level": source["module"],
                "software_or_package": source["software_or_package"],
                "version": source["version"],
                "parameter_or_required_evidence": source["main_parameters"],
                "reference_or_supporting_evidence": source[
                    "reference_panel_or_resource"
                ],
                "interpretation_boundary": source["genome_build"],
                "notes": source["notes"],
            }
        )
    rows.extend(
        [
            {
                "section": "Software and reference resources",
                "module_or_level": "LCV",
                "software_or_package": "LCV R implementation",
                "version": "Git commit 39950a8b9c67d812a5f7ec8a1df6951942cb258e",
                "parameter_or_required_evidence": (
                    "RunLCV on signed summary statistics; ancestry-matched LD scores; "
                    "MHC excluded; 20-100 block jackknife blocks according to merged SNP count"
                ),
                "reference_or_supporting_evidence": (
                    "O'Connor and Price LCV implementation; European LD-score reference"
                ),
                "interpretation_boundary": (
                    "GCP is secondary directional evidence and is not a direct causal-effect estimate"
                ),
                "notes": "All 31 final disease pairs completed successfully.",
            },
            {
                "section": "Software and reference resources",
                "module_or_level": "Bidirectional MR",
                "software_or_package": "Custom R workflow; MRPRESSO; PLINK",
                "version": (
                    "R 4.5.3; data.table 1.17.8; MRPRESSO 1.0; "
                    "PLINK v1.9.0-b.8"
                ),
                "parameter_or_required_evidence": (
                    "P <= 5e-8; F statistic >= 10; European LD clumping r2 < 0.001 "
                    "within 10 Mb; palindromic variants removed; >=4 harmonized instruments"
                ),
                "reference_or_supporting_evidence": (
                    "IVW primary estimate; MR-Egger, weighted median and MR-PRESSO "
                    "sensitivity estimates; MR-PRESSO 2,000 simulations, seed 1"
                ),
                "interpretation_boundary": (
                    "Exploratory secondary directionality only; heterogeneous source-GWAS "
                    "effect scales prevent cross-direction effect-size ranking"
                ),
                "notes": "58 of 62 directions were estimable; four were instrument-limited.",
            },
            {
                "section": "Software and reference resources",
                "module_or_level": "Approved drug-target annotation",
                "software_or_package": "DGIdb GraphQL API",
                "version": "Public DGIdb GraphQL API queried 5 June 2026",
                "parameter_or_required_evidence": (
                    "Gene-A or Gene-B prioritized genes queried for drug-gene interactions; "
                    "display restricted to records with DGIdb approved=True"
                ),
                "reference_or_supporting_evidence": (
                    "DGIdb approval flag and drug-gene interaction records"
                ),
                "interpretation_boundary": (
                    "Exploratory target annotation only; not evidence of efficacy, safety "
                    "or a nasal-disease indication"
                ),
                "notes": "127 unique approved-drug annotations after deduplication.",
            },
        ]
    )
    for source in testing:
        if source["module"] in {"FUMA GENE2FUNC", "Cell marker ORA/WebCSEA"}:
            continue
        rows.append(
            {
                "section": "Multiple-testing correction",
                "module_or_level": source["module"],
                "analysis_stratum": source["analysis_stratum"],
                "testing_universe": source["testing_universe"],
                "number_of_tests": source["number_of_tests"],
                "correction_method": source["correction_method"],
                "significance_threshold": source["significance_threshold"],
                "fdr_threshold": source["fdr_threshold"],
                "notes": source["notes"],
            }
        )
    rows.extend(
        [
            {
                "section": "Multiple-testing correction",
                "module_or_level": "LCV",
                "analysis_stratum": "31 final disease pairs",
                "testing_universe": "one GCP test per final disease pair",
                "number_of_tests": "31",
                "correction_method": "Bonferroni used for strict summary",
                "significance_threshold": "nominal P < 0.05; strict P < 0.05/31",
                "fdr_threshold": "Not applicable",
                "interpretation_boundary": "Secondary directional evidence only",
                "notes": "Nominal and Bonferroni results are distinguished.",
            },
            {
                "section": "Multiple-testing correction",
                "module_or_level": "Bidirectional MR",
                "analysis_stratum": "62 prespecified exposure-outcome directions",
                "testing_universe": "two directions for each of 31 final disease pairs",
                "number_of_tests": "62",
                "correction_method": "Bonferroni used for strict IVW summary",
                "significance_threshold": "nominal P < 0.05; strict IVW P < 0.05/62",
                "fdr_threshold": "Not applicable",
                "interpretation_boundary": (
                    "Sensitivity flags and instrument eligibility must accompany P values"
                ),
                "notes": (
                    "MR-Egger intercept, Cochran heterogeneity and MR-PRESSO global tests "
                    "are reported as sensitivity diagnostics."
                ),
            },
        ]
    )
    for source in grading:
        if source["grade"] == "MAGMA-only sensitivity gene" or source[
            "evidence_level"
        ] == "annotation":
            continue
        required = source["required_evidence"].replace("FUMA ", "post hoc ")
        supporting = source["supporting_evidence"].replace(
            "FUMA concordance, and pathway/tissue/cell annotation",
            "post hoc functional annotation",
        )
        not_allowed = (
            source["not_allowed_for_upgrade"]
            .replace("FUMA mapped-gene status, ", "")
            .replace("FUMA mapping, ", "")
            .replace("FUMA annotation, ", "Post hoc functional annotation, ")
        )
        rows.append(
            {
                "section": "Evidence definitions",
                "module_or_level": source["evidence_level"],
                "item_or_grade": source["grade"],
                "parameter_or_required_evidence": required,
                "reference_or_supporting_evidence": supporting,
                "interpretation_boundary": source["interpretation_boundary"],
                "notes": "Not allowed for upgrade: "
                + not_allowed
                + (
                    "; " + source["notes"]
                    if source.get("notes") and source["notes"] != "NA"
                    else ""
                ),
            }
        )
    normalized = []
    for row in rows:
        normalized.append(
            {
                field: row.get(field, "") or "Not applicable"
                for field in fields
            }
        )
    return normalized, fields


def build_st2_qc() -> tuple[list[dict[str, str]], list[str]]:
    rows = read_tsv(
        ROOT
        / "method_table"
        / "Supplementary_Table_3_Summary_statistics_harmonization_and_QC.tsv"
    )
    for row in rows:
        for field, value in row.items():
            row[field] = value.replace(
                "GRCh37/hg19 for GWAS/locus/FUMA coordinate analyses",
                "GRCh37/hg19 for GWAS and locus-coordinate analyses",
            )
    return rows, list(rows[0])


def build_st3() -> tuple[list[dict[str, str]], list[str]]:
    pairs = read_tsv(ROOT / "figure_final" / "source_data" / "pair_summary.tsv")
    trait_metadata = read_tsv(
        ROOT
        / "method_table"
        / "Supplementary_Table_1_GWAS_source_and_phenotype_information.tsv"
    )
    trait_labels = {row["trait_id"]: row["trait_label"] for row in trait_metadata}
    display_labels = {
        "ALLERG_ASTHMA": "Allergic asthma",
        "ASTHMA_EOSINOPHIL_SUGG": "Eosinophilic asthma",
        "AUTOIMMUNE_NONTHYROID": "Non-thyroid autoimmune disease",
        "H7_ALLERGICCONJUNCTIVITIS": "Allergic conjunctivitis",
        "H8_EUSTSALP": "Eustachian salpingitis and obstruction",
        "H8_MIDDLEMASTOID": "Middle-ear and mastoid disease",
        "H8_SUP_ACUTE": "Acute suppurative otitis media",
        "J10_ASTHMA_EXMORE": "Asthma",
        "J10_BRONCHITIS": "Acute bronchitis",
        "J10_CHRONTONSADEN": "Chronic tonsil and adenoid disease",
        "J10_COPD": "Chronic obstructive pulmonary disease",
        "J10_PNEUMOBACT": "Bacterial pneumonia",
        "J10_PNEUMONIA": "Pneumonia",
        "K11_CD_STRICT2": "Crohn disease",
        "K11_COELIAC": "Coeliac disease",
        "K11_IBD_STRICT": "Inflammatory bowel disease",
        "L12_ATOPIC": "Atopic dermatitis",
        "L12_DERMATITISECZEMA": "Dermatitis and eczema",
        "L12_PSORIASIS": "Psoriasis",
        "M13_RHEUMA": "Rheumatoid arthritis",
        "NONALLERG_ASTHMA_EXMORE": "Non-allergic asthma",
        "RHEUMA_SEROPOS_WIDE": "Seropositive rheumatoid arthritis",
    }
    for row in pairs:
        row["pair_name"] = (
            f"{trait_labels.get(row['trait1'], row['base_group'])} - "
            f"{display_labels.get(row['trait2'], trait_labels.get(row['trait2'], row['trait2']))}"
        )
    manifest = read_tsv(
        ROOT
        / "method_table"
        / "Supplementary_Table_2_Disease_pair_analysis_manifest.tsv"
    )
    append_fields = [
        "ldsc_sig_bh05",
        "placo_genomewide_support",
        "cpassoc_genomewide_support",
    ]
    rows = merge_left(pairs, manifest, "pair_id", append_fields)
    fields = [
        "P",
        "pair_order",
        "pair_name",
        "base_group",
        "phenotype_domain",
        "pair_id",
        "trait1",
        "trait2",
        "rg",
        "rg_se",
        "rg_lo",
        "rg_hi",
        "rg_p",
        "rg_q",
        "lava_n_loci",
        "lava_sig_regions",
        "magma_sig_gene_count",
        "route",
        *append_fields,
    ]
    return rows, fields


def build_st4_loci() -> tuple[list[dict[str, str]], list[str]]:
    rows = read_tsv(ROOT / "figure_final" / "source_data" / "locus_genes.tsv")
    for row in rows:
        if not row["nearest_gene"] or row["nearest_gene"] == "NA":
            row["nearest_gene"] = "Not assigned"
        if not row["gene_a_genes"] or row["gene_a_genes"] == "NA":
            row["gene_a_genes"] = "No Gene-A gene assigned"
        if not row["n_gene_a"] or row["n_gene_a"] == "NA":
            row["n_gene_a"] = "0"
        if not row["best_mtag_score"] or row["best_mtag_score"] == "NA":
            row["best_mtag_score"] = "Not available"
    fields = [
        "P",
        "locus_id",
        "chr",
        "cytoband",
        "cytoband_name",
        "mid",
        "genome_pos",
        "locus_grade",
        "support_n",
        "n_pairs_at_locus",
        "nearest_gene",
        "positive_genes",
        "gene_a_genes",
        "n_gene_a",
        "best_mtag_score",
        "pair_id",
    ]
    return rows, fields


def build_st4_genes() -> tuple[list[dict[str, str]], list[str]]:
    rows = read_tsv(
        RESULTS
        / "final_gtexv8_twas_total_domain_fuma_enrichment"
        / "final_gene_tables"
        / "final_twas_supported_unique_genes.tsv"
    )
    return rows, list(rows[0])


def selected_pathway_files() -> list[tuple[str, str, Path]]:
    selected = []
    for analysis in ("primary", "strict"):
        for database in ("GO_BP", "KEGG", "Reactome"):
            path = (
                ENRICHMENT
                / "global"
                / f"enrichment_global_{analysis}{'_123' if analysis == 'primary' else ''}_vs_AB3782_{database}.tsv"
            )
            selected.append((f"global_{analysis}", "global", path))
    domain_dir = ENRICHMENT / "domain"
    pattern = "enrichment_domain_*_vs_global_AB3782_*.tsv"
    for path in sorted(domain_dir.glob(pattern)):
        match = re.match(
            r"enrichment_domain_(.+)_vs_global_AB3782_(GO_BP|KEGG|Reactome)\.tsv",
            path.name,
        )
        if match:
            selected.append(("domain", match.group(1), path))
    return selected


def build_pathways() -> tuple[list[dict[str, str]], list[str]]:
    rows = []
    fields = [
        "analysis_scope",
        "phenotype_domain",
        "database",
        "term_id",
        "term_name",
        "FDR",
        "p_value",
        "odds_ratio",
        "overlap_count",
        "input_gene_count",
        "background_gene_count",
        "term_size_in_background",
        "gene_ratio",
        "background_ratio",
        "overlap_genes",
        "source_analysis",
    ]
    for scope, domain, path in selected_pathway_files():
        if not path.exists():
            raise FileNotFoundError(path)
        for source in read_tsv(path):
            try:
                significant = float(source.get("FDR", "nan")) < 0.05
            except ValueError:
                significant = False
            if not significant:
                continue
            row = {field: source.get(field, "") for field in fields}
            row["analysis_scope"] = scope
            row["phenotype_domain"] = domain
            row["source_analysis"] = source.get("analysis_name", path.stem)
            rows.append(row)
    rows.sort(
        key=lambda row: (
            row["analysis_scope"],
            row["phenotype_domain"],
            row["database"],
            float(row["FDR"]),
        )
    )
    return rows, fields


def build_drugs() -> tuple[list[dict[str, str]], list[str]]:
    source = read_tsv(RESULTS / "network_pharmacology_res" / "drug_recurrence_summary.tsv")
    grouped: OrderedDict[str, list[dict[str, str]]] = OrderedDict()
    for row in source:
        if row["approved_any"] == "True":
            grouped.setdefault(row["drug_name"], []).append(row)
    fields = [
        "drug_name",
        "annotation_source",
        "analysis_date",
        "annotation_scope",
        "n_pairs",
        "pair_list",
        "n_target_genes_total",
        "target_gene_list",
        "approved_drug_annotation",
        "interpretation_boundary",
    ]
    rows = []
    for drug_name, records in grouped.items():
        pairs = unique_tokens([row["pair_list"] for row in records])
        genes = unique_tokens([row["target_gene_list"] for row in records])
        rows.append(
            {
                "drug_name": drug_name,
                "annotation_source": "DGIdb GraphQL API",
                "analysis_date": "2026-06-05",
                "annotation_scope": "Gene-A or Gene-B prioritized targets",
                "n_pairs": str(len(pairs)),
                "pair_list": ";".join(pairs),
                "n_target_genes_total": str(len(genes)),
                "target_gene_list": ";".join(genes),
                "approved_drug_annotation": "True",
                "interpretation_boundary": (
                    "Exploratory target annotation; does not imply therapeutic efficacy "
                    "or clinical indication for nasal disease."
                ),
            }
        )
    rows.sort(key=lambda row: (-int(row["n_pairs"]), row["drug_name"]))
    return rows, fields


def build_directional() -> tuple[list[dict[str, str]], list[str]]:
    pair_summary = read_tsv(ROOT / "figure_final" / "source_data" / "pair_summary.tsv")
    pair_key = {row["pair_id"]: row for row in pair_summary}
    lcv = read_tsv(ROOT / "result" / "Table" / "Table_R14_LCVMR_final31_summary.tsv")
    lcv_key = {row["pair_id"]: row for row in lcv}
    direction_qc = read_tsv(
        ROOT / "result" / "Table" / "Table_R15_MR_direction_QC_final31.tsv"
    )
    sensitivity = read_tsv(
        ROOT / "result" / "Table" / "Table_R16_MR_sensitivity_final31.tsv"
    )
    sensitivity_key = {
        (row["pair_id"], row["direction"]): row for row in sensitivity
    }
    mr_method_rows = read_tsv(
        ROOT / "LCVMR_data" / "results_final31" / "MR" / "mr_results_all.tsv"
    )
    mr_method_key = {
        (row["pair_id"], row["direction"], row["method"]): row
        for row in mr_method_rows
        if row["method"]
    }
    fields = [
        "P",
        "pair_id",
        "pair_name",
        "exposure",
        "outcome",
        "direction",
        "LCV_gcp",
        "LCV_gcp_se",
        "LCV_gcp_p",
        "LCV_status",
        "MR_IVW_beta",
        "MR_IVW_se",
        "MR_IVW_p",
        "MR_IVW_OR",
        "MR_IVW_OR_lci95",
        "MR_IVW_OR_uci95",
        "MR_status",
        "MR_Egger_beta",
        "MR_Egger_se",
        "MR_Egger_p",
        "MR_Egger_OR",
        "MR_Egger_OR_lci95",
        "MR_Egger_OR_uci95",
        "weighted_median_beta",
        "weighted_median_se",
        "weighted_median_p",
        "weighted_median_OR",
        "weighted_median_OR_lci95",
        "weighted_median_OR_uci95",
        "MR_PRESSO_beta",
        "MR_PRESSO_se",
        "MR_PRESSO_p",
        "MR_PRESSO_OR",
        "MR_PRESSO_OR_lci95",
        "MR_PRESSO_OR_uci95",
        "n_instruments_raw",
        "n_instruments_after_strength_filter",
        "n_instruments_after_clump",
        "n_instruments_after_harmonise",
        "n_palindromic_removed",
        "egger_intercept",
        "egger_p",
        "heterogeneity_q",
        "heterogeneity_p",
        "mr_presso_global_p",
        "n_outlier",
        "sensitivity_status",
        "mr_presso_analysis",
        "interpretation_note",
    ]
    rows = []
    for qc in direction_qc:
        pair = pair_key[qc["pair_id"]]
        summary = lcv_key[qc["pair_id"]]
        sens = sensitivity_key.get((qc["pair_id"], qc["direction"]), {})
        method_data = {
            method: mr_method_key.get((qc["pair_id"], qc["direction"], method), {})
            for method in ("IVW", "MR-Egger", "Weighted Median", "MR-PRESSO")
        }
        forward = (
            qc["exposure"] == summary["trait1"] and qc["outcome"] == summary["trait2"]
        )
        prefix = "MR_trait1_to_trait2" if forward else "MR_trait2_to_trait1"
        row = {
            "P": pair["P"],
            "pair_id": qc["pair_id"],
            "pair_name": pair["pair_name"],
            "exposure": qc["exposure"],
            "outcome": qc["outcome"],
            "direction": qc["direction"],
            "LCV_gcp": summary["LCV_gcp"],
            "LCV_gcp_se": summary["LCV_gcp_se"],
            "LCV_gcp_p": summary["LCV_gcp_p"],
            "LCV_status": summary["LCV_status"],
            "MR_IVW_beta": summary[f"{prefix}_IVW_beta"],
            "MR_IVW_se": method_data["IVW"].get("se", ""),
            "MR_IVW_p": summary[f"{prefix}_IVW_p"],
            "MR_IVW_OR": summary[f"{prefix}_IVW_OR"],
            "MR_IVW_OR_lci95": method_data["IVW"].get("or_lci95", ""),
            "MR_IVW_OR_uci95": method_data["IVW"].get("or_uci95", ""),
            "MR_status": summary[f"{prefix}_status"],
            "MR_Egger_beta": method_data["MR-Egger"].get("beta", ""),
            "MR_Egger_se": method_data["MR-Egger"].get("se", ""),
            "MR_Egger_p": method_data["MR-Egger"].get("pval", ""),
            "MR_Egger_OR": method_data["MR-Egger"].get("or", ""),
            "MR_Egger_OR_lci95": method_data["MR-Egger"].get("or_lci95", ""),
            "MR_Egger_OR_uci95": method_data["MR-Egger"].get("or_uci95", ""),
            "weighted_median_beta": method_data["Weighted Median"].get("beta", ""),
            "weighted_median_se": method_data["Weighted Median"].get("se", ""),
            "weighted_median_p": method_data["Weighted Median"].get("pval", ""),
            "weighted_median_OR": method_data["Weighted Median"].get("or", ""),
            "weighted_median_OR_lci95": method_data["Weighted Median"].get(
                "or_lci95", ""
            ),
            "weighted_median_OR_uci95": method_data["Weighted Median"].get(
                "or_uci95", ""
            ),
            "MR_PRESSO_beta": method_data["MR-PRESSO"].get("beta", ""),
            "MR_PRESSO_se": method_data["MR-PRESSO"].get("se", ""),
            "MR_PRESSO_p": method_data["MR-PRESSO"].get("pval", ""),
            "MR_PRESSO_OR": method_data["MR-PRESSO"].get("or", ""),
            "MR_PRESSO_OR_lci95": method_data["MR-PRESSO"].get("or_lci95", ""),
            "MR_PRESSO_OR_uci95": method_data["MR-PRESSO"].get("or_uci95", ""),
            "n_instruments_raw": qc["n_instruments_raw"],
            "n_instruments_after_strength_filter": qc[
                "n_instruments_after_strength_filter"
            ],
            "n_instruments_after_clump": qc["n_instruments_after_clump"],
            "n_instruments_after_harmonise": qc[
                "n_instruments_after_harmonise"
            ],
            "n_palindromic_removed": qc["n_palindromic_removed"],
            "egger_intercept": sens.get("egger_intercept", ""),
            "egger_p": sens.get("egger_p", ""),
            "heterogeneity_q": sens.get("heterogeneity_q", ""),
            "heterogeneity_p": sens.get("heterogeneity_p", ""),
            "mr_presso_global_p": sens.get("mr_presso_global_p", ""),
            "n_outlier": sens.get("n_outlier", ""),
            "sensitivity_status": sens.get("status", ""),
            "mr_presso_analysis": sens.get("mr_presso_analysis", ""),
            "interpretation_note": summary["interpretation_note"],
        }
        if qc["status"] == "insufficient_instruments":
            for field in fields:
                if (
                    field.startswith("MR_")
                    or field.startswith("weighted_median_")
                    or field
                    in {
                        "egger_intercept",
                        "egger_p",
                        "heterogeneity_q",
                        "heterogeneity_p",
                        "mr_presso_global_p",
                        "n_outlier",
                        "mr_presso_analysis",
                    }
                ) and field not in {"MR_status"}:
                    if not row.get(field):
                        row[field] = "Not estimable: insufficient instruments"
        rows.append(row)
    rows.sort(key=lambda row: (int(row["P"][1:]), row["direction"]))
    return rows, fields


def excel_column(index: int) -> str:
    label = ""
    while index:
        index, remainder = divmod(index - 1, 26)
        label = chr(65 + remainder) + label
    return label


NUMERIC_RE = re.compile(r"^-?(?:\d+\.?\d*|\.\d+)(?:[eE][+-]?\d+)?$")


def numeric_value(value: str) -> str | None:
    text = value.strip()
    if not text or text in {"NA", "NaN", "Inf", "-Inf"}:
        return None
    if not NUMERIC_RE.match(text):
        return None
    if len(text) > 1 and text[0] == "0" and text[1].isdigit() and "." not in text:
        return None
    try:
        number = float(text)
    except ValueError:
        return None
    return text if math.isfinite(number) else None


def xml_text(value: str) -> str:
    cleaned = "".join(
        char for char in str(value) if char in "\t\n\r" or ord(char) >= 32
    )
    return escape(cleaned)


def sheet_xml(
    rows: list[dict[str, str]], fields: list[str], widths: list[float]
) -> str:
    row_xml = []
    header_cells = []
    for column, field in enumerate(fields, 1):
        ref = f"{excel_column(column)}1"
        header_cells.append(
            f'<c r="{ref}" s="1" t="inlineStr"><is><t>{xml_text(field)}</t></is></c>'
        )
    row_xml.append(f'<row r="1" ht="34" customHeight="1">{"".join(header_cells)}</row>')
    for row_number, row in enumerate(rows, 2):
        cells = []
        for column, field in enumerate(fields, 1):
            value = str(row.get(field, "") or "")
            if not value:
                continue
            ref = f"{excel_column(column)}{row_number}"
            number = numeric_value(value)
            if number is not None:
                cells.append(f'<c r="{ref}" s="3"><v>{number}</v></c>')
            else:
                preserve = ' xml:space="preserve"' if value != value.strip() else ""
                cells.append(
                    f'<c r="{ref}" s="2" t="inlineStr"><is><t{preserve}>'
                    f"{xml_text(value)}</t></is></c>"
                )
        row_xml.append(f'<row r="{row_number}">{"".join(cells)}</row>')
    columns = "".join(
        f'<col min="{i}" max="{i}" width="{width:.1f}" customWidth="1"/>'
        for i, width in enumerate(widths, 1)
    )
    last_cell = f"{excel_column(len(fields))}{len(rows) + 1}"
    return (
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
        '<worksheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main">'
        f'<dimension ref="A1:{last_cell}"/>'
        '<sheetViews><sheetView workbookViewId="0">'
        '<pane ySplit="1" topLeftCell="A2" activePane="bottomLeft" state="frozen"/>'
        "</sheetView></sheetViews>"
        '<sheetFormatPr defaultRowHeight="15"/>'
        f"<cols>{columns}</cols><sheetData>{''.join(row_xml)}</sheetData>"
        f'<autoFilter ref="A1:{last_cell}"/>'
        '<pageMargins left="0.25" right="0.25" top="0.5" bottom="0.5" '
        'header="0.2" footer="0.2"/>'
        '<pageSetup orientation="landscape" fitToWidth="1" fitToHeight="0"/>'
        "</worksheet>"
    )


def column_widths(rows: list[dict[str, str]], fields: list[str]) -> list[float]:
    widths = []
    for field in fields:
        lengths = [len(field)] + [
            max((len(part) for part in str(row.get(field, "")).splitlines()), default=0)
            for row in rows[:500]
        ]
        widths.append(min(45.0, max(10.0, max(lengths) * 1.05 + 2)))
    return widths


def write_xlsx(
    path: Path,
    sheets: list[tuple[str, list[dict[str, str]], list[str]]],
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    content_types = [
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>',
        '<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types">',
        '<Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/>',
        '<Default Extension="xml" ContentType="application/xml"/>',
        '<Override PartName="/xl/workbook.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet.main+xml"/>',
        '<Override PartName="/xl/styles.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.styles+xml"/>',
        '<Override PartName="/docProps/core.xml" ContentType="application/vnd.openxmlformats-package.core-properties+xml"/>',
        '<Override PartName="/docProps/app.xml" ContentType="application/vnd.openxmlformats-officedocument.extended-properties+xml"/>',
    ]
    for index in range(1, len(sheets) + 1):
        content_types.append(
            f'<Override PartName="/xl/worksheets/sheet{index}.xml" '
            'ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.worksheet+xml"/>'
        )
    content_types.append("</Types>")
    workbook_sheets = "".join(
        f'<sheet name="{xml_text(name)}" sheetId="{index}" r:id="rId{index}"/>'
        for index, (name, _, _) in enumerate(sheets, 1)
    )
    workbook = (
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
        '<workbook xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main" '
        'xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships">'
        f"<sheets>{workbook_sheets}</sheets></workbook>"
    )
    workbook_rels = [
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>',
        '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">',
    ]
    for index in range(1, len(sheets) + 1):
        workbook_rels.append(
            f'<Relationship Id="rId{index}" '
            'Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/worksheet" '
            f'Target="worksheets/sheet{index}.xml"/>'
        )
    workbook_rels.append(
        f'<Relationship Id="rId{len(sheets) + 1}" '
        'Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/styles" '
        'Target="styles.xml"/>'
    )
    workbook_rels.append("</Relationships>")
    styles = (
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
        '<styleSheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main">'
        '<fonts count="2"><font><sz val="10"/><name val="Arial"/></font>'
        '<font><b/><color rgb="FFFFFFFF"/><sz val="10"/><name val="Arial"/></font></fonts>'
        '<fills count="3"><fill><patternFill patternType="none"/></fill>'
        '<fill><patternFill patternType="gray125"/></fill>'
        '<fill><patternFill patternType="solid"><fgColor rgb="FF365F7D"/>'
        '<bgColor indexed="64"/></patternFill></fill></fills>'
        '<borders count="1"><border><left/><right/><top/><bottom/><diagonal/></border></borders>'
        '<cellStyleXfs count="1"><xf numFmtId="0" fontId="0" fillId="0" borderId="0"/></cellStyleXfs>'
        '<cellXfs count="4">'
        '<xf numFmtId="0" fontId="0" fillId="0" borderId="0" xfId="0"/>'
        '<xf numFmtId="0" fontId="1" fillId="2" borderId="0" xfId="0" applyAlignment="1">'
        '<alignment wrapText="1" vertical="center"/></xf>'
        '<xf numFmtId="0" fontId="0" fillId="0" borderId="0" xfId="0" applyAlignment="1">'
        '<alignment vertical="top"/></xf>'
        '<xf numFmtId="0" fontId="0" fillId="0" borderId="0" xfId="0" applyAlignment="1">'
        '<alignment vertical="top"/></xf>'
        "</cellXfs>"
        '<cellStyles count="1"><cellStyle name="Normal" xfId="0" builtinId="0"/></cellStyles>'
        "</styleSheet>"
    )
    now = datetime.now(timezone.utc).replace(microsecond=0).isoformat().replace("+00:00", "Z")
    core = (
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
        '<cp:coreProperties xmlns:cp="http://schemas.openxmlformats.org/package/2006/metadata/core-properties" '
        'xmlns:dc="http://purl.org/dc/elements/1.1/" '
        'xmlns:dcterms="http://purl.org/dc/terms/" '
        'xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance">'
        "<dc:creator>Study authors</dc:creator>"
        "<dc:title>Nasal disease shared genetics supplementary tables</dc:title>"
        f'<dcterms:created xsi:type="dcterms:W3CDTF">{now}</dcterms:created>'
        "</cp:coreProperties>"
    )
    app = (
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
        '<Properties xmlns="http://schemas.openxmlformats.org/officeDocument/2006/extended-properties" '
        'xmlns:vt="http://schemas.openxmlformats.org/officeDocument/2006/docPropsVTypes">'
        "<Application>Python standard library</Application></Properties>"
    )
    root_rels = (
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
        '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">'
        '<Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/officeDocument" Target="xl/workbook.xml"/>'
        '<Relationship Id="rId2" Type="http://schemas.openxmlformats.org/package/2006/relationships/metadata/core-properties" Target="docProps/core.xml"/>'
        '<Relationship Id="rId3" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/extended-properties" Target="docProps/app.xml"/>'
        "</Relationships>"
    )
    with zipfile.ZipFile(path, "w", zipfile.ZIP_DEFLATED) as archive:
        archive.writestr("[Content_Types].xml", "".join(content_types))
        archive.writestr("_rels/.rels", root_rels)
        archive.writestr("docProps/core.xml", core)
        archive.writestr("docProps/app.xml", app)
        archive.writestr("xl/workbook.xml", workbook)
        archive.writestr("xl/_rels/workbook.xml.rels", "".join(workbook_rels))
        archive.writestr("xl/styles.xml", styles)
        for index, (_, rows, fields) in enumerate(sheets, 1):
            archive.writestr(
                f"xl/worksheets/sheet{index}.xml",
                sheet_xml(rows, fields, column_widths(rows, fields)),
            )


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    SOURCE_OUT.mkdir(parents=True, exist_ok=True)
    QA_DIR.mkdir(parents=True, exist_ok=True)
    for old_tsv in SOURCE_OUT.glob("*.tsv"):
        old_tsv.unlink()

    st1_rows, st1_fields = build_st1()
    st2_param_rows, st2_param_fields = build_st2_parameters()
    st2_qc_rows, st2_qc_fields = build_st2_qc()
    st3_rows, st3_fields = build_st3()
    st4_loci_rows, st4_loci_fields = build_st4_loci()
    st4_gene_rows, st4_gene_fields = build_st4_genes()
    pathway_rows, pathway_fields = build_pathways()
    drug_rows, drug_fields = build_drugs()
    direction_rows, direction_fields = build_directional()

    workbooks = [
        (
            "Supplementary_Table_1_GWAS_datasets_and_QC.xlsx",
            [("GWAS_and_power_QC", st1_rows, st1_fields)],
        ),
        (
            "Supplementary_Table_2_Analysis_parameters_and_definitions.xlsx",
            [
                ("Analysis_parameters", st2_param_rows, st2_param_fields),
                ("Harmonization_QC", st2_qc_rows, st2_qc_fields),
            ],
        ),
        (
            "Supplementary_Table_3_Final_31_disease_pairs.xlsx",
            [("Final_31_pairs", st3_rows, st3_fields)],
        ),
        (
            "Supplementary_Table_4_Prioritized_loci_and_genes.xlsx",
            [
                ("Prioritized_loci", st4_loci_rows, st4_loci_fields),
                ("Prioritized_genes", st4_gene_rows, st4_gene_fields),
            ],
        ),
        (
            "Supplementary_Table_5_Downstream_supporting_evidence.xlsx",
            [
                ("Significant_pathways", pathway_rows, pathway_fields),
                ("Approved_drug_targets", drug_rows, drug_fields),
                ("Directional_evidence", direction_rows, direction_fields),
            ],
        ),
    ]

    qa_rows = []
    for workbook_name, sheets in workbooks:
        write_xlsx(OUT_DIR / workbook_name, sheets)
        table_id = f"ST{workbook_name.split('_')[2]}"
        for sheet_name, rows, fields in sheets:
            tsv_name = f"{table_id}_{sheet_name}.tsv"
            write_tsv(SOURCE_OUT / tsv_name, rows, fields)
            qa_rows.append(
                {
                    "workbook": workbook_name,
                    "sheet": sheet_name,
                    "rows": str(len(rows)),
                    "columns": str(len(fields)),
                    "source_data_tsv": str((SOURCE_OUT / tsv_name).relative_to(ROOT)),
                }
            )
    write_tsv(
        QA_DIR / "final_table_inventory.tsv",
        qa_rows,
        ["workbook", "sheet", "rows", "columns", "source_data_tsv"],
    )
    print(f"Built {len(workbooks)} workbooks and {len(qa_rows)} worksheets.")
    for row in qa_rows:
        print(
            f"{row['workbook']} :: {row['sheet']} "
            f"({row['rows']} rows x {row['columns']} columns)"
        )


if __name__ == "__main__":
    main()
