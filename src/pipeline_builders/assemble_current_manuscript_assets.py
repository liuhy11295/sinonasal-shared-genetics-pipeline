#!/usr/bin/env python3
"""Assemble manuscript-ready assets that already exist after the partial phase0 run."""

from __future__ import annotations

import argparse
import csv
import math
import shutil
from pathlib import Path


DEFAULT_SOURCE = Path(r"D:\phase0_server_downloads\extracted_phase0_partial_T1_T4_20260525_082056")


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


def write_text(path: Path, content: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content.rstrip() + "\n", encoding="utf-8")


def parse_float(value: str | None) -> float:
    try:
        result = float(str(value))
    except (TypeError, ValueError):
        return math.nan
    return result


def finite(value: float) -> bool:
    return math.isfinite(value)


def fmt(value: float) -> str:
    return "NA" if not finite(value) else f"{value:.8g}"


def bh_adjust(p_values: list[float]) -> list[float]:
    n = len(p_values)
    order = sorted(range(n), key=lambda index: p_values[index])
    adjusted = [1.0] * n
    previous = 1.0
    for rank in range(n, 0, -1):
        index = order[rank - 1]
        current = min(previous, p_values[index] * n / rank, 1.0)
        adjusted[index] = current
        previous = current
    return adjusted


def description_path(data_path: Path) -> Path:
    return data_path.with_name(f"{data_path.stem}.DESCRIPTION.md")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source-root", type=Path, default=DEFAULT_SOURCE)
    parser.add_argument("--outdir", type=Path)
    args = parser.parse_args()

    project_root = Path(__file__).resolve().parents[1]
    source_root = args.source_root.resolve()
    outdir = (args.outdir or project_root / "results" / "manuscript_assets_current_20260525").resolve()

    server = source_root / "results" / "phase0_server"
    metadata = source_root / "results" / "phase0_shared_genetics" / "figure_metadata"
    disease = source_root / "results" / "phase0_shared_genetics" / "disease_tables"
    contract = project_root / "results" / "phase0_server" / "plot_contract_package" / "phase0_plot_output_contract.tsv"

    required = [
        server / "ldsc_package" / "figure_global_rg_summary.tsv",
        server / "ldsc_package" / "figure_sinonasal_ldsc_h2_qc.tsv",
        server / "phase0_v3_munge_qc.tsv",
        server / "download_package" / "phase0_v3_download_verify.tsv",
        metadata / "figure_workflow_nodes.tsv",
        metadata / "figure_workflow_edges.tsv",
        metadata / "figure1_workflow.mmd",
        disease / "phase0_active_v3_panel.tsv",
        disease / "phase0_nasal_vs_other_pair_manifest.tsv",
        contract,
    ]
    missing = [str(path) for path in required if not path.exists()]
    if missing:
        raise SystemExit("Required input missing:\n" + "\n".join(missing))

    manifest: list[dict[str, str]] = []

    def register(
        data_path: Path,
        manuscript_role: str,
        asset_type: str,
        status: str,
        source: str,
        explanation: str,
        limitations: str,
    ) -> None:
        desc = description_path(data_path)
        relative = data_path.relative_to(outdir).as_posix()
        write_text(
            desc,
            f"""# {data_path.name}

## 用途

{explanation}

## 来源与处理

来源：`{source}`  
状态：`{status}`  
文章位置：`{manuscript_role}`

## 使用限制

{limitations}
""",
        )
        manifest.append(
            {
                "asset_file": relative,
                "description_file": desc.relative_to(outdir).as_posix(),
                "manuscript_role": manuscript_role,
                "asset_type": asset_type,
                "status": status,
                "source": source,
            }
        )

    def copy_registered(
        source: Path,
        destination: Path,
        manuscript_role: str,
        asset_type: str,
        status: str,
        explanation: str,
        limitations: str,
    ) -> None:
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, destination)
        register(
            destination,
            manuscript_role,
            asset_type,
            status,
            str(source),
            explanation,
            limitations,
        )

    # Figure 1: design/workflow assets already exist independently of the failed downstream step.
    fig1 = outdir / "main_figures" / "Fig1_workflow"
    copy_registered(
        metadata / "figure_workflow_nodes.tsv",
        fig1 / "figure_workflow_nodes.tsv",
        "Main Fig1",
        "figure_input",
        "available",
        "定义工作流图中的阶段节点、输入资源、分析方法和产物，用于绘制研究设计流程图。",
        "这是研究流程和产物规划信息，不是效应估计结果。",
    )
    copy_registered(
        metadata / "figure_workflow_edges.tsv",
        fig1 / "figure_workflow_edges.tsv",
        "Main Fig1",
        "figure_input",
        "available",
        "定义 Fig1 节点之间的数据流连接关系，与节点表共同用于生成流程图。",
        "仅表达流程关系；不能用来判断某个下游分析已经成功完成。",
    )
    copy_registered(
        metadata / "figure1_workflow.mmd",
        fig1 / "figure1_workflow.mmd",
        "Main Fig1",
        "render_specification",
        "available",
        "已有 Mermaid 流程图定义，可作为 Fig1 的可编辑渲染草稿。",
        "其中列出的后续模块属于设计范围；模块是否完成应以可用性清单为准。",
    )

    # Global rg: construct an audit-ready table without overwriting the original server output.
    rg_fields, rg_rows = read_tsv(server / "ldsc_package" / "figure_global_rg_summary.tsv")
    rg_log_dir = server / "ldsc_package" / "rg"
    planned_n = len(rg_rows)
    reviewed: list[dict[str, object]] = []
    p_for_family: list[float] = []
    for row in rg_rows:
        pair_id = f"{row['trait_a']}__{row['trait_b']}"
        log_path = rg_log_dir / f"{pair_id}.log"
        log_text = log_path.read_text(encoding="utf-8", errors="replace") if log_path.exists() else ""
        rg = parse_float(row.get("rg"))
        se = parse_float(row.get("se"))
        p = parse_float(row.get("p"))
        if "ERROR computing rg" in log_text or "Traceback" in log_text:
            status = "analysis_error"
            reason = "LDSC log reports computation error/traceback"
        elif not all(finite(value) for value in (rg, se, p)):
            status = "invalid_rg_nonfinite"
            reason = "rg, se or p is missing/non-finite"
        elif abs(rg) > 1:
            status = "invalid_rg_out_of_bounds"
            reason = "absolute rg exceeds 1"
        elif se < 0 or p < 0 or p > 1:
            status = "invalid_statistic"
            reason = "se or p is outside valid range"
        else:
            status = "valid"
            reason = ""
        valid = status == "valid"
        p_for_family.append(p if valid else 1.0)
        outrow: dict[str, object] = dict(row)
        outrow.update(
            {
                "pair_id": pair_id,
                "validation_status": status,
                "exclusion_reason": reason,
                "manuscript_include": str(valid),
                "multiple_testing_family_n": planned_n,
            }
        )
        reviewed.append(outrow)
    adjusted = bh_adjust(p_for_family)
    for row, q_value in zip(reviewed, adjusted):
        row["p_adj_bh_planned_family"] = fmt(q_value)

    reviewed_fields = rg_fields + [
        "pair_id",
        "validation_status",
        "exclusion_reason",
        "manuscript_include",
        "multiple_testing_family_n",
        "p_adj_bh_planned_family",
    ]
    valid_rg = [row for row in reviewed if row["validation_status"] == "valid"]
    significant_rg = [
        row for row in valid_rg if parse_float(str(row["p_adj_bh_planned_family"])) < 0.05
    ]
    excluded_rg = [row for row in reviewed if row["validation_status"] != "valid"]

    fig2 = outdir / "main_figures" / "Fig2_global_ldsc_available_layer"
    reviewed_path = fig2 / "figure_global_rg_summary_reviewed.tsv"
    significant_path = fig2 / "figure_global_rg_significant.tsv"
    write_tsv(reviewed_path, reviewed_fields, reviewed)
    write_tsv(significant_path, reviewed_fields, significant_rg)
    register(
        reviewed_path,
        "Main Fig2 (global LDSC layer) / Supplementary full data",
        "reviewed_result_table",
        "available_after_qc_review",
        str(server / "ldsc_package" / "figure_global_rg_summary.tsv"),
        "全局遗传相关结果审阅版。保留全部计划配对，并新增有效性状态、排除原因和按全部计划配对校正的 BH q 值，可用于绘制 Fig2 的 LDSC 层及作为补充全量结果表。",
        "Fig2 计划中的候选位点或局部相关层尚未产生；原始输出的 `p_adj` 全为 NA，文章引用应使用新增的 `p_adj_bh_planned_family`。无效结果不得作生物学解释。",
    )
    register(
        significant_path,
        "Main Fig2 (global LDSC layer)",
        "plotting_subset",
        "available_after_qc_review",
        str(reviewed_path),
        "经过有效性检查且在全部计划配对检验族中满足 BH q < 0.05 的全局遗传相关配对，可直接作为 Fig2 显著边或热图标注输入。",
        "该表仅代表全局 LDSC 显著信号，不代表已通过局部相关、共定位或其他下游验证。",
    )

    supp_rg = outdir / "supplementary_tables" / "SuppTable_S3_ldsc_rg_qc"
    excluded_path = supp_rg / "ldsc_rg_excluded_pairs.tsv"
    counts_path = supp_rg / "ldsc_rg_qc_counts.tsv"
    write_tsv(excluded_path, reviewed_fields, excluded_rg)
    status_order = [
        "valid",
        "analysis_error",
        "invalid_rg_nonfinite",
        "invalid_rg_out_of_bounds",
        "invalid_statistic",
    ]
    qc_counts = [
        {"validation_status": status, "n_pairs": sum(row["validation_status"] == status for row in reviewed)}
        for status in status_order
        if any(row["validation_status"] == status for row in reviewed)
    ]
    qc_counts.append({"validation_status": "significant_valid_bh_q_lt_0.05", "n_pairs": len(significant_rg)})
    write_tsv(counts_path, ["validation_status", "n_pairs"], qc_counts)
    register(
        excluded_path,
        "Supplementary Table S3",
        "qc_exclusion_table",
        "available_after_qc_review",
        str(reviewed_path),
        "列出不能进入正文解释的 LDSC `rg` 配对，并记录具体排除原因，供审稿和复核。",
        "被排除结果可以报告为 QC 事件，但不得作为相关方向或显著性的证据。",
    )
    register(
        counts_path,
        "Supplementary Table S3",
        "qc_summary_table",
        "available_after_qc_review",
        str(reviewed_path),
        "汇总 LDSC 配对的可用、异常及显著结果数量，便于补充材料中透明报告分析过滤流程。",
        "数量对应本次已下载的 T1-T4 产物；重跑后需要重新生成。",
    )

    # SNP-heritability QC: retain all traits for QC display while distinguishing interpretable estimates.
    h2_fields, h2_rows = read_tsv(server / "ldsc_package" / "figure_sinonasal_ldsc_h2_qc.tsv")
    h2_reviewed: list[dict[str, object]] = []
    for row in h2_rows:
        h2 = parse_float(row.get("h2"))
        se = parse_float(row.get("h2_se"))
        z_approx = h2 / se if finite(h2) and finite(se) and se > 0 else math.nan
        if not finite(h2) or h2 <= 0:
            status = "invalid_h2_nonpositive"
            reason = "non-positive or non-finite h2"
            interpretable = False
        elif not finite(se) or se <= 0:
            status = "invalid_h2_se"
            reason = "non-positive or non-finite standard error"
            interpretable = False
        elif z_approx < 2:
            status = "valid_numeric_low_power"
            reason = "h2 estimate is positive but z_approx < 2"
            interpretable = False
        else:
            status = "valid_numeric"
            reason = ""
            interpretable = True
        outrow = dict(row)
        outrow.update(
            {
                "h2_z_approx": fmt(z_approx),
                "validation_status": status,
                "interpretation_caution": reason,
                "qc_display_include": "True",
                "association_interpretable": str(interpretable),
            }
        )
        h2_reviewed.append(outrow)
    h2_review_fields = h2_fields + [
        "h2_z_approx",
        "validation_status",
        "interpretation_caution",
        "qc_display_include",
        "association_interpretable",
    ]
    supp_fig1 = outdir / "supplementary_figures" / "SuppFig_S1_ldsc_h2_qc"
    h2_path = supp_fig1 / "figure_sinonasal_ldsc_h2_qc_reviewed.tsv"
    write_tsv(h2_path, h2_review_fields, h2_reviewed)
    register(
        h2_path,
        "Supplementary Figure S1",
        "reviewed_figure_input",
        "available_after_qc_review",
        str(server / "ldsc_package" / "figure_sinonasal_ldsc_h2_qc.tsv"),
        "SNP 遗传力与 QC 参数审阅版，保留全部表型用于展示，同时新增低功效或非正遗传力警示列。",
        "标记为 `association_interpretable=False` 的表型可以在 QC 图中展示，但不应据此声明遗传相关缺失或存在。",
    )

    # Supporting tables and execution/QC evidence.
    supp_panel = outdir / "supplementary_tables" / "SuppTable_S1_trait_panel"
    copy_registered(
        disease / "phase0_active_v3_panel.tsv",
        supp_panel / "phase0_active_v3_panel.tsv",
        "Supplementary Table S1",
        "study_design_table",
        "available",
        "记录纳入分析的疾病表型、FinnGen accession、样本量、病例/对照数与纳入依据。",
        "这是纳入设计表；最终结果可用性还受每个分析步骤 QC 影响。",
    )
    supp_pairs = outdir / "supplementary_tables" / "SuppTable_S2_pair_design"
    copy_registered(
        disease / "phase0_nasal_vs_other_pair_manifest.tsv",
        supp_pairs / "phase0_nasal_vs_other_pair_manifest.tsv",
        "Supplementary Table S2",
        "study_design_table",
        "available",
        "定义以鼻部表型为锚点的预设跨表型比较集合，明确多重检验族及结果表的配对来源。",
        "该表定义计划比较，不等同于所有配对都产生了可解释结果。",
    )
    supp_qc = outdir / "supplementary_tables" / "SuppTable_S4_processing_qc"
    copy_registered(
        server / "download_package" / "phase0_v3_download_verify.tsv",
        supp_qc / "phase0_v3_download_verify.tsv",
        "Supplementary Table S4",
        "qc_evidence_table",
        "available",
        "记录预设 GWAS 汇总统计文件的下载存在性与文件大小检查，可作为输入数据到位证据。",
        "检查证明文件存在，不替代内容语义或统计质量审查。",
    )
    copy_registered(
        server / "phase0_v3_munge_qc.tsv",
        supp_qc / "phase0_v3_munge_qc.tsv",
        "Supplementary Table S4",
        "qc_evidence_table",
        "available",
        "记录各表型清洗、去重、排除歧义变异/MHC 以及写出行数，用于补充报告数据预处理过程。",
        "这些计数对应当前服务器运行版本与输入文件；重跑或规则变化后应更新。",
    )

    # Provenance needed to understand or reproduce the manuscript-facing package.
    provenance = outdir / "provenance"
    provenance_files = [
        (metadata / "figure_catalog.tsv", "原有图件目录及其预期定位。"),
        (metadata / "figure_data_schema.tsv", "原有图件输入字段模式定义。"),
        (metadata / "phase0_main_figure_plan.tsv", "正文图件规划及内容范围。"),
        (server / "manifest_validation_report.tsv", "输入 manifest 验证报告。"),
        (server / "manifest_validation_summary.md", "输入 manifest 验证摘要。"),
        (contract, "各图件输出文件的字段合同。"),
    ]
    for source, explanation in provenance_files:
        copy_registered(
            source,
            provenance / source.name,
            "Provenance / audit trail",
            "provenance",
            "available",
            explanation,
            "用于复核和追溯，不作为直接生物学结果表解读。",
        )

    availability_rows = [
        {
            "manuscript_item": "Fig1",
            "component": "workflow design",
            "status": "available",
            "available_file": "main_figures/Fig1_workflow/figure_workflow_nodes.tsv",
            "reason_or_gap": "workflow metadata and rendering specification exist",
        },
        {
            "manuscript_item": "Fig2",
            "component": "global LDSC rg layer",
            "status": "available_after_qc_review",
            "available_file": "main_figures/Fig2_global_ldsc_available_layer/figure_global_rg_summary_reviewed.tsv",
            "reason_or_gap": "reviewed global rg output is available; invalid estimates excluded from interpretation",
        },
        {
            "manuscript_item": "Fig2",
            "component": "candidate variant / local evidence layer",
            "status": "not_yet_available",
            "available_file": "",
            "reason_or_gap": "requires downstream local-rg/MTAG/coloc completion",
        },
        {
            "manuscript_item": "SuppFig_S1",
            "component": "LDSC h2 QC",
            "status": "available_after_qc_review",
            "available_file": "supplementary_figures/SuppFig_S1_ldsc_h2_qc/figure_sinonasal_ldsc_h2_qc_reviewed.tsv",
            "reason_or_gap": "available with low-power/non-positive h2 annotations",
        },
        {
            "manuscript_item": "SuppFig_S2-S5",
            "component": "local/partitioned rg",
            "status": "not_yet_available",
            "available_file": "",
            "reason_or_gap": "T5 LAVA failed before output; corrected pipeline requires later rerun",
        },
        {
            "manuscript_item": "SuppFig_S6",
            "component": "MTAG Manhattan",
            "status": "not_yet_available",
            "available_file": "",
            "reason_or_gap": "downstream of incomplete local-rg stage",
        },
        {
            "manuscript_item": "SuppFig_S7-S10",
            "component": "sensitivity/replication/candidate locus panels",
            "status": "not_yet_available",
            "available_file": "",
            "reason_or_gap": "required downstream analyses have not completed",
        },
        {
            "manuscript_item": "Fig3 / SuppFig_S11-S15",
            "component": "functional interpretation",
            "status": "not_yet_available",
            "available_file": "",
            "reason_or_gap": "functional analysis outputs have not been generated",
        },
        {
            "manuscript_item": "Fig4 / SuppFig_S16-S18",
            "component": "intervention/MR evidence",
            "status": "not_yet_available",
            "available_file": "",
            "reason_or_gap": "intervention analysis outputs have not been generated",
        },
    ]
    availability = outdir / "availability_status.tsv"
    write_tsv(
        availability,
        ["manuscript_item", "component", "status", "available_file", "reason_or_gap"],
        availability_rows,
    )
    register(
        availability,
        "Package overview",
        "status_inventory",
        "available",
        "figure output contract plus current partial-run evidence",
        "逐项声明哪些正文或补充图件现在具备可用信息、哪些只能在后续重跑成功后补入文章。",
        "标记为 `not_yet_available` 的项目不得在当前稿件中作为已有结果描述。",
    )

    unavailable = outdir / "unavailable_planned_assets" / "planned_but_not_yet_available.tsv"
    unavailable_rows = [row for row in availability_rows if row["status"] == "not_yet_available"]
    write_tsv(
        unavailable,
        ["manuscript_item", "component", "status", "available_file", "reason_or_gap"],
        unavailable_rows,
    )
    register(
        unavailable,
        "Gap register",
        "status_inventory",
        "available",
        str(availability),
        "集中列出规划中但目前没有正式可用产物的图表模块，防止在写作阶段误将计划结果当作已完成结果。",
        "该清单是缺口管理材料，而非待补图件本身。",
    )

    summary = [
        {"metric": "planned_global_rg_pairs", "value": str(planned_n)},
        {"metric": "valid_global_rg_pairs", "value": str(len(valid_rg))},
        {"metric": "excluded_global_rg_pairs", "value": str(len(excluded_rg))},
        {"metric": "significant_global_rg_pairs_bh_q_lt_0.05", "value": str(len(significant_rg))},
        {"metric": "h2_traits_total", "value": str(len(h2_reviewed))},
        {
            "metric": "h2_traits_invalid_nonpositive",
            "value": str(sum(row["validation_status"] == "invalid_h2_nonpositive" for row in h2_reviewed)),
        },
        {
            "metric": "h2_traits_low_power_positive",
            "value": str(sum(row["validation_status"] == "valid_numeric_low_power" for row in h2_reviewed)),
        },
    ]
    metrics = outdir / "manuscript_current_result_counts.tsv"
    write_tsv(metrics, ["metric", "value"], summary)
    register(
        metrics,
        "Package overview",
        "qc_summary_table",
        "available_after_qc_review",
        "derived from reviewed LDSC tables in this package",
        "本次可用结果的核心计数摘要，适合写作时快速核对正文陈述与补充材料数量。",
        "计数仅涵盖当前已完成的全局 LDSC 层，不能外推到未完成的下游分析。",
    )

    manifest_path = outdir / "manuscript_asset_manifest.tsv"
    manifest_fields = [
        "asset_file",
        "description_file",
        "manuscript_role",
        "asset_type",
        "status",
        "source",
    ]
    write_tsv(manifest_path, manifest_fields, manifest)
    write_text(
        description_path(manifest_path),
        """# manuscript_asset_manifest.tsv

## 用途

该文件是当前文章可用材料的总目录。每一行指向一个实际数据/定义/审计文件及其专属说明文件，可用于写作移交和复核。

## 使用限制

目录只纳入本次明确存在的文件或由现有 LDSC 结果审阅派生的表；缺失的后续分析单独列在 `availability_status.tsv` 中。
""",
    )

    readme = f"""# Current Manuscript Assets: Phase0 Partial Run

本目录只收集截至 2026-05-25 已经真实产生、能够支撑文章正文或补充材料组织工作的文件。服务器 full run 在 LAVA 阶段失败，因此本包不会把 LAVA、MTAG、coloc 或功能/干预分析规划伪装成现有结果。

## 当前可用内容

- `main_figures/Fig1_workflow/`：研究流程图节点、连接和 Mermaid 定义。
- `main_figures/Fig2_global_ldsc_available_layer/`：Fig2 当前已完成的全局 LDSC 遗传相关层，包含经 QC 的全表和显著结果子集。
- `supplementary_figures/SuppFig_S1_ldsc_h2_qc/`：经有效性/低功效标注的 SNP 遗传力 QC 输入。
- `supplementary_tables/`：表型纳入设计、预设配对、LDSC 排除/QC、下载与清洗 QC 证据。
- `provenance/`：图件合同和 manifest 审计依据。

## 关键 QC 结论

- 全局 `rg` 计划配对：{planned_n}；QC 后可解释：{len(valid_rg)}；排除：{len(excluded_rg)}。
- 按全部 {planned_n} 个计划配对构成的检验族重新计算 BH 校正后，`q < 0.05` 的有效配对：{len(significant_rg)}。
- `h2` 表中非正值结果：{sum(row["validation_status"] == "invalid_h2_nonpositive" for row in h2_reviewed)}；正值但低功效警示结果：{sum(row["validation_status"] == "valid_numeric_low_power" for row in h2_reviewed)}。
- 原服务器 `figure_global_rg_summary.tsv` 的 `p_adj` 未填充，文章统计陈述应使用审阅版新增列 `p_adj_bh_planned_family`。

## 每项说明

每个纳入的 `.tsv` 或 `.mmd` 文件旁均有同名 `.DESCRIPTION.md`，说明其来源、文章用途和使用限制。总目录见 `manuscript_asset_manifest.tsv`；未产生图件的状态见 `availability_status.tsv`。

## 来源

- 已下载服务器部分运行产物：`{source_root}`
- 本地图件输出合同：`{contract}`
"""
    write_text(outdir / "README.md", readme)

    print(f"Created manuscript assets: {outdir}")
    print(f"Manifest assets: {len(manifest)}")
    print(f"Global rg planned={planned_n} valid={len(valid_rg)} excluded={len(excluded_rg)} significant={len(significant_rg)}")


if __name__ == "__main__":
    main()
