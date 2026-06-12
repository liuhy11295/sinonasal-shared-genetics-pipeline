#!/usr/bin/env python3
from __future__ import annotations

import csv
import gzip
import math
import os
import subprocess
import sys
import time
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import hypergeom
from statsmodels.stats.multitest import multipletests


project_root = os.environ.get("PROJECT_ROOT", "").strip()
if not project_root:
    raise SystemExit("Set PROJECT_ROOT to the analysis project root.")
ROOT = Path(project_root)
OUT = ROOT / "results_end/twas"
REF = OUT / "reference"
FUSION = REF / "fusion_twas/FUSION.assoc_test.R"
LDREF = REF / "LDREF/1000G.EUR."
WEIGHTS = REF / "weights"
conda_root = os.environ.get("CONDA_ROOT", "").strip()
conda_env = os.environ.get("PHASE0_CONDA_ENV", "").strip()
if not conda_root or not conda_env:
    raise SystemExit("Set CONDA_ROOT and PHASE0_CONDA_ENV before running the TWAS pipeline.")
CONDA_ROOT = Path(conda_root)
ENV = Path(conda_env)
R_BIN = ENV / "bin/Rscript"
PY_BIN = ENV / "bin/python"
TIME = "/usr/bin/time"

PRIMARY = ["Lung", "Whole_Blood"]
SECONDARY = [
    "Spleen",
    "Cells_EBV-transformed_lymphocytes",
    "Skin_Sun_Exposed_Lower_leg",
    "Skin_Not_Sun_Exposed_Suprapubic",
]
TISSUES = PRIMARY + SECONDARY


def open_text(path: Path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path, "rt")


def write_tsv(rows, path: Path, fields=None):
    path.parent.mkdir(parents=True, exist_ok=True)
    if fields is None:
        fields = sorted({k for row in rows for k in row})
    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, delimiter="\t", fieldnames=fields, extrasaction="ignore", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def fnum(x):
    try:
        if x is None or str(x).strip() == "":
            return math.nan
        return float(x)
    except Exception:
        return math.nan


def read_pairs():
    manifest = ROOT / "results/step8_mtag/evidence_pair_manifest.tsv"
    rows = pd.read_csv(manifest, sep="\t", dtype=str).fillna("")
    rows = rows.loc[~rows.astype(str).apply(lambda c: c.str.contains("CHRONIC_RHINITIS_PANUKB_J31", regex=False, na=False)).any(axis=1)]
    return rows.to_dict("records")


def full_mtag(pair_id: str) -> Path:
    return ROOT / "results/step8_mtag/full" / pair_id / "full_mtag.tsv.gz"


def inspect_and_convert(pair: dict) -> dict:
    pair_id = pair["pair_id"]
    src = full_mtag(pair_id)
    dest = OUT / "sumstats" / f"{pair_id}.fusion.sumstats"
    info = {
        "pair_id": pair_id,
        "trait1": pair.get("trait1", ""),
        "trait2": pair.get("trait2", ""),
        "full_mtag_file": str(src),
        "fusion_sumstats": str(dest),
        "status": "PASS",
        "reason": "",
        "n_input_rows": 0,
        "n_output_rows": 0,
        "has_SNP": False,
        "has_A1": False,
        "has_A2": False,
        "has_Z": False,
        "has_P": False,
        "has_N": False,
        "z_source": "",
    }
    if not src.exists():
        info.update(status="failed_missing_full_mtag", reason="full_mtag.tsv.gz not found")
        return info
    with open_text(src) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        cols = set(reader.fieldnames or [])
        info.update(
            has_SNP="SNP" in cols,
            has_A1="A1" in cols,
            has_A2="A2" in cols,
            has_Z=("Z" in cols) or ("mtag_z" in cols),
            has_P=("P" in cols) or ("mtag_pval" in cols),
            has_N="N" in cols,
        )
        has_beta_se = ("BETA" in cols and "SE" in cols) or ("mtag_beta" in cols and "mtag_se" in cols)
        if not (info["has_SNP"] and info["has_A1"] and info["has_A2"] and info["has_P"] and info["has_N"]):
            info.update(status="failed_missing_required_fields", reason="missing one of SNP,A1,A2,P,N")
            return info
        if not info["has_Z"] and not has_beta_se:
            info.update(status="failed_missing_z", reason="missing Z and cannot compute from BETA/SE")
            return info
        if dest.exists() and dest.stat().st_size > 0:
            with open(dest) as existing:
                info["n_output_rows"] = max(sum(1 for _ in existing) - 1, 0)
            info["n_input_rows"] = ""
            info["z_source"] = "existing_converted_sumstats"
            return info
        dest.parent.mkdir(parents=True, exist_ok=True)
        with open(dest, "w", newline="") as out:
            writer = csv.DictWriter(out, delimiter="\t", fieldnames=["SNP", "A1", "A2", "Z", "P", "N"], lineterminator="\n")
            writer.writeheader()
            for row in reader:
                info["n_input_rows"] += 1
                snp = row.get("SNP")
                a1 = row.get("A1")
                a2 = row.get("A2")
                p = row.get("mtag_pval") or row.get("P")
                n = row.get("N")
                z = row.get("mtag_z") or row.get("Z")
                z_source = "existing_Z"
                if not z:
                    beta = fnum(row.get("mtag_beta") or row.get("BETA"))
                    se = fnum(row.get("mtag_se") or row.get("SE"))
                    if math.isfinite(beta) and math.isfinite(se) and se > 0:
                        z = str(beta / se)
                        z_source = "BETA_over_SE"
                if not all([snp, a1, a2, z, p, n]):
                    continue
                writer.writerow({"SNP": snp, "A1": a1, "A2": a2, "Z": z, "P": p, "N": n})
                info["n_output_rows"] += 1
                info["z_source"] = info["z_source"] or z_source
    if info["n_output_rows"] == 0:
        info.update(status="failed_no_valid_rows", reason="no rows after field filtering")
    return info


def input_qc():
    rows = [inspect_and_convert(pair) for pair in read_pairs()]
    write_tsv(rows, OUT / "twas_input_qc.tsv")
    return rows


def benchmark_pair() -> str:
    qc_path = OUT / "twas_input_qc.tsv"
    if not qc_path.exists():
        input_qc()
    qc = pd.read_csv(qc_path, sep="\t", dtype=str).fillna("")
    ok = qc.loc[qc["status"] == "PASS"].copy()
    ok["n_output_rows_num"] = pd.to_numeric(ok["n_output_rows"], errors="coerce").fillna(0)
    if ok.empty:
        raise SystemExit("No PASS full MTAG inputs for TWAS benchmark")
    return ok.sort_values("n_output_rows_num", ascending=False).iloc[0]["pair_id"]


def run_cmd(cmd, log: Path, cwd: Path | None = None):
    log.parent.mkdir(parents=True, exist_ok=True)
    with open(log, "w") as fh:
        fh.write("COMMAND\t" + " ".join(map(str, cmd)) + "\n")
        fh.flush()
        start = time.time()
        proc = subprocess.run(cmd, cwd=str(cwd) if cwd else None, stdout=fh, stderr=subprocess.STDOUT, text=True)
        fh.write(f"\nRETURN_CODE\t{proc.returncode}\nWALLTIME_SECONDS\t{time.time() - start:.3f}\n")
    return proc.returncode


def run_one(pair_id: str, tissue: str, benchmark: bool = False) -> dict:
    pair_lookup = {p["pair_id"]: p for p in read_pairs()}
    pair = pair_lookup[pair_id]
    sumstats = OUT / "sumstats" / f"{pair_id}.fusion.sumstats"
    if not sumstats.exists():
        inspect_and_convert(pair)
    pair_dir = OUT / "per_pair" / pair_id
    pair_dir.mkdir(parents=True, exist_ok=True)
    out_file = pair_dir / f"{tissue}.fusion.tsv"
    log = pair_dir / f"{tissue}.log"
    time_log = pair_dir / f"{tissue}.time.txt"
    for stale in [out_file, log, time_log]:
        if stale.exists():
            stale.unlink()
    weights = WEIGHTS / f"{tissue}.pos"
    if not weights.exists():
        return {"pair_id": pair_id, "tissue": tissue, "status": "failed_missing_weights", "peak_memory_kb": "", "out": str(out_file)}
    cmd = [
        TIME, "-v", "-o", str(time_log),
        str(R_BIN), str(FUSION),
        "--sumstats", str(sumstats),
        "--weights", str(weights),
        "--weights_dir", str(WEIGHTS),
        "--ref_ld_chr", str(LDREF),
        "--chr", "1",
        "--out", str(out_file),
    ]
    # FUSION handles one chromosome at a time; run all chromosomes and concatenate.
    parts = []
    statuses = []
    total_peak = 0
    for chrom in range(1, 23):
        chr_out = pair_dir / f"{tissue}.chr{chrom}.fusion.tsv"
        chr_log = pair_dir / f"{tissue}.chr{chrom}.log"
        chr_time = pair_dir / f"{tissue}.chr{chrom}.time.txt"
        for stale in [chr_out, chr_log, chr_time]:
            if stale.exists():
                stale.unlink()
        chr_cmd = [
            TIME, "-v", "-o", str(chr_time),
            str(R_BIN), str(FUSION),
            "--sumstats", str(sumstats),
            "--weights", str(weights),
            "--weights_dir", str(WEIGHTS),
            "--ref_ld_chr", str(LDREF),
            "--chr", str(chrom),
            "--out", str(chr_out),
        ]
        rc = run_cmd(chr_cmd, chr_log, cwd=FUSION.parent)
        statuses.append(rc)
        if chr_out.exists() and chr_out.stat().st_size > 0:
            parts.append(chr_out)
        peak = parse_peak(chr_time)
        if peak and peak > total_peak:
            total_peak = peak
    concat_fusion(parts, out_file)
    return {"pair_id": pair_id, "tissue": tissue, "status": "PASS" if all(s == 0 for s in statuses) else "failed_some_chrom", "peak_memory_kb": total_peak, "out": str(out_file)}


def parse_peak(time_file: Path):
    if not time_file.exists():
        return 0
    for line in time_file.read_text(errors="ignore").splitlines():
        if "Maximum resident set size" in line:
            try:
                return int(line.split(":")[-1].strip())
            except Exception:
                return 0
    return 0


def concat_fusion(parts, out_file: Path):
    out_file.parent.mkdir(parents=True, exist_ok=True)
    wrote = False
    with open(out_file, "w", newline="") as out:
        for part in parts:
            with open(part) as fh:
                for i, line in enumerate(fh):
                    if i == 0 and wrote:
                        continue
                    out.write(line)
                wrote = True
    if not wrote:
        out_file.write_text("")


def run_benchmark():
    if not (OUT / "twas_input_qc.tsv").exists():
        input_qc()
    pair_id = benchmark_pair()
    rows = []
    for tissue in TISSUES:
        rows.append(run_one(pair_id, tissue, benchmark=True))
    fields = ["pair_id", "tissue", "status", "peak_memory_kb", "out"]
    write_tsv(rows, OUT / "resource_benchmark.tsv", fields)
    peak = max([int(r["peak_memory_kb"] or 0) for r in rows] or [0])
    safe_gb = peak / 1024 / 1024 * 1.5
    jobs = max(1, min(40, int(512 // safe_gb) if safe_gb > 0 else 1))
    write_tsv([{"peak_memory_kb": peak, "safe_memory_gb": safe_gb, "final_parallel_jobs": jobs}], OUT / "parallelization_plan.tsv")


def run_array_task():
    task = int(os.environ.get("SLURM_ARRAY_TASK_ID", "1"))
    pairs = [r for r in pd.read_csv(OUT / "twas_input_qc.tsv", sep="\t", dtype=str).fillna("").to_dict("records") if r["status"] == "PASS"]
    jobs = [(p["pair_id"], t) for p in pairs for t in TISSUES]
    pair_id, tissue = jobs[task - 1]
    result = run_one(pair_id, tissue)
    write_tsv([result], OUT / "logs" / f"array_{task}.status.tsv")


def aggregate_twas():
    pairs = {p["pair_id"]: p for p in read_pairs()}
    all_rows = []
    qc_rows = []
    for pair_id, pair in pairs.items():
        for tissue in TISSUES:
            path = OUT / "per_pair" / pair_id / f"{tissue}.fusion.tsv"
            if not path.exists() or path.stat().st_size == 0:
                qc_rows.append({"pair_id": pair_id, "tissue": tissue, "status": "missing_output", "n_genes": 0})
                continue
            df = pd.read_csv(path, sep="\t", dtype=str)
            if df.empty:
                qc_rows.append({"pair_id": pair_id, "tissue": tissue, "status": "empty_output", "n_genes": 0})
                continue
            p = pd.to_numeric(df["TWAS.P"], errors="coerce")
            fdr = np.full(len(df), np.nan)
            ok = p.notna()
            if ok.any():
                fdr[ok] = multipletests(p[ok], method="fdr_bh")[1]
            n_tested = int(ok.sum())
            bonf_thr = 0.05 / n_tested if n_tested else np.nan
            df["pair_id"] = pair_id
            df["trait1"] = pair.get("trait1", "")
            df["trait2"] = pair.get("trait2", "")
            df["tissue"] = tissue
            df["gene_symbol"] = df.get("ID", "")
            df["FDR"] = fdr
            df["n_genes_tested_per_tissue"] = n_tested
            df["bonferroni_threshold"] = bonf_thr
            df["bonferroni_significant"] = p <= bonf_thr
            df["fdr_significant"] = df["FDR"].astype(float) <= 0.05
            keep = ["pair_id", "trait1", "trait2", "tissue", "gene_symbol", "TWAS.Z", "TWAS.P", "FDR", "NSNP", "MODEL", "bonferroni_significant", "fdr_significant", "n_genes_tested_per_tissue"]
            all_rows.append(df[[c for c in keep if c in df.columns]])
            qc_rows.append({"pair_id": pair_id, "tissue": tissue, "status": "PASS", "n_genes": n_tested})
    all_df = pd.concat(all_rows, ignore_index=True) if all_rows else pd.DataFrame(columns=["pair_id", "trait1", "trait2", "tissue", "gene_symbol", "TWAS.Z", "TWAS.P", "FDR", "NSNP", "MODEL"])
    all_df.to_csv(OUT / "twas_all_results.tsv", sep="\t", index=False)
    bonf = all_df.loc[all_df.get("bonferroni_significant", False).astype(str).isin(["True", "true", "1"])].copy()
    fdr = all_df.loc[all_df.get("fdr_significant", False).astype(str).isin(["True", "true", "1"])].copy()
    bonf.to_csv(OUT / "twas_significant_bonferroni.tsv", sep="\t", index=False)
    fdr.to_csv(OUT / "twas_significant_fdr.tsv", sep="\t", index=False)
    top = all_df.assign(_p=pd.to_numeric(all_df.get("TWAS.P"), errors="coerce")).sort_values(["pair_id", "tissue", "_p"]).groupby(["pair_id", "tissue"]).head(20).drop(columns=["_p"])
    top.to_csv(OUT / "top_twas_genes.tsv", sep="\t", index=False)
    write_tsv(qc_rows, OUT / "twas_qc_summary.tsv")
    overlap_and_ora(all_df, bonf, fdr)


def overlap_and_ora(all_df, bonf, fdr):
    magma = pd.read_csv(ROOT / "results/magma/gene_results.tsv", sep="\t", dtype=str).fillna("")
    magma["magma_p_num"] = pd.to_numeric(magma["p"], errors="coerce")
    magma["magma_bonf_p_num"] = pd.to_numeric(magma["gene_bonf_p"], errors="coerce")
    magma_sig = magma.loc[magma["magma_bonf_p_num"] <= 0.05].copy()
    make_overlap(magma_sig, bonf.loc[bonf["tissue"].isin(PRIMARY)].copy(), OUT / "magma_twas_overlap_bonferroni.tsv")
    make_overlap(magma_sig, fdr.copy(), OUT / "magma_twas_overlap_fdr.tsv")
    run_overlap_ora()


def make_overlap(magma_sig, twas_sig, path: Path):
    if twas_sig.empty:
        write_tsv([], path, ["pair_id", "trait1", "trait2", "gene_symbol", "magma_p", "magma_bonf_p", "best_twas_tissue", "best_twas_p", "best_twas_fdr", "n_supporting_tissues", "supporting_tissues"])
        return
    twas_sig = twas_sig.copy()
    twas_sig["TWAS.P.num"] = pd.to_numeric(twas_sig["TWAS.P"], errors="coerce")
    twas_sig["FDR.num"] = pd.to_numeric(twas_sig["FDR"], errors="coerce")
    rows = []
    for (pair, gene), tw in twas_sig.groupby(["pair_id", "gene_symbol"]):
        mg = magma_sig.loc[(magma_sig["pair_id"] == pair) & (magma_sig["gene_symbol"] == gene)]
        if mg.empty:
            continue
        best = tw.sort_values("TWAS.P.num").iloc[0]
        mrow = mg.sort_values("magma_bonf_p_num").iloc[0]
        tissues = sorted(tw["tissue"].dropna().unique())
        rows.append({"pair_id": pair, "trait1": mrow["trait1"], "trait2": mrow["trait2"], "gene_symbol": gene, "magma_p": mrow["p"], "magma_bonf_p": mrow["gene_bonf_p"], "best_twas_tissue": best["tissue"], "best_twas_p": best["TWAS.P"], "best_twas_fdr": best["FDR"], "n_supporting_tissues": len(tissues), "supporting_tissues": ",".join(tissues)})
    write_tsv(rows, path, ["pair_id", "trait1", "trait2", "gene_symbol", "magma_p", "magma_bonf_p", "best_twas_tissue", "best_twas_p", "best_twas_fdr", "n_supporting_tissues", "supporting_tissues"])


def load_gmt(path: Path):
    sets = []
    if not path.exists():
        return sets
    with open(path) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 3:
                sets.append((parts[0], parts[1], set(parts[2:])))
    return sets


def ora_for_db(fg_by_pair, bg_by_pair, sets, out_path):
    rows = []
    for pair, fg in fg_by_pair.items():
        bg = bg_by_pair.get(pair, set())
        if not fg or not bg:
            continue
        M = len(bg)
        N = len(fg & bg)
        for sid, desc, genes in sets:
            pathway = genes & bg
            n = len(pathway)
            k = len(fg & pathway)
            if n < 5 or k == 0:
                continue
            p = hypergeom.sf(k - 1, M, n, N)
            rows.append({"pair_id": pair, "pathway_id": sid, "pathway_description": desc, "overlap_genes": ",".join(sorted(fg & pathway)), "k": k, "set_size": n, "foreground_size": N, "background_size": M, "pvalue": p})
    columns = ["pair_id", "pathway_id", "pathway_description", "overlap_genes", "k", "set_size", "foreground_size", "background_size", "pvalue", "p.adjust"]
    df = pd.DataFrame(rows)
    if not df.empty:
        df["p.adjust"] = multipletests(pd.to_numeric(df["pvalue"], errors="coerce"), method="fdr_bh")[1]
        df = df.sort_values(["p.adjust", "pvalue"])
    else:
        df = pd.DataFrame(columns=columns)
    df.to_csv(out_path, sep="\t", index=False)
    return df


def run_overlap_ora():
    overlap = pd.read_csv(OUT / "magma_twas_overlap_bonferroni.tsv", sep="\t", dtype=str).fillna("")
    magma = pd.read_csv(ROOT / "results/magma/gene_results.tsv", sep="\t", dtype=str).fillna("")
    fg_by_pair = defaultdict(set)
    for _, r in overlap.iterrows():
        fg_by_pair[r["pair_id"]].add(r["gene_symbol"])
    bg_by_pair = defaultdict(set)
    for _, r in magma.iterrows():
        bg_by_pair[r["pair_id"]].add(r["gene_symbol"])
    gmt = ROOT / "results/phase0_extension/step8_magma_gene_pathway/resources/gmt"
    ora_dir = OUT / "ora"
    ora_dir.mkdir(exist_ok=True)
    dbs = {
        "go_bp": gmt / "c5.go.bp.v2023.2.Hs.symbols.gmt",
        "kegg": gmt / "c2.cp.kegg.v2023.2.Hs.symbols.gmt",
        "reactome": gmt / "c2.cp.reactome.v2023.2.Hs.symbols.gmt",
    }
    summaries = []
    for name, path in dbs.items():
        df = ora_for_db(fg_by_pair, bg_by_pair, load_gmt(path), ora_dir / f"{name}.tsv")
        summaries.append({"database": name, "n_rows": len(df), "n_sig_fdr05": int((pd.to_numeric(df.get("p.adjust"), errors="coerce") <= 0.05).sum()) if not df.empty else 0})
    write_tsv(summaries, ora_dir / "pathway_summary.tsv")
    all_ora = []
    for name in dbs:
        p = ora_dir / f"{name}.tsv"
        if p.exists() and p.stat().st_size:
            try:
                df = pd.read_csv(p, sep="\t", dtype=str)
            except pd.errors.EmptyDataError:
                df = pd.DataFrame()
            if not df.empty:
                df["database"] = name
                all_ora.append(df)
    if all_ora:
        all_df = pd.concat(all_ora, ignore_index=True)
        rec = all_df.loc[pd.to_numeric(all_df.get("p.adjust"), errors="coerce") <= 0.05].groupby(["database", "pathway_id", "pathway_description"], as_index=False).agg(n_pairs=("pair_id", "nunique"))
        rec.sort_values(["n_pairs", "database"], ascending=[False, True]).to_csv(ora_dir / "recurrent_pathways.tsv", sep="\t", index=False)
    else:
        write_tsv([], ora_dir / "recurrent_pathways.tsv", ["database", "pathway_id", "pathway_description", "n_pairs"])
    write_notes()


def write_notes():
    qc = pd.read_csv(OUT / "twas_qc_summary.tsv", sep="\t", dtype=str) if (OUT / "twas_qc_summary.tsv").exists() else pd.DataFrame()
    twas = pd.read_csv(OUT / "twas_all_results.tsv", sep="\t", dtype=str) if (OUT / "twas_all_results.tsv").exists() else pd.DataFrame()
    bonf = pd.read_csv(OUT / "twas_significant_bonferroni.tsv", sep="\t", dtype=str) if (OUT / "twas_significant_bonferroni.tsv").exists() else pd.DataFrame()
    ov = pd.read_csv(OUT / "magma_twas_overlap_bonferroni.tsv", sep="\t", dtype=str) if (OUT / "magma_twas_overlap_bonferroni.tsv").exists() else pd.DataFrame()
    readme = OUT / "README.md"
    readme.write_text(
        "# Step 9 TWAS results\n\n"
        "All Step 9 outputs are stored under `results_end/twas/`.\n\n"
        f"Tissues: primary={', '.join(PRIMARY)}; secondary={', '.join(SECONDARY)}.\n\n"
        f"Successful pair-tissue TWAS outputs: {int((qc.get('status') == 'PASS').sum()) if not qc.empty else 0}.\n\n"
        f"TWAS result rows: {len(twas)}.\n\n"
        f"Bonferroni-significant TWAS rows: {len(bonf)}.\n\n"
        f"MAGMA-TWAS Bonferroni overlap rows: {len(ov)}.\n",
        encoding="utf-8",
    )
    (OUT / "method_notes.md").write_text(
        "FUSION TWAS was run using full MTAG cross-trait summary statistics. "
        "Inputs were checked for SNP, A1, A2, Z, P, N. If needed, Z can be computed as BETA/SE. "
        "MAGMA overlap used existing MAGMA gene results only; MAGMA and GSEA were not rerun. "
        "ORA background was the MAGMA tested gene set per pair, not all human genes.\n",
        encoding="utf-8",
    )
    with open(OUT / "update_log.md", "a", encoding="utf-8") as fh:
        fh.write(f"- {time.strftime('%F %T')} Step9 aggregation/notes updated.\\n")


def main():
    cmd = sys.argv[1] if len(sys.argv) > 1 else "input_qc"
    if cmd == "input_qc":
        input_qc()
    elif cmd == "benchmark":
        run_benchmark()
    elif cmd == "array_task":
        run_array_task()
    elif cmd == "aggregate":
        aggregate_twas()
    elif cmd == "benchmark_pair":
        print(benchmark_pair())
    else:
        raise SystemExit(f"unknown command: {cmd}")


if __name__ == "__main__":
    main()
