#!/usr/bin/env python3
from __future__ import annotations

import csv
import gzip
import math
import os
import random
import shutil
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from joblib import Parallel, delayed
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests


ROOT = Path(os.environ.get("PROJECT_ROOT", "/platform_data/p_user/p010/phase0"))
PAIR_ID = os.environ.get("PAIR_ID", "").strip()
STEP8_MTAG = ROOT / "results/step8_mtag"
STEP8_MAGMA = ROOT / "results/step8_magma"
STEP8_PATHWAY = ROOT / "results/step8_pathway"
BENCH = ROOT / "results/resource_benchmark"
MANIFEST = STEP8_MTAG / "evidence_pair_manifest.tsv"
MTAG_PY = ROOT / "tools/mtag/mtag.py"
PY2 = Path(os.environ.get("PHASE0_PY2", "/platform_data/p_user/p010/conda_envs/ldsc_mtag_py2/bin/python"))
MAGMA = ROOT / "results/phase0_extension/step8_magma_gene_pathway/resources/magma_official/magma"
GENELOC = ROOT / "results/phase0_extension/step8_magma_gene_pathway/resources/magma_official/NCBI37.3.gene.loc"
BFILE = ROOT / "results/phase0_extension/step8_magma_gene_pathway/resources/plink_merged/1000G.EUR.QC"
GMT_DIR = ROOT / "results/phase0_extension/step8_magma_gene_pathway/resources/gmt"


def open_text(path: Path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path, "rt")


def write_tsv(rows, path: Path, fields):
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, delimiter="\t", fieldnames=fields, extrasaction="ignore")
        w.writeheader()
        for row in rows:
            w.writerow(row)


def read_manifest_row() -> dict:
    with open(MANIFEST, newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    if not rows:
        raise SystemExit(f"empty manifest: {MANIFEST}")
    if PAIR_ID:
        for row in rows:
            if row["pair_id"] == PAIR_ID:
                return row
        raise SystemExit(f"PAIR_ID not found in manifest: {PAIR_ID}")
    task_id = os.environ.get("SLURM_ARRAY_TASK_ID", "").strip()
    if task_id:
        idx = int(task_id) - 1
        if idx < 0 or idx >= len(rows):
            raise SystemExit(f"SLURM_ARRAY_TASK_ID {task_id} out of range 1..{len(rows)}")
        return rows[idx]
    rows.sort(key=lambda r: int(float(r.get("n_snps_estimated") or 0)), reverse=True)
    return rows[0]


def float_or_nan(x):
    try:
        if x is None or str(x).strip() == "":
            return math.nan
        return float(x)
    except Exception:
        return math.nan


def standard_to_mtag_input(src: Path, dest: Path) -> int:
    dest.parent.mkdir(parents=True, exist_ok=True)
    n = 0
    with open_text(src) as fh, gzip.open(dest, "wt", newline="") as out:
        reader = csv.DictReader(fh, delimiter="\t")
        fields = ["SNP", "CHR", "BP", "A1", "A2", "Z", "P", "N", "FRQ"]
        writer = csv.DictWriter(out, delimiter="\t", fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        for row in reader:
            snp = row.get("SNP") or row.get("rsid")
            chr_ = row.get("CHR") or row.get("chrom") or row.get("chromosome")
            bp = row.get("BP") or row.get("POS") or row.get("position")
            a1 = row.get("A1") or row.get("EA") or row.get("ALT")
            a2 = row.get("A2") or row.get("OA") or row.get("REF")
            p = row.get("P") or row.get("p") or row.get("PVAL")
            nval = row.get("N") or row.get("n")
            frq = row.get("FRQ") or row.get("EAF") or row.get("eaf")
            z = row.get("Z") or row.get("z")
            if not z:
                beta = float_or_nan(row.get("BETA"))
                se = float_or_nan(row.get("SE"))
                if math.isfinite(beta) and math.isfinite(se) and se > 0:
                    z = str(beta / se)
            if not all([snp, chr_, bp, a1, a2, z, p, nval]):
                continue
            writer.writerow({"SNP": snp, "CHR": chr_, "BP": bp, "A1": a1, "A2": a2, "Z": z, "P": p, "N": nval, "FRQ": frq or "NA"})
            n += 1
    return n


def run(cmd, log: Path, cwd: Path | None = None):
    log.parent.mkdir(parents=True, exist_ok=True)
    with open(log, "w") as fh:
        fh.write("COMMAND\t" + " ".join(map(str, cmd)) + "\n")
        fh.flush()
        start = time.time()
        proc = subprocess.run(cmd, cwd=str(cwd) if cwd else None, stdout=fh, stderr=subprocess.STDOUT, text=True)
        fh.write(f"\nRETURN_CODE\t{proc.returncode}\nWALLTIME_SECONDS\t{time.time() - start:.3f}\n")
    if proc.returncode != 0:
        raise RuntimeError(f"command failed ({proc.returncode}): {' '.join(map(str, cmd))}")


def gzip_copy(src: Path, dest: Path):
    dest.parent.mkdir(parents=True, exist_ok=True)
    with open(src, "rt") as fh, gzip.open(dest, "wt") as out:
        shutil.copyfileobj(fh, out)


def summarize_full_mtag(full_gz: Path, qc_path: Path) -> dict:
    pvals = []
    chrs = set()
    n = 0
    with gzip.open(full_gz, "rt") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            n += 1
            if row.get("CHR"):
                chrs.add(str(row["CHR"]).replace("chr", ""))
            p = float_or_nan(row.get("mtag_pval") or row.get("P"))
            if math.isfinite(p):
                pvals.append(p)
    arr = np.array(pvals, dtype=float)
    summary = {
        "n_snps": n,
        "min_p": np.nanmin(arr) if arr.size else math.nan,
        "max_p": np.nanmax(arr) if arr.size else math.nan,
        "median_p": np.nanmedian(arr) if arr.size else math.nan,
        "prop_p_lt_5e8": float(np.mean(arr < 5e-8)) if arr.size else math.nan,
        "prop_p_lt_1e5": float(np.mean(arr < 1e-5)) if arr.size else math.nan,
        "n_chr": len(chrs),
    }
    write_tsv([summary], qc_path, list(summary.keys()))
    return summary


def make_magma_inputs(full_gz: Path, pair_id: str):
    input_dir = STEP8_MAGMA / "input"
    input_dir.mkdir(parents=True, exist_ok=True)
    pfile = input_dir / f"{pair_id}.magma_input.tsv"
    locfile = input_dir / f"{pair_id}.snps_loc.tsv"
    seen = set()
    with gzip.open(full_gz, "rt") as fh, open(pfile, "w", newline="") as pfh, open(locfile, "w", newline="") as lfh:
        reader = csv.DictReader(fh, delimiter="\t")
        pw = csv.DictWriter(pfh, delimiter="\t", fieldnames=["SNP", "P", "N"], lineterminator="\n")
        pw.writeheader()
        for row in reader:
            snp = row.get("SNP")
            if not snp or snp in seen:
                continue
            p = row.get("mtag_pval") or row.get("P")
            n = row.get("N") or "50000"
            chr_ = row.get("CHR")
            bp = row.get("BP")
            if not p or not chr_ or not bp:
                continue
            seen.add(snp)
            pw.writerow({"SNP": snp, "P": p, "N": n})
            lfh.write(f"{snp}\t{chr_}\t{bp}\n")
    return pfile, locfile


def run_magma(pair_id: str, pfile: Path, locfile: Path):
    out_dir = STEP8_MAGMA / "per_pair"
    out_dir.mkdir(parents=True, exist_ok=True)
    annot_prefix = out_dir / f"{pair_id}.annot"
    out_prefix = out_dir / pair_id
    run([str(MAGMA), "--annotate", "window=5,5", "--snp-loc", str(locfile), "--gene-loc", str(GENELOC), "--out", str(annot_prefix)], out_dir / f"{pair_id}.annot.wrapper.log")
    gene_annot = Path(str(annot_prefix) + ".genes.annot")
    run([str(MAGMA), "--bfile", str(BFILE), "--pval", str(pfile), "ncol=N", "--gene-annot", str(gene_annot), "--out", str(out_prefix)], out_dir / f"{pair_id}.magma.wrapper.log")
    return Path(str(out_prefix) + ".genes.out")


def gene_symbol_map():
    out = {}
    with open(GENELOC) as fh:
        for line in fh:
            parts = line.strip().split()
            if len(parts) >= 6:
                out[str(parts[0])] = parts[5]
    return out


def normalize_gene_results(pair_row: dict, genes_out: Path) -> pd.DataFrame:
    df = pd.read_csv(genes_out, sep=r"\s+", dtype=str)
    mapper = gene_symbol_map()
    p = pd.to_numeric(df["P"], errors="coerce")
    n_tested = int(p.notna().sum())
    out = pd.DataFrame(
        {
            "pair_id": pair_row["pair_id"],
            "trait1": pair_row["trait1"],
            "trait2": pair_row["trait2"],
            "nasal_trait": pair_row["nasal_trait"],
            "partner_trait": pair_row["partner_trait"],
            "gene_id": df["GENE"].astype(str),
            "gene_symbol": df["GENE"].astype(str).map(mapper).fillna(df["GENE"].astype(str)),
            "chr": df.get("CHR", ""),
            "start": df.get("START", ""),
            "stop": df.get("STOP", ""),
            "nsnps": df.get("NSNPS", ""),
            "nparam": df.get("NPARAM", ""),
            "z": df.get("ZSTAT", ""),
            "p": df.get("P", ""),
            "n_tested_genes": n_tested,
        }
    )
    out["gene_bonf_p"] = (p * n_tested).clip(upper=1.0)
    out["significant_gene"] = out["gene_bonf_p"] <= 0.05
    out["analysis_type"] = "standard_full_mtag_magma"
    out["magma_input_source"] = str(STEP8_MAGMA / "input" / f"{pair_row['pair_id']}.magma_input.tsv")
    dest = STEP8_MAGMA / "per_pair" / f"{pair_row['pair_id']}.gene_results.tsv"
    out.to_csv(dest, sep="\t", index=False)
    return out


def parse_gmt(path: Path):
    with open(path, errors="ignore") as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 3:
                yield parts[0], parts[1], [g for g in parts[2:] if g]


def enrichment_score(ranks: dict[str, float], genes: set[str]):
    ordered = sorted(ranks.items(), key=lambda x: x[1], reverse=True)
    hits = np.array([g in genes for g, _ in ordered], dtype=bool)
    nh = hits.sum()
    n = len(hits)
    if nh == 0 or nh == n:
        return 0.0, []
    weights = np.array([abs(v) for _, v in ordered], dtype=float)
    hit_weights = hits * weights
    hit_norm = hit_weights.sum()
    miss_norm = n - nh
    running = np.cumsum(np.where(hits, hit_weights / hit_norm, -1.0 / miss_norm))
    if abs(running.max()) >= abs(running.min()):
        idx = int(running.argmax())
        core = [g for g, h in zip([x[0] for x in ordered[: idx + 1]], hits[: idx + 1]) if h]
        return float(running[idx]), core
    idx = int(running.argmin())
    core = [g for g, h in zip([x[0] for x in ordered[idx:]], hits[idx:]) if h]
    return float(running[idx]), core


def run_gsea_for_db(pair_row: dict, gene_df: pd.DataFrame, db_name: str, gmt_name: str, n_perm: int = 100):
    gmt = GMT_DIR / gmt_name
    fields = ["pair_id", "trait1", "trait2", "database", "pathway_id", "pathway_description", "setSize", "NES", "pvalue", "p.adjust", "qvalue", "core_enrichment"]
    if not gmt.exists():
        return pd.DataFrame(columns=fields)
    vals = pd.to_numeric(gene_df["z"], errors="coerce")
    if vals.notna().sum() == 0:
        vals = -np.log10(pd.to_numeric(gene_df["p"], errors="coerce").clip(lower=1e-300))
    ranks = dict(zip(gene_df["gene_symbol"].astype(str), vals.fillna(0.0)))
    universe = list(ranks.keys())
    rows = []
    rng = random.Random(20260531)
    for pid, desc, genes in parse_gmt(gmt):
        gs = set(genes) & set(universe)
        if len(gs) < 10 or len(gs) > 500:
            continue
        es, core = enrichment_score(ranks, gs)
        null = []
        size = len(gs)
        for _ in range(n_perm):
            perm = set(rng.sample(universe, size))
            null_es, _ = enrichment_score(ranks, perm)
            null.append(null_es)
        null = np.array(null, dtype=float)
        if es >= 0:
            denom = np.mean(null[null >= 0]) if np.any(null >= 0) else 1.0
            pval = (np.sum(null >= es) + 1) / (len(null) + 1)
        else:
            denom = abs(np.mean(null[null < 0])) if np.any(null < 0) else 1.0
            pval = (np.sum(null <= es) + 1) / (len(null) + 1)
        nes = es / denom if denom else 0.0
        rows.append(
            {
                "pair_id": pair_row["pair_id"],
                "trait1": pair_row["trait1"],
                "trait2": pair_row["trait2"],
                "database": db_name,
                "pathway_id": pid,
                "pathway_description": desc,
                "setSize": len(gs),
                "NES": nes,
                "pvalue": pval,
                "core_enrichment": "/".join(core[:200]),
            }
        )
    df = pd.DataFrame(rows, columns=fields[:-2] + ["core_enrichment"])
    if df.empty:
        return pd.DataFrame(columns=fields)
    q = multipletests(pd.to_numeric(df["pvalue"], errors="coerce").fillna(1.0), method="fdr_bh")[1]
    df["p.adjust"] = q
    df["qvalue"] = q
    return df[fields]


def run_pathways(pair_row: dict, gene_df: pd.DataFrame):
    STEP8_PATHWAY.mkdir(parents=True, exist_ok=True)
    per_pair_dir = STEP8_PATHWAY / "per_pair"
    per_pair_dir.mkdir(parents=True, exist_ok=True)
    outputs = [
        ("GO_BP", "GO_Biological_Process_2023.gmt", "pair_specific_gsea_go_bp.tsv"),
        ("KEGG", "KEGG_2021_Human.gmt", "pair_specific_gsea_kegg.tsv"),
        ("Reactome", "Reactome_2022.gmt", "pair_specific_gsea_reactome.tsv"),
    ]
    for db, gmt, outname in outputs:
        df = run_gsea_for_db(pair_row, gene_df, db, gmt)
        df.to_csv(per_pair_dir / f"{pair_row['pair_id']}.{outname}", sep="\t", index=False)
        if os.environ.get("STEP8_WRITE_GLOBAL_PATHWAY", "0") == "1":
            df.to_csv(STEP8_PATHWAY / outname, sep="\t", index=False)


def gsea_one_set(pair_row: dict, ranks: dict[str, float], universe: list[str], pid: str, desc: str, genes: list[str], db_name: str, n_perm: int):
    gs = set(genes) & set(universe)
    if len(gs) < 10 or len(gs) > 500:
        return None
    set_vals = np.array([ranks[g] for g in gs], dtype=float)
    bg_vals = np.array([v for g, v in ranks.items() if g not in gs], dtype=float)
    if not len(set_vals) or not len(bg_vals):
        return None
    try:
        pval = float(mannwhitneyu(set_vals, bg_vals, alternative="two-sided").pvalue)
    except Exception:
        pval = 1.0
    all_vals = np.array(list(ranks.values()), dtype=float)
    sd = float(np.nanstd(all_vals)) or 1.0
    nes = float((np.nanmean(set_vals) - np.nanmean(all_vals)) / sd * math.sqrt(len(set_vals)))
    if nes >= 0:
        core = [g for g, _ in sorted(((g, ranks[g]) for g in gs), key=lambda x: x[1], reverse=True)[: min(200, len(gs))]]
    else:
        core = [g for g, _ in sorted(((g, ranks[g]) for g in gs), key=lambda x: x[1])[: min(200, len(gs))]]
    return {
        "pair_id": pair_row["pair_id"],
        "trait1": pair_row["trait1"],
        "trait2": pair_row["trait2"],
        "database": db_name,
        "pathway_id": pid,
        "pathway_description": desc,
        "setSize": len(gs),
        "NES": nes,
        "pvalue": pval,
        "core_enrichment": "/".join(core[:200]),
    }


def run_gsea_for_db(pair_row: dict, gene_df: pd.DataFrame, db_name: str, gmt_name: str, n_perm: int = 100):
    gmt = GMT_DIR / gmt_name
    fields = ["pair_id", "trait1", "trait2", "database", "pathway_id", "pathway_description", "setSize", "NES", "pvalue", "p.adjust", "qvalue", "core_enrichment"]
    if not gmt.exists():
        return pd.DataFrame(columns=fields)
    vals = pd.to_numeric(gene_df["z"], errors="coerce")
    if vals.notna().sum() == 0:
        vals = -np.log10(pd.to_numeric(gene_df["p"], errors="coerce").clip(lower=1e-300))
    ranks = dict(zip(gene_df["gene_symbol"].astype(str), vals.fillna(0.0)))
    universe = list(ranks.keys())
    terms = list(parse_gmt(gmt))
    n_jobs = max(1, min(int(os.environ.get("SLURM_CPUS_PER_TASK", "1")), len(terms) or 1))
    rows = Parallel(n_jobs=n_jobs, backend="loky", batch_size=8)(
        delayed(gsea_one_set)(pair_row, ranks, universe, pid, desc, genes, db_name, n_perm)
        for pid, desc, genes in terms
    )
    rows = [r for r in rows if r is not None]
    if not rows:
        return pd.DataFrame(columns=fields)
    df = pd.DataFrame(rows)
    q = multipletests(pd.to_numeric(df["pvalue"], errors="coerce").fillna(1.0), method="fdr_bh")[1]
    df["p.adjust"] = q
    df["qvalue"] = q
    return df[fields]


def main():
    row = read_manifest_row()
    pair_id = row["pair_id"]
    start = time.time()
    print(f"STEP8_PAIR\t{pair_id}")

    mtag_in_dir = STEP8_MTAG / "input" / pair_id
    in1 = mtag_in_dir / f"{row['trait1']}.mtag_input.tsv.gz"
    in2 = mtag_in_dir / f"{row['trait2']}.mtag_input.tsv.gz"
    if in1.exists() and in1.stat().st_size > 0:
        n1 = -1
    else:
        n1 = standard_to_mtag_input(Path(row["trait1_standard_gwas"]), in1)
    if in2.exists() and in2.stat().st_size > 0:
        n2 = -1
    else:
        n2 = standard_to_mtag_input(Path(row["trait2_standard_gwas"]), in2)
    print(f"MTAG_INPUT_N\t{n1}\t{n2}")

    full_dir = STEP8_MTAG / "full" / pair_id
    full_dir.mkdir(parents=True, exist_ok=True)
    out_prefix = full_dir / f"{pair_id}_full"
    full_gz = full_dir / "full_mtag.tsv.gz"
    if full_gz.exists() and full_gz.stat().st_size > 0 and (full_dir / "snp_qc.tsv").exists():
        print(f"SKIP_FULL_MTAG\t{full_gz}")
    else:
        cmd = [
            str(PY2),
            str(MTAG_PY),
            "--sumstats",
            f"{in1},{in2}",
            "--out",
            str(out_prefix),
            "--snp_name",
            "SNP",
            "--chr_name",
            "CHR",
            "--bpos_name",
            "BP",
            "--a1_name",
            "A1",
            "--a2_name",
            "A2",
            "--eaf_name",
            "FRQ",
            "--z_name",
            "Z",
            "--p_name",
            "P",
        "--n_name",
        "N",
        "--force",
        "--stream_stdout",
    ]
        run(cmd, full_dir / "mtag.log", cwd=ROOT / "tools/mtag")
        trait1_out = Path(str(out_prefix) + "_trait_1.txt")
        if not trait1_out.exists():
            raise FileNotFoundError(trait1_out)
        gzip_copy(trait1_out, full_gz)
        for maxfdr in full_dir.glob("*maxFDR*"):
            if maxfdr.name != "maxFDR.txt":
                shutil.copyfile(maxfdr, full_dir / "maxFDR.txt")
                break
        summarize_full_mtag(full_gz, full_dir / "snp_qc.tsv")

    gene_result_file = STEP8_MAGMA / "per_pair" / f"{pair_id}.gene_results.tsv"
    if gene_result_file.exists() and gene_result_file.stat().st_size > 0:
        print(f"SKIP_MAGMA\t{gene_result_file}")
        gene_df = pd.read_csv(gene_result_file, sep="\t", dtype=str)
    else:
        pfile, locfile = make_magma_inputs(full_gz, pair_id)
        genes_out = run_magma(pair_id, pfile, locfile)
        gene_df = normalize_gene_results(row, genes_out)
    run_pathways(row, gene_df)

    BENCH.mkdir(parents=True, exist_ok=True)
    with open(BENCH / f"{pair_id}.step8_runtime.tsv", "w") as fh:
        fh.write("pair_id\twalltime_minutes\tn_snps_input_trait1\tn_snps_input_trait2\tn_genes\n")
        fh.write(f"{pair_id}\t{(time.time() - start) / 60:.3f}\t{n1}\t{n2}\t{len(gene_df)}\n")
    print(f"STEP8_DONE\t{pair_id}\twalltime_minutes={(time.time() - start) / 60:.3f}")


if __name__ == "__main__":
    main()
