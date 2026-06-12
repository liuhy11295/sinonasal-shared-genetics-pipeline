#!/usr/bin/env python3
import csv
import gzip
import math
import argparse
import os
import re
import shutil
import subprocess
from collections import Counter, defaultdict
from pathlib import Path


project_root = os.environ.get("NASAL_PROJECT_ROOT", "").strip()
if not project_root:
    raise SystemExit("Set NASAL_PROJECT_ROOT to the parent directory containing results_end and resources.")
BASE = Path(project_root)
RESULTS = BASE / "results_end"
OUT = RESULTS / "twas_fusion_gtexv8_magma_candidates"
RES = BASE / "resources/fusion_gtexv8"
FUSION = RES / "fusion_twas"
WEIGHTS = RES / "weights/GTEx_v8"
LDREF = RES / "LDREF"
TWAS = RESULTS / "twas"
MTAG_FULL = RESULTS / "mtag/full_data"
OFFICIAL_PAGE = OUT / "audit/fusion_official_page_http.html"


GTEX49 = [
    "Adipose_Subcutaneous",
    "Adipose_Visceral_Omentum",
    "Adrenal_Gland",
    "Artery_Aorta",
    "Artery_Coronary",
    "Artery_Tibial",
    "Brain_Amygdala",
    "Brain_Anterior_cingulate_cortex_BA24",
    "Brain_Caudate_basal_ganglia",
    "Brain_Cerebellar_Hemisphere",
    "Brain_Cerebellum",
    "Brain_Cortex",
    "Brain_Frontal_Cortex_BA9",
    "Brain_Hippocampus",
    "Brain_Hypothalamus",
    "Brain_Nucleus_accumbens_basal_ganglia",
    "Brain_Putamen_basal_ganglia",
    "Brain_Spinal_cord_cervical_c-1",
    "Brain_Substantia_nigra",
    "Breast_Mammary_Tissue",
    "Cells_Cultured_fibroblasts",
    "Cells_EBV-transformed_lymphocytes",
    "Colon_Sigmoid",
    "Colon_Transverse",
    "Esophagus_Gastroesophageal_Junction",
    "Esophagus_Mucosa",
    "Esophagus_Muscularis",
    "Heart_Atrial_Appendage",
    "Heart_Left_Ventricle",
    "Kidney_Cortex",
    "Liver",
    "Lung",
    "Minor_Salivary_Gland",
    "Muscle_Skeletal",
    "Nerve_Tibial",
    "Ovary",
    "Pancreas",
    "Pituitary",
    "Prostate",
    "Skin_Not_Sun_Exposed_Suprapubic",
    "Skin_Sun_Exposed_Lower_leg",
    "Small_Intestine_Terminal_Ileum",
    "Spleen",
    "Stomach",
    "Testis",
    "Thyroid",
    "Uterus",
    "Vagina",
    "Whole_Blood",
]


def write_tsv(path, rows, fields):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, delimiter="\t", fieldnames=fields, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


def fnum(x):
    try:
        if x is None or str(x).strip() == "":
            return None
        return float(x)
    except Exception:
        return None


def rsid(x):
    return isinstance(x, str) and x.startswith("rs")


def open_text(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else path.open()


def run(cmd, cwd=None):
    return subprocess.run(cmd, cwd=cwd, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)


def check_fusion():
    rows = []
    if not (FUSION / ".git").exists() and not (FUSION / "FUSION.assoc_test.R").exists():
        if any(FUSION.iterdir()):
            pass
        else:
            res = run(["git", "clone", "--depth", "1", "https://github.com/gusevlab/fusion_twas.git", str(FUSION)])
            rows.append({"resource": "fusion_clone", "status": "ok" if res.returncode == 0 else "failed", "detail": res.stderr.strip()})
    for rel in ["FUSION.assoc_test.R", "FUSION.post_process.R", "utils/plink_utils.R"]:
        rows.append({
            "resource": rel,
            "status": "present" if (FUSION / rel).exists() else "missing",
            "detail": str(FUSION / rel),
        })
    write_tsv(OUT / "audit/fusion_resource_check.tsv", rows, ["resource", "status", "detail"])
    return rows


def check_weights():
    rows = []
    commands = []
    pos_index = {}
    for p in list(WEIGHTS.glob("*/*.pos")) + list(WEIGHTS.glob("*/*.pos.gz")):
        if ".nofilter." in p.name:
            continue
        m = re.match(r"GTExv8\.EUR\.(.+)\.pos(?:\.gz)?$", p.name)
        if m:
            pos_index[m.group(1)] = p
        else:
            pos_index[p.parent.name] = p
    for tissue in GTEX49:
        candidates = []
        for root in [WEIGHTS, WEIGHTS / tissue]:
            candidates.extend([
                root / f"GTExv8.EUR.{tissue}.pos",
                root / f"GTExv8.EUR.{tissue}.pos.gz",
                root / f"{tissue}.pos",
                root / f"{tissue}.pos.gz",
                root / f"GTEx.{tissue}.pos",
                root / f"GTEx.{tissue}.pos.gz",
            ])
        pos = pos_index.get(tissue) or next((p for p in candidates if p.exists()), None)
        wdir = None
        if pos:
            wdir = pos.parent
        n_wgt = len(list(wdir.glob("*.wgt.RDat"))) if wdir else 0
        if wdir and n_wgt == 0:
            nested = wdir / f"GTExv8.EUR.{tissue}"
            n_wgt = len(list(nested.glob("*.wgt.RDat"))) if nested.exists() else 0
        rows.append({
            "tissue": tissue,
            "pos_status": "present" if pos else "missing",
            "pos_path": str(pos or ""),
            "weights_dir": str(wdir or ""),
            "n_wgt_RDat": str(n_wgt),
            "weights_status": "present" if n_wgt > 0 else "missing",
        })
        if not pos or n_wgt == 0:
            url = f"https://s3.us-west-1.amazonaws.com/gtex.v8.fusion/EUR/GTExv8.EUR.{tissue}.tar.gz"
            archive = WEIGHTS / "archives" / f"GTExv8.EUR.{tissue}.tar.gz"
            target = WEIGHTS / tissue
            commands.append({
                "tissue": tissue,
                "url": url,
                "archive": str(archive),
                "target_dir": str(target),
                "command": f"mkdir -p {archive.parent} {target} && wget -c -O {archive} {url} && tar -xzf {archive} -C {target}",
            })
    write_tsv(OUT / "audit/weight_resource_check.tsv", rows, ["tissue", "pos_status", "pos_path", "weights_dir", "n_wgt_RDat", "weights_status"])
    write_tsv(OUT / "download_manifest.tsv", commands, ["tissue", "url", "archive", "target_dir", "command"])
    with (OUT / "download_commands.sh").open("w") as fh:
        fh.write("#!/usr/bin/env bash\nset -euo pipefail\n")
        fh.write(f"mkdir -p {WEIGHTS / 'archives'}\n")
        for r in commands:
            fh.write(r["command"] + "\n")
    return rows


def check_ldref():
    rows = []
    commands = []
    for chrom in range(1, 23):
        prefix = LDREF / f"1000G.EUR.{chrom}"
        bed = Path(str(prefix) + ".bed")
        bim = Path(str(prefix) + ".bim")
        fam = Path(str(prefix) + ".fam")
        ok = all(p.exists() for p in [bed, bim, fam])
        snp_id_sample = ""
        if bim.exists():
            with bim.open() as fh:
                for line in fh:
                    parts = line.split()
                    if len(parts) > 1:
                        snp_id_sample = parts[1]
                        break
        rows.append({
            "chr": str(chrom),
            "bed": str(bed.exists()),
            "bim": str(bim.exists()),
            "fam": str(fam.exists()),
            "complete": str(ok),
            "snp_id_sample": snp_id_sample,
            "rsid_like": str(snp_id_sample.startswith("rs")),
        })
    if not all(r["complete"] == "True" for r in rows):
        url = "https://data.broadinstitute.org/alkesgroup/FUSION/LDREF.tar.bz2"
        commands.append({"resource": "FUSION_LDREF", "url": url, "command": f"wget -c -O {RES / 'LDREF.tar.bz2'} {url} && tar -xjf {RES / 'LDREF.tar.bz2'} -C {RES}"})
    write_tsv(OUT / "audit/ldref_resource_check.tsv", rows, ["chr", "bed", "bim", "fam", "complete", "snp_id_sample", "rsid_like"])
    write_tsv(OUT / "ldref_download_commands.tsv", commands, ["resource", "url", "command"])
    return rows


def prepare_candidates():
    src = TWAS / "gene_grade_comparison.tsv"
    if not src.exists():
        raise SystemExit(f"Missing {src}")
    rows = []
    for r in csv.DictReader(src.open(), delimiter="\t"):
        if r.get("locus_grade") not in {"A", "B"}:
            continue
        if r.get("locus_ab_magma_bonf_positive") != "True":
            continue
        rows.append({
            "pair_id": r["pair_id"],
            "gene_symbol": r["gene_symbol"],
            "locus_id": r.get("best_locus_id", ""),
            "locus_grade": r.get("locus_grade", ""),
            "magma_p": r.get("magma_p", ""),
            "magma_bonf_threshold": r.get("locus_ab_magma_bonf_threshold", ""),
            "magma_bonf_sig": "True",
            "chr": r.get("chr", ""),
            "gene_start": r.get("start", ""),
            "gene_end": r.get("stop", ""),
            "gene_id": r.get("gene_id", ""),
        })
    fields = ["pair_id", "gene_symbol", "locus_id", "locus_grade", "magma_p", "magma_bonf_threshold", "magma_bonf_sig", "chr", "gene_start", "gene_end", "gene_id"]
    write_tsv(OUT / "input/magma_bonf_candidate_pair_genes.tsv", rows, fields)
    unique = sorted({r["gene_symbol"] for r in rows})
    (OUT / "input/magma_bonf_candidate_unique_genes.txt").write_text("\n".join(unique) + "\n")
    audit = [{
        "n_pair_gene_records": str(len(rows)),
        "n_unique_genes": str(len(unique)),
        "n_pairs": str(len({r["pair_id"] for r in rows})),
        "all_locus_AB": str(all(r["locus_grade"] in {"A", "B"} for r in rows)),
        "all_magma_bonf_positive": str(all(r["magma_bonf_sig"] == "True" for r in rows)),
    }]
    write_tsv(OUT / "audit/candidate_gene_audit.tsv", audit, list(audit[0].keys()))
    return rows


def normalize_sumstats(candidate_rows):
    pairs = sorted({r["pair_id"] for r in candidate_rows})
    audit = []
    missing = []
    for pair in pairs:
        src = MTAG_FULL / pair / "full_mtag.tsv.gz"
        if not src.exists():
            missing.append({"pair_id": pair, "reason": "missing_mtag_full_data", "path": str(src)})
            audit.append({"pair_id": pair, "source": "", "input_snp_count": "0", "usable_snp_count": "0", "missing_z_count": "0", "ambiguous_snp_count": "0", "has_N": "False", "status": "stop_missing_source"})
            continue
        with open_text(src) as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            cols = reader.fieldnames or []
            zcol = next((c for c in ["Z", "z", "Z_mtag", "mtag_z", "Z_META", "z_mtag"] if c in cols), None)
            snpcol = next((c for c in ["SNP", "snp", "rsid", "RSID"] if c in cols), None)
            a1col = next((c for c in ["A1", "EA", "effect_allele", "ALT"] if c in cols), None)
            a2col = next((c for c in ["A2", "NEA", "other_allele", "REF"] if c in cols), None)
            betacol = next((c for c in ["BETA", "beta", "b"] if c in cols), None)
            secol = next((c for c in ["SE", "se"] if c in cols), None)
            ncol = next((c for c in ["N", "n"] if c in cols), None)
            if not all([zcol, snpcol, a1col, a2col]):
                missing.append({"pair_id": pair, "reason": f"missing_required_cols cols={','.join(cols)}", "path": str(src)})
                audit.append({"pair_id": pair, "source": str(src), "input_snp_count": "0", "usable_snp_count": "0", "missing_z_count": "0", "ambiguous_snp_count": "0", "has_N": str(bool(ncol)), "status": "stop_missing_required_columns"})
                continue
            out = OUT / "input/sumstats" / f"{pair}.fusion.sumstats"
            input_n = usable = missz = ambig = 0
            with out.open("w", newline="") as ofh:
                fields = ["SNP", "A1", "A2", "Z"] + ([x for x in ["BETA", "SE", "N"] if (x == "BETA" and betacol) or (x == "SE" and secol) or (x == "N" and ncol)])
                writer = csv.DictWriter(ofh, delimiter="\t", fieldnames=fields)
                writer.writeheader()
                for r in reader:
                    input_n += 1
                    snp, a1, a2, z = r.get(snpcol, ""), r.get(a1col, ""), r.get(a2col, ""), r.get(zcol, "")
                    alleles = {a1.upper(), a2.upper()}
                    if not snp or not rsid(snp) or not a1 or not a2:
                        continue
                    if alleles in [{"A", "T"}, {"C", "G"}]:
                        ambig += 1
                        continue
                    if fnum(z) is None:
                        missz += 1
                        continue
                    rec = {"SNP": snp, "A1": a1.upper(), "A2": a2.upper(), "Z": z}
                    if betacol:
                        rec["BETA"] = r.get(betacol, "")
                    if secol:
                        rec["SE"] = r.get(secol, "")
                    if ncol:
                        rec["N"] = r.get(ncol, "")
                    writer.writerow(rec)
                    usable += 1
            audit.append({"pair_id": pair, "source": str(src), "input_snp_count": str(input_n), "usable_snp_count": str(usable), "missing_z_count": str(missz), "ambiguous_snp_count": str(ambig), "has_N": str(bool(ncol)), "status": "run" if usable else "stop_no_usable_snps"})
    write_tsv(OUT / "audit/sumstats_audit.tsv", audit, ["pair_id", "source", "input_snp_count", "usable_snp_count", "missing_z_count", "ambiguous_snp_count", "has_N", "status"])
    write_tsv(OUT / "missing_z_summary.tsv", missing, ["pair_id", "reason", "path"])
    return audit


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--normalize-sumstats",
        action="store_true",
        help="Generate FUSION sumstats files from MTAG full_data. Off by default to keep the first-stage input MAGMA-candidate restricted.",
    )
    args = parser.parse_args()
    for d in [OUT / "audit", OUT / "input/sumstats", OUT / "input/restricted_weights", OUT / "scripts"]:
        d.mkdir(parents=True, exist_ok=True)
    fusion_rows = check_fusion()
    weight_rows = check_weights()
    ld_rows = check_ldref()
    candidates = prepare_candidates()
    if args.normalize_sumstats:
        sumstats = normalize_sumstats(candidates)
    else:
        sumstats = []
        write_tsv(
            OUT / "audit/sumstats_audit.tsv",
            [{"pair_id": "", "source": "", "input_snp_count": "", "usable_snp_count": "", "missing_z_count": "", "ambiguous_snp_count": "", "has_N": "", "status": "not_generated_run_prepare_with_--normalize-sumstats"}],
            ["pair_id", "source", "input_snp_count", "usable_snp_count", "missing_z_count", "ambiguous_snp_count", "has_N", "status"],
        )
        write_tsv(OUT / "missing_z_summary.tsv", [], ["pair_id", "reason", "path"])

    missing_weights = [r["tissue"] for r in weight_rows if r["pos_status"] != "present" or r["weights_status"] != "present"]
    missing_ld = [r["chr"] for r in ld_rows if r["complete"] != "True"]
    missing_fusion = [r["resource"] for r in fusion_rows if r["status"] in {"missing", "failed"}]
    blockers = []
    if missing_fusion:
        blockers.append(f"FUSION missing/failed: {', '.join(missing_fusion)}")
    if missing_weights:
        blockers.append(f"GTEx v8 weights incomplete: {len(missing_weights)} tissues missing")
    if missing_ld:
        blockers.append(f"LDREF incomplete: chromosomes {', '.join(missing_ld)}")
    if args.normalize_sumstats and any(r["status"].startswith("stop") for r in sumstats):
        blockers.append("Some pair sumstats missing usable SNP/A1/A2/Z")
    build_mismatch = OUT / "audit/build_mismatch_report.tsv"
    if build_mismatch.exists():
        blockers.append("Build mismatch detected; see audit/build_mismatch_report.tsv")
    status = "blocked_missing_resources" if blockers else "ready_for_weight_mapping_and_trial"
    report = [
        "# Missing Resource Report",
        "",
        f"Status: `{status}`",
        "",
        "## Blockers",
        "",
    ]
    report.extend([f"- {b}" for b in blockers] or ["- None"])
    report.extend([
        "",
        "## Candidate Genes",
        "",
        f"- pair-gene records: {len(candidates)}",
        f"- unique genes: {len({r['gene_symbol'] for r in candidates})}",
        f"- sumstats normalized: {args.normalize_sumstats}",
        "",
        "## Download files",
        "",
        "- `download_manifest.tsv`",
        "- `download_commands.sh`",
        "- `ldref_download_commands.tsv`",
    ])
    if build_mismatch.exists():
        report.extend([
            "",
            "## Build Mismatch Stop",
            "",
            "Status: `blocked_build_mismatch`",
            "",
            "GTEx v8 EUR FUSION weight SNP coordinates were observed to be hg38/GRCh38-like while the current FUSION LDREF uses hg19/GRCh37-like coordinates.",
            "Per the analysis constraint, FUSION TWAS jobs must not be launched until matched resources are supplied or a documented liftover/matched-LDREF strategy is approved.",
            "",
            "See `audit/build_mismatch_report.tsv`.",
        ])
    (OUT / "missing_resource_report.md").write_text("\n".join(report) + "\n")
    print(status)
    print("candidate_pair_gene_records", len(candidates))
    print("candidate_unique_genes", len({r["gene_symbol"] for r in candidates}))
    print("missing_weight_tissues", len(missing_weights))
    print("missing_ld_chromosomes", len(missing_ld))


if __name__ == "__main__":
    main()
