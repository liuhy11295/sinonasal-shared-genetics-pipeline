#!/usr/bin/env python3
import csv
import gzip
import urllib.request
from collections import defaultdict
from pathlib import Path


BASE = Path("/home/lhy/nasal")
ROOT = BASE / "results_end/twas_fusion_gtexv8_magma_candidates"
RES = BASE / "resources/fusion_gtexv8"
WEIGHTS = RES / "weights/GTEx_v8"
HGNC = RES / "hgnc_complete_set.txt"
OUT_DIR = ROOT / "input/restricted_weights"
PAIR_OUT_DIR = ROOT / "input/restricted_weights_by_pair"
AUDIT = ROOT / "audit"

HGNC_URL = "https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt"


def open_text(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else path.open()


def read_tsv(path):
    with open_text(path) as fh:
        yield from csv.DictReader(fh, delimiter="\t")


def write_tsv(path, rows, fields):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, delimiter="\t", fieldnames=fields, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


def ensure_hgnc():
    if HGNC.exists() and HGNC.stat().st_size > 0:
        return
    HGNC.parent.mkdir(parents=True, exist_ok=True)
    urllib.request.urlretrieve(HGNC_URL, HGNC)


def load_candidate_symbols():
    path = ROOT / "input/magma_bonf_candidate_pair_genes.tsv"
    rows = list(read_tsv(path))
    return rows, sorted({r["gene_symbol"] for r in rows if r.get("gene_symbol")})


def split_symbols(x):
    if not x:
        return []
    return [v.strip() for v in x.replace("|", ",").split(",") if v.strip()]


def load_hgnc_mapping(candidate_symbols):
    ensure_hgnc()
    wanted = set(candidate_symbols)
    by_symbol = {}
    alias_hits = defaultdict(set)
    for r in read_tsv(HGNC):
        approved = r.get("symbol", "")
        ensembl = (r.get("ensembl_gene_id", "") or "").split(".")[0]
        entrez = r.get("entrez_id", "")
        if approved in wanted and ensembl:
            by_symbol[approved] = {"gene_symbol": approved, "ensembl_id": ensembl, "entrez_id": entrez, "mapping_type": "approved_symbol"}
        for alias in split_symbols(r.get("alias_symbol", "")) + split_symbols(r.get("prev_symbol", "")):
            if alias in wanted and ensembl:
                alias_hits[alias].add((approved, ensembl, entrez))
    for symbol in wanted - set(by_symbol):
        hits = sorted(alias_hits.get(symbol, []))
        if len(hits) == 1:
            approved, ensembl, entrez = hits[0]
            by_symbol[symbol] = {"gene_symbol": symbol, "ensembl_id": ensembl, "entrez_id": entrez, "mapping_type": f"alias_or_previous_symbol:{approved}"}
    return by_symbol, alias_hits


def find_pos_files():
    pos_files = []
    for path in sorted(WEIGHTS.glob("*/*.pos")) + sorted(WEIGHTS.glob("*/*.pos.gz")):
        if ".nofilter." in path.name:
            continue
        tissue = path.parent.name
        if path.name.startswith("GTExv8.EUR."):
            tissue = path.name.replace("GTExv8.EUR.", "").replace(".pos.gz", "").replace(".pos", "")
        pos_files.append((tissue, path))
    return pos_files


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    PAIR_OUT_DIR.mkdir(parents=True, exist_ok=True)
    AUDIT.mkdir(parents=True, exist_ok=True)
    candidate_rows, candidate_symbols = load_candidate_symbols()
    pair_symbols = defaultdict(set)
    for r in candidate_rows:
        if r.get("pair_id") and r.get("gene_symbol"):
            pair_symbols[r["pair_id"]].add(r["gene_symbol"])
    mapping, alias_hits = load_hgnc_mapping(candidate_symbols)
    ensembl_to_symbol = defaultdict(list)
    for symbol, rec in mapping.items():
        ensembl_to_symbol[rec["ensembl_id"]].append(symbol)

    pos_files = find_pos_files()
    tissue_rows = []
    pair_tissue_rows = []
    matched_by_symbol = defaultdict(int)
    duplicated_rows = []
    for tissue, pos in pos_files:
        out = OUT_DIR / f"{tissue}.magma_candidates.pos"
        total = kept = dup = 0
        pair_handles = {}
        pair_writers = {}
        pair_counts = defaultdict(int)
        with open_text(pos) as fh:
            reader0 = csv.DictReader(fh, delimiter="\t")
            fields0 = reader0.fieldnames or []
        for pair_id in pair_symbols:
            pair_dir = PAIR_OUT_DIR / pair_id
            pair_dir.mkdir(parents=True, exist_ok=True)
            fhp = (pair_dir / f"{tissue}.magma_candidates.pos").open("w", newline="")
            pair_handles[pair_id] = fhp
            writer = csv.DictWriter(fhp, delimiter="\t", fieldnames=fields0, extrasaction="ignore")
            writer.writeheader()
            pair_writers[pair_id] = writer
        with open_text(pos) as fh, out.open("w", newline="") as ofh:
            reader = csv.DictReader(fh, delimiter="\t")
            fields = reader.fieldnames or []
            writer = csv.DictWriter(ofh, delimiter="\t", fieldnames=fields, extrasaction="ignore")
            writer.writeheader()
            for row in reader:
                total += 1
                ensembl = (row.get("ID", "") or "").split(".")[0]
                symbols = ensembl_to_symbol.get(ensembl, [])
                if not symbols:
                    continue
                if len(symbols) > 1:
                    dup += 1
                    duplicated_rows.append({"tissue": tissue, "ensembl_id": ensembl, "gene_symbols": ",".join(symbols), "pos_id": row.get("ID", "")})
                writer.writerow(row)
                kept += 1
                for symbol in symbols:
                    matched_by_symbol[symbol] += 1
                for pair_id, symbols_for_pair in pair_symbols.items():
                    if any(symbol in symbols_for_pair for symbol in symbols):
                        pair_writers[pair_id].writerow(row)
                        pair_counts[pair_id] += 1
        for pair_id, fhp in pair_handles.items():
            fhp.close()
            pair_tissue_rows.append({
                "pair_id": pair_id,
                "tissue": tissue,
                "restricted_pos": str(PAIR_OUT_DIR / pair_id / f"{tissue}.magma_candidates.pos"),
                "n_candidate_weight_rows_kept": str(pair_counts[pair_id]),
            })
        tissue_rows.append({
            "tissue": tissue,
            "source_pos": str(pos),
            "restricted_pos": str(out),
            "n_weight_rows_total": str(total),
            "n_candidate_weight_rows_kept": str(kept),
            "n_duplicated_symbol_mappings": str(dup),
        })

    gene_rows = []
    for symbol in candidate_symbols:
        rec = mapping.get(symbol)
        gene_rows.append({
            "gene_symbol": symbol,
            "mapping_status": "mapped" if rec else "unmapped",
            "ensembl_id": rec["ensembl_id"] if rec else "",
            "entrez_id": rec["entrez_id"] if rec else "",
            "mapping_type": rec["mapping_type"] if rec else "",
            "n_available_weight_rows_across_downloaded_tissues": str(matched_by_symbol.get(symbol, 0)),
            "ambiguous_alias_hits": ";".join([":".join(x) for x in sorted(alias_hits.get(symbol, []))]) if symbol not in mapping else "",
        })

    write_tsv(AUDIT / "weight_mapping_audit.tsv", gene_rows, [
        "gene_symbol", "mapping_status", "ensembl_id", "entrez_id", "mapping_type",
        "n_available_weight_rows_across_downloaded_tissues", "ambiguous_alias_hits",
    ])
    write_tsv(AUDIT / "restricted_weight_tissue_audit.tsv", tissue_rows, [
        "tissue", "source_pos", "restricted_pos", "n_weight_rows_total",
        "n_candidate_weight_rows_kept", "n_duplicated_symbol_mappings",
    ])
    write_tsv(AUDIT / "restricted_weight_pair_tissue_audit.tsv", pair_tissue_rows, [
        "pair_id", "tissue", "restricted_pos", "n_candidate_weight_rows_kept",
    ])
    write_tsv(AUDIT / "duplicated_weight_mappings.tsv", duplicated_rows, ["tissue", "ensembl_id", "gene_symbols", "pos_id"])
    write_tsv(AUDIT / "unmapped_genes.tsv", [r for r in gene_rows if r["mapping_status"] == "unmapped"], list(gene_rows[0].keys()))
    print(f"candidate_symbols={len(candidate_symbols)}")
    print(f"mapped_symbols={sum(1 for r in gene_rows if r['mapping_status'] == 'mapped')}")
    print(f"available_tissues={len(pos_files)}")


if __name__ == "__main__":
    main()
