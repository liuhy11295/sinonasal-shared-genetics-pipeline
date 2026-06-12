#!/usr/bin/env bash
set -euo pipefail
ROOT="${PROJECT_ROOT:?Set PROJECT_ROOT}"
STEP8="$ROOT/results/phase0_extension/step8_magma_gene_pathway"
REF="$ROOT/data/reference/lava/1000G_Phase3_plinkfiles/1000G_EUR_Phase3_plink"
OUT="$STEP8/resources/plink_merged/1000G.EUR.QC"
mkdir -p "$(dirname "$OUT")"
if [ -s "$OUT.bed" ] && [ -s "$OUT.bim" ] && [ -s "$OUT.fam" ]; then
  echo "merged plink exists: $OUT"
  exit 0
fi
LIST="$STEP8/resources/plink_merged/merge_list.txt"
: > "$LIST"
for chr in $(seq 2 22); do echo "$REF/1000G.EUR.QC.$chr" >> "$LIST"; done
plink --bfile "$REF/1000G.EUR.QC.1" --merge-list "$LIST" --make-bed --out "$OUT"
