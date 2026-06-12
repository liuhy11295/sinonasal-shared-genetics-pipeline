#!/usr/bin/env bash
set -euo pipefail
ROOT="${PROJECT_ROOT:?Set PROJECT_ROOT}"
STEP8="$ROOT/results/phase0_extension/step8_magma_gene_pathway"
RES="$STEP8/resources"
mkdir -p "$RES" "$STEP8/logs"
cd "$RES"
log="$STEP8/logs/prepare_resources.log"
: > "$log"
run(){ echo "[RUN] $*" | tee -a "$log"; "$@" 2>&1 | tee -a "$log"; }
# MAGMA binary and hg19 gene locations; all under p010 project only.
if [ ! -x "$RES/magma/magma" ]; then
  mkdir -p "$RES/magma"
  for url in \
    https://ctg.cncr.nl/software/MAGMA/prog/magma_v1.10_static.zip \
    https://ctg.cncr.nl/software/MAGMA/prog/magma_v1.10.zip \
    https://ctg.cncr.nl/software/MAGMA/prog/magma_v1.09_static.zip; do
    echo "Trying $url" | tee -a "$log"
    if wget -q -O "$RES/magma/magma.zip" "$url"; then
      unzip -o "$RES/magma/magma.zip" -d "$RES/magma" >> "$log" 2>&1 || true
      found=$(find "$RES/magma" -type f -name magma -perm -u+x | head -1 || true)
      if [ -z "$found" ]; then found=$(find "$RES/magma" -type f -name magma | head -1 || true); fi
      if [ -n "$found" ]; then chmod +x "$found"; cp "$found" "$RES/magma/magma"; break; fi
    fi
  done
fi
if [ ! -f "$RES/NCBI37.3.gene.loc" ]; then
  if wget -q -O "$RES/NCBI37.3.zip" https://ctg.cncr.nl/software/MAGMA/aux_files/NCBI37.3.zip; then
    unzip -o "$RES/NCBI37.3.zip" -d "$RES" >> "$log" 2>&1 || true
  fi
fi
# Public GMT gene sets from Enrichr mirrors.
mkdir -p "$RES/gmt"
for spec in \
  "GO_Biological_Process_2023.gmt|https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=GO_Biological_Process_2023" \
  "KEGG_2021_Human.gmt|https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=KEGG_2021_Human" \
  "Reactome_2022.gmt|https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=Reactome_2022"; do
  name=${spec%%|*}; url=${spec#*|}
  if [ ! -s "$RES/gmt/$name" ]; then
    echo "Downloading $name" | tee -a "$log"
    wget -q -O "$RES/gmt/$name" "$url" || true
  fi
done
{
  echo -e "resource\tpath\texists\tnote"
  echo -e "magma\t$RES/magma/magma\t$([ -x "$RES/magma/magma" ] && echo true || echo false)\tMAGMA binary"
  echo -e "gene_loc\t$RES/NCBI37.3.gene.loc\t$([ -s "$RES/NCBI37.3.gene.loc" ] && echo true || echo false)\thg19/GRCh37 gene locations"
  for f in "$RES"/gmt/*.gmt; do [ -e "$f" ] && echo -e "gmt\t$f\t$([ -s "$f" ] && echo true || echo false)\tgene set"; done
} > "$STEP8/resource_manifest.tsv"
cat "$STEP8/resource_manifest.tsv"
