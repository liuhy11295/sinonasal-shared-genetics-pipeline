#!/bin/bash
set -euo pipefail
BASE=/platform_data/p_user/p010
PROJECT=$BASE/phase0
PKG=$PROJECT/results/phase0_server/nasal4_vs_other_package
LDSC_PKG=$PROJECT/results/phase0_server/ldsc_package
OUT=$PROJECT/results/phase0_extension/step0_1
INPUTS=$OUT/inputs
MAN=$OUT/manifests/data_manifest.tsv
mkdir -p "$INPUTS" "$INPUTS/standardized_gwas" "$INPUTS/ldsc_sumstats" "$INPUTS/ldsc_results" "$INPUTS/lava_results" "$INPUTS/mtag_results" "$INPUTS/coloc_results" "$INPUTS/reference" "$INPUTS/reference/ld_plink" "$INPUTS/reference/ldscore" "$OUT/manifests" "$OUT/checks"
: > "$MAN"
printf 'data_type\ttrait\tpair_id\toriginal_path\tsymlink_path\tfile_exists\tis_symlink\tkey_fields\tnote\n' >> "$MAN"

safe_ln() {
  local src="$1" dst="$2"
  if [ -e "$dst" ] || [ -L "$dst" ]; then
    if [ "$(readlink "$dst" 2>/dev/null || true)" != "$src" ]; then
      echo "ERROR: existing path is not expected symlink: $dst" >&2
      return 2
    fi
  else
    ln -s "$src" "$dst"
  fi
}

row() {
  local data_type="$1" trait="$2" pair_id="$3" src="$4" dst="$5" fields="$6" note="$7"
  local exists=FALSE issym=FALSE
  [ -e "$dst" ] && exists=TRUE
  [ -L "$dst" ] && issym=TRUE
  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$data_type" "$trait" "$pair_id" "$src" "$dst" "$exists" "$issym" "$fields" "$note" >> "$MAN"
}

# standardized GWAS summary statistics
find "$PROJECT/data/phase0_v3_munged" -maxdepth 1 -type f -name '*.common.tsv.gz' | sort | while read -r src; do
  b=$(basename "$src")
  trait=${b%.common.tsv.gz}
  dst="$INPUTS/standardized_gwas/$b"
  safe_ln "$src" "$dst"
  row standardized_gwas "$trait" NA "$src" "$dst" 'SNP,CHR,BP,A1,A2,BETA,SE,P,EAF,N' 'A1/A2 correspond to effect/other allele fields used by current pipeline.'
done

# LDSC munged summary statistics
find "$LDSC_PKG/sumstats" -maxdepth 1 -type f -name '*.sumstats.gz' | sort | while read -r src; do
  b=$(basename "$src")
  trait=${b%.sumstats.gz}
  dst="$INPUTS/ldsc_sumstats/$b"
  safe_ln "$src" "$dst"
  row ldsc_sumstats "$trait" NA "$src" "$dst" 'SNP,A1,A2,Z,N' 'Existing LDSC munged sumstats.'
done

# Existing LDSC results
for src in \
  "$PKG/nasal4_vs_other_ldsc_rg_summary.tsv" \
  "$PKG/nasal4_vs_other_ldsc_h2_summary.tsv" \
  "$PKG/nasal4_vs_other_ldsc_verify.tsv"; do
  [ -e "$src" ] || continue
  b=$(basename "$src")
  dst="$INPUTS/ldsc_results/$b"
  safe_ln "$src" "$dst"
  row ldsc_result NA NA "$src" "$dst" 'pair_id,trait_a,trait_b,rg,p,h2_obs' 'Existing LDSC result or verification table.'
done

# Existing LAVA results
for src in \
  "$PKG/figure_local_rg_edges.tsv" \
  "$PKG/phase0_v3_local_rg_locus_qc.tsv" \
  "$PKG/nasal4_vs_other_lava_verify.tsv"; do
  [ -e "$src" ] || continue
  b=$(basename "$src")
  dst="$INPUTS/lava_results/$b"
  safe_ln "$src" "$dst"
  row lava_result NA NA "$src" "$dst" 'region_id,chr,start,end,trait_pair,local_rg,p,p_adj' 'Existing LAVA result; trait_pair can be split into trait1/trait2; p_adj is q-like adjusted P.'
done

# Existing MTAG results: link manifest/figure table plus per-pair output directory as directory symlink
for src in \
  "$PKG/figure_mtag_manhattan.tsv" \
  "$PKG/nasal4_vs_other_mtag_verify.tsv"; do
  [ -e "$src" ] || continue
  b=$(basename "$src")
  dst="$INPUTS/mtag_results/$b"
  safe_ln "$src" "$dst"
  row mtag_result NA NA "$src" "$dst" 'chr,pos,rsid,trait_pair,p_mtag,beta_a,beta_b' 'Existing MTAG figure/filter table.'
done
if [ -d "$PKG/mtag" ]; then
  dst="$INPUTS/mtag_results/mtag"
  safe_ln "$PKG/mtag" "$dst"
  row mtag_result_dir NA NA "$PKG/mtag" "$dst" 'SNP,CHR,BP,A1,A2,Z,N,FRQ,mtag_beta,mtag_se,mtag_pval' 'Directory symlink to existing per-pair MTAG outputs.'
fi

# Existing coloc results
for src in \
  "$PKG/phase0_v3_coloc_locus_manifest.tsv" \
  "$PKG/figure_candidate_variant_gene_map.tsv" \
  "$PKG/nasal4_vs_other_coloc_verify.tsv"; do
  [ -e "$src" ] || continue
  b=$(basename "$src")
  dst="$INPUTS/coloc_results/$b"
  safe_ln "$src" "$dst"
  row coloc_result NA NA "$src" "$dst" 'pair_id,locus_id,chr,start,end,PP.H3,PP.H4,SNP' 'Existing coloc manifest/candidate/verify table; PP columns are in coloc_merged outputs.'
done
if [ -d "$PKG/coloc_merged" ]; then
  dst="$INPUTS/coloc_results/coloc_merged"
  safe_ln "$PKG/coloc_merged" "$dst"
  row coloc_result_dir NA NA "$PKG/coloc_merged" "$dst" 'nsnps,PP.H3.abf,PP.H4.abf,locus_id,pair_id' 'Directory symlink to existing coloc outputs.'
fi

# References
LDPLINK="$PROJECT/data/reference/lava/1000G_Phase3_plinkfiles/1000G_EUR_Phase3_plink"
if [ -d "$LDPLINK" ]; then
  dst="$INPUTS/reference/ld_plink/1000G_EUR_Phase3_plink"
  safe_ln "$LDPLINK" "$dst"
  row ld_reference_plink NA NA "$LDPLINK" "$dst" '.bed,.bim,.fam' 'Directory symlink to existing 1000G EUR PLINK reference.'
fi

for src in \
  "$PROJECT/data/reference/ldsc/1000G_Phase3_ldscores/LDscore" \
  "$PROJECT/data/reference/ldsc/1000G_Phase3_weights_hm3_no_MHC/1000G_Phase3_weights_hm3_no_MHC" \
  "$PROJECT/data/reference/ldsc/w_hm3.snplist"; do
  [ -e "$src" ] || continue
  b=$(basename "$src")
  dst="$INPUTS/reference/ldscore/$b"
  safe_ln "$src" "$dst"
  row ldscore_reference NA NA "$src" "$dst" '.l2.ldscore.gz,.M,.M_5_50' 'Existing LDSC LD score/weight reference.'
done

# Validate symlinks only inside current input tree.
find "$INPUTS" -type l -print | sort > "$OUT/checks/symlink_paths.txt"
{
  printf 'symlink_path\ttarget\ttarget_exists\n'
  while read -r lnk; do
    tgt=$(readlink "$lnk")
    if [ -e "$lnk" ]; then ok=TRUE; else ok=FALSE; fi
    printf '%s\t%s\t%s\n' "$lnk" "$tgt" "$ok"
  done < "$OUT/checks/symlink_paths.txt"
} > "$OUT/checks/symlink_validation.tsv"

awk -F'\t' 'NR>1 && $3!="TRUE" {bad++} END {if (bad>0) {print "BROKEN_SYMLINKS="bad; exit 1} else print "BROKEN_SYMLINKS=0"}' "$OUT/checks/symlink_validation.tsv" > "$OUT/checks/symlink_validation.summary.txt"

echo "STEP1_DONE=$(date -Is)" > "$OUT/checks/step1.done"
