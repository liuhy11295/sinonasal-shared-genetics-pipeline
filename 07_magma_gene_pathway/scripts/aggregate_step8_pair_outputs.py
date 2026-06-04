#!/usr/bin/env python3
from pathlib import Path
import pandas as pd
import zipfile

root = Path('/platform_data/p_user/p010/phase0')
step8 = root / 'results/phase0_extension/step8_magma_gene_pathway'
pair_file = step8 / 'pair_ids_step8.tsv'
per_pair = step8 / 'per_pair'
fuma_dir = root / 'results/fuma_upload'
fuma_dir.mkdir(parents=True, exist_ok=True)

pairs_df = pd.read_csv(pair_file, sep='\t')
pairs = pairs_df['pair_id'].dropna().astype(str).tolist() if 'pair_id' in pairs_df.columns else [x.strip() for x in pair_file.read_text().splitlines() if x.strip() and x.strip() != 'pair_id']

def read_tsv(path, cols=None):
    if not path.exists() or path.stat().st_size == 0:
        return pd.DataFrame(columns=cols or [])
    try:
        return pd.read_csv(path, sep='\t')
    except Exception:
        return pd.DataFrame(columns=cols or [])

def add_pair(df, pid):
    if df.empty:
        return df
    if 'pair_id' not in df.columns:
        df.insert(0, 'pair_id', pid)
    else:
        df['pair_id'] = df['pair_id'].fillna(pid)
    return df

outputs = {
    'magma_gene_results.tsv': [],
    'variant_gene_mapping.tsv': [],
    'gene_priority.tsv': [],
    'gene_set_enrichment.tsv': [],
    'key_shared_genes.tsv': [],
    'magma_qc_report.tsv': [],
    'fuma_input_snps_by_pair.tsv': [],
}
missing = []
for pid in pairs:
    pdir = per_pair
    files = {
        'magma_gene_results.tsv': pdir / f'{pid}.magma_gene_results.tsv',
        'variant_gene_mapping.tsv': pdir / f'{pid}.variant_gene_mapping.tsv',
        'gene_priority.tsv': pdir / f'{pid}.gene_priority.tsv',
        'gene_set_enrichment.tsv': pdir / f'{pid}.gene_set_enrichment.tsv',
        'key_shared_genes.tsv': pdir / f'{pid}.key_shared_genes.tsv',
        'magma_qc_report.tsv': pdir / f'{pid}.magma_qc_report.tsv',
        'fuma_input_snps_by_pair.tsv': pdir / f'{pid}.fuma_input_snps.tsv',
    }
    for out, path in files.items():
        df = add_pair(read_tsv(path), pid)
        if df.empty:
            missing.append({'pair_id': pid, 'file': str(path), 'problem': 'missing_or_empty'})
        outputs[out].append(df)

required_cols = {
    'magma_gene_results.tsv': ['pair_id','gene','CHR','START','END','N_SNP','ZSTAT','P','FDR'],
    'variant_gene_mapping.tsv': ['pair_id','SNP','CHR','BP','nearest_gene','distance_to_gene','source','SuSiE_support','PLACO_support','CPASSOC_support'],
    'gene_priority.tsv': ['pair_id','gene','evidence_source','min_p','n_supporting_methods','priority_tier'],
    'gene_set_enrichment.tsv': ['pair_id','gene_set','source_database','n_genes','overlap_genes','P','FDR'],
    'key_shared_genes.tsv': ['pair_id','gene','evidence_source','supporting_pair','supporting_SNP','priority_tier'],
    'magma_qc_report.tsv': ['pair_id','step','status','reason'],
    'fuma_input_snps_by_pair.tsv': ['pair_id','SNP','CHR','BP','P','source'],
}

for out, frames in outputs.items():
    df = pd.concat(frames, ignore_index=True, sort=False) if frames else pd.DataFrame()
    for c in required_cols[out]:
        if c not in df.columns:
            df[c] = ''
    # keep required first, then extras
    ordered = required_cols[out] + [c for c in df.columns if c not in required_cols[out]]
    df = df[ordered]
    if out == 'gene_priority.tsv' and not df.empty:
        if 'evidence_source' in df.columns and 'evidence_type' in df.columns:
            df['evidence_source'] = df['evidence_source'].where(df['evidence_source'].astype(str).str.len() > 0, df['evidence_type'])
        if 'priority_tier' in df.columns and 'evidence_type' in df.columns:
            def _tier(ev):
                parts = set(str(ev).split(';'))
                if 'SuSiE' in parts and ('PLACO' in parts or 'CPASSOC' in parts):
                    return 'Tier1'
                if 'SuSiE' in parts:
                    return 'Tier2'
                if parts and parts != {''}:
                    return 'Tier3'
                return ''
            df['priority_tier'] = df['priority_tier'].where(df['priority_tier'].astype(str).str.len() > 0, df['evidence_type'].map(_tier))
        tier_order = {'Tier1': 1, 'Tier2': 2, 'Tier3': 3}
        df['_tier_order'] = df['priority_tier'].map(tier_order).fillna(99)
        sort_cols = ['_tier_order'] + (['min_p'] if 'min_p' in df.columns else [])
        df = df.sort_values(sort_cols, kind='mergesort').drop(columns=['_tier_order'])
    if out == 'key_shared_genes.tsv' and not df.empty:
        if 'priority_tier' in df.columns and 'priority_level' in df.columns:
            df['priority_tier'] = df['priority_tier'].where(df['priority_tier'].astype(str).str.len() > 0, df['priority_level'])
        if 'evidence_source' in df.columns and 'support_evidence' in df.columns:
            df['evidence_source'] = df['evidence_source'].where(df['evidence_source'].astype(str).str.len() > 0, df['support_evidence'])
        if 'supporting_pair' in df.columns:
            df['supporting_pair'] = df['supporting_pair'].where(df['supporting_pair'].astype(str).str.len() > 0, df['pair_id'])
    if out == 'gene_set_enrichment.tsv' and not df.empty and 'FDR' in df.columns:
        df['_fdr_num'] = pd.to_numeric(df['FDR'], errors='coerce')
        df = df.sort_values(['pair_id','_fdr_num'], kind='mergesort').drop(columns=['_fdr_num'])
    df.to_csv(step8 / out, sep='\t', index=False)

fuma_pair = read_tsv(step8 / 'fuma_input_snps_by_pair.tsv')
for c in ['SNP','CHR','BP','P','source']:
    if c not in fuma_pair.columns:
        fuma_pair[c] = ''
fuma = fuma_pair[['SNP','CHR','BP','P','source']].copy()
fuma = fuma.dropna(subset=['SNP']).drop_duplicates(['SNP','CHR','BP'], keep='first')
fuma.to_csv(step8 / 'fuma_input_snps.tsv', sep='\t', index=False)
fuma.to_csv(fuma_dir / 'fuma_input_snps.tsv', sep='\t', index=False)

readme = fuma_dir / 'README_FUMA_UPLOAD.txt'
readme.write_text('FUMA input SNP file generated from Step 8 per-pair shared SNP sources. Upload fuma_input_snps.tsv to FUMA SNP2GENE; fields: SNP, CHR, BP, P, source. All source files remain under /platform_data/p_user/p010/.\n')
zip_path = fuma_dir / 'fuma_upload_package.zip'
with zipfile.ZipFile(zip_path, 'w', compression=zipfile.ZIP_DEFLATED) as z:
    z.write(fuma_dir / 'fuma_input_snps.tsv', arcname='fuma_input_snps.tsv')
    z.write(readme, arcname='README_FUMA_UPLOAD.txt')

miss_df = pd.DataFrame(missing)
miss_df.to_csv(step8 / 'step8_missing_or_empty_per_pair_files.tsv', sep='\t', index=False)

magma = read_tsv(step8 / 'magma_gene_results.tsv')
pathway = read_tsv(step8 / 'gene_set_enrichment.tsv')
priority = read_tsv(step8 / 'gene_priority.tsv')
fuma_final = read_tsv(step8 / 'fuma_input_snps.tsv')
qc = read_tsv(step8 / 'magma_qc_report.tsv')
summary = pd.DataFrame([
    {'metric': 'pair_count', 'value': len(pairs)},
    {'metric': 'magma_gene_rows', 'value': len(magma)},
    {'metric': 'magma_significant_gene_rows_FDR_lt_0_05', 'value': int((pd.to_numeric(magma.get('FDR', pd.Series(dtype=float)), errors='coerce') < 0.05).sum())},
    {'metric': 'pathway_rows', 'value': len(pathway)},
    {'metric': 'pathway_significant_rows_FDR_lt_0_05', 'value': int((pd.to_numeric(pathway.get('FDR', pd.Series(dtype=float)), errors='coerce') < 0.05).sum())},
    {'metric': 'tier1_gene_rows', 'value': int((priority.get('priority_tier', pd.Series(dtype=str)) == 'Tier1').sum())},
    {'metric': 'tier2_gene_rows', 'value': int((priority.get('priority_tier', pd.Series(dtype=str)) == 'Tier2').sum())},
    {'metric': 'tier3_gene_rows', 'value': int((priority.get('priority_tier', pd.Series(dtype=str)) == 'Tier3').sum())},
    {'metric': 'fuma_input_snp_rows', 'value': len(fuma_final)},
    {'metric': 'qc_fail_rows', 'value': int((qc.get('status', pd.Series(dtype=str)) == 'FAIL').sum())},
    {'metric': 'qc_warn_rows', 'value': int((qc.get('status', pd.Series(dtype=str)) == 'WARN').sum())},
    {'metric': 'missing_or_empty_per_pair_files', 'value': len(miss_df)},
])
summary.to_csv(step8 / 'step8_summary.tsv', sep='\t', index=False)

print(summary.to_csv(sep='\t', index=False), end='')
