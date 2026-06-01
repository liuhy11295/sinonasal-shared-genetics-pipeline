#!/usr/bin/env python3
from pathlib import Path
import pandas as pd
import numpy as np

ROOT = Path('/platform_data/p_user/p010/phase0')
RESULTS = ROOT / 'results'
MAGMA = RESULTS / 'step8_magma' / 'pair_specific_gene_results.tsv'
OUTDIR = RESULTS / 'magma'
OUTDIR.mkdir(parents=True, exist_ok=True)

EXCLUDE = 'CHRONIC_RHINITIS_PANUKB_J31'

def read_tsv(path):
    if not path.exists() or path.stat().st_size == 0:
        return pd.DataFrame()
    return pd.read_csv(path, sep='\t', dtype=str, low_memory=False)

def num(s):
    return pd.to_numeric(s, errors='coerce')

def filter_exclude(df):
    if df.empty:
        return df
    mask = pd.Series(False, index=df.index)
    for c in ['pair_id','trait1','trait2','nasal_trait','partner_trait']:
        if c in df.columns:
            mask |= df[c].fillna('').astype(str).str.contains(EXCLUDE, regex=False)
    return df.loc[~mask].copy()

genes = filter_exclude(read_tsv(MAGMA))
if genes.empty:
    raise SystemExit(f'Missing or empty MAGMA results: {MAGMA}')

for col in ['pair_id','trait1','trait2']:
    if col not in genes.columns:
        genes[col] = ''

p_col = 'p' if 'p' in genes.columns else ('P' if 'P' in genes.columns else None)
if p_col is None:
    raise SystemExit('MAGMA gene results lack p/P column')

genes['_p'] = num(genes[p_col])
if 'gene_bonf_p' in genes.columns:
    genes['_bonf'] = num(genes['gene_bonf_p'])
elif 'n_tested_genes' in genes.columns:
    genes['_bonf'] = genes['_p'] * num(genes['n_tested_genes'])
else:
    genes['_bonf'] = np.nan

if 'gene_id' in genes.columns:
    gene_key = 'gene_id'
elif 'gene_symbol' in genes.columns:
    gene_key = 'gene_symbol'
else:
    gene_key = genes.columns[0]

# Pair-specific BH/FDR on MAGMA gene P-values.
genes['_fdr'] = np.nan
for pair, idx in genes.groupby('pair_id').groups.items():
    vals = genes.loc[idx, '_p']
    valid = vals.dropna().sort_values()
    if valid.empty:
        continue
    m = len(valid)
    q = valid.to_numpy() * m / np.arange(1, m + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    q = np.minimum(q, 1.0)
    genes.loc[valid.index, '_fdr'] = q

summary = genes.groupby(['pair_id','trait1','trait2'], dropna=False).agg(
    tested_genes=(gene_key, 'nunique'),
    bonferroni_sig_genes=('_bonf', lambda x: int((x <= 0.05).sum())),
    fdr05_genes=('_fdr', lambda x: int((x <= 0.05).sum())),
    nominal_p05_genes=('_p', lambda x: int((x < 0.05).sum())),
    min_p=('_p', 'min'),
    min_bonf_p=('_bonf', 'min'),
    min_fdr=('_fdr', 'min'),
).reset_index().sort_values(['bonferroni_sig_genes','fdr05_genes','nominal_p05_genes'], ascending=[False,False,False])

summary.to_csv(OUTDIR / 'pair_gene_count_summary.tsv', sep='\t', index=False)

threshold_rows = []
for label, mask in [
    ('Bonferroni<=0.05', genes['_bonf'] <= 0.05),
    ('FDR<0.05_pairwise_BH', genes['_fdr'] < 0.05),
    ('Nominal_P<0.05', genes['_p'] < 0.05),
]:
    sub = genes.loc[mask.fillna(False)].copy()
    threshold_rows.append({
        'threshold': label,
        'n_gene_rows': len(sub),
        'n_unique_genes': sub[gene_key].nunique() if gene_key in sub.columns else len(sub),
        'n_pairs_with_any_gene': sub['pair_id'].nunique() if 'pair_id' in sub.columns else 0,
        'top_pair_by_count': sub['pair_id'].value_counts().index[0] if not sub.empty else '',
        'top_pair_gene_count': int(sub['pair_id'].value_counts().iloc[0]) if not sub.empty else 0,
    })
threshold = pd.DataFrame(threshold_rows)
threshold.to_csv(OUTDIR / 'magma_threshold_comparison.tsv', sep='\t', index=False)

# Check GSEA inputs include all tested genes by comparing per-pair tested genes to pathway files and method notes/readme.
pathway_dir = RESULTS / 'step8_pathway'
readme_paths = [RESULTS / 'step8_magma' / 'README.md', RESULTS / 'step8_magma' / 'method_notes.md']
notes = []
for p in readme_paths:
    if p.exists():
        notes.append(p.read_text(errors='ignore'))
notes_text = '\n'.join(notes).lower()
path_files = [
    pathway_dir / 'pair_specific_gsea_go_bp.tsv',
    pathway_dir / 'pair_specific_gsea_kegg.tsv',
    pathway_dir / 'pair_specific_gsea_reactome.tsv',
]
pathway_pair_counts = {}
for p in path_files:
    df = filter_exclude(read_tsv(p))
    pathway_pair_counts[p.name] = df['pair_id'].nunique() if not df.empty and 'pair_id' in df.columns else 0

bonf_by_pair = summary[['pair_id','bonferroni_sig_genes']].copy()
bonf_total_rows = int((genes['_bonf'] <= 0.05).sum())
bonf_pairs = int((summary['bonferroni_sig_genes'] > 0).sum())
top_bonf = summary.loc[summary['bonferroni_sig_genes'] > 0, ['pair_id','trait1','trait2','bonferroni_sig_genes']].head(20)

print('PAIR_GENE_COUNT_SUMMARY', OUTDIR / 'pair_gene_count_summary.tsv', len(summary))
print('MAGMA_THRESHOLD_COMPARISON', OUTDIR / 'magma_threshold_comparison.tsv', len(threshold))
print('TOTAL_TESTED_GENE_ROWS', len(genes))
print('TOTAL_BONFERRONI_SIG_GENE_ROWS', bonf_total_rows)
print('PAIRS_WITH_BONFERRONI_SIG_GENES', bonf_pairs)
print('GSEA_PAIR_COUNTS', pathway_pair_counts)
print('GSEA_NOTES_MENTION_ALL_TESTED', ('all magma tested genes' in notes_text) or ('all tested genes' in notes_text))
print('TOP_BONFERRONI_PAIR_CONTRIBUTORS')
print(top_bonf.to_csv(sep='\t', index=False))
