#!/usr/bin/env python3
from pathlib import Path
import csv, gzip, math, os, subprocess, sys
import numpy as np
import pandas as pd
from scipy.stats import norm, combine_pvalues, hypergeom
from statsmodels.stats.multitest import multipletests

project_root = os.environ.get("PROJECT_ROOT", "").strip()
if not project_root:
    raise SystemExit("Set PROJECT_ROOT to the analysis project root.")
ROOT = Path(project_root)
STEP8=ROOT/'results/phase0_extension/step8_magma_gene_pathway'
RES=STEP8/'resources'
MAGMA=RES/'magma_official/magma'
GENELOC=RES/'magma_official/NCBI37.3.gene.loc'
PAIR_LIST=STEP8/'pair_ids_step8.tsv'
TASK=int(os.environ.get('SLURM_ARRAY_TASK_ID','1'))
OUTDIR=STEP8/'per_pair'
OUTDIR.mkdir(parents=True, exist_ok=True)

def read(path):
    p=Path(path)
    if not p.exists() or p.stat().st_size==0: return pd.DataFrame()
    return pd.read_csv(p, sep='\t', dtype=str, low_memory=False)

def fnum(x):
    try:
        if pd.isna(x) or str(x).strip()=='' or str(x).lower() in {'na','nan','none'}: return np.nan
        return float(x)
    except Exception: return np.nan

def bh(vals):
    arr=np.array([fnum(x) for x in vals], dtype=float); out=np.full(len(arr), np.nan); mask=np.isfinite(arr)
    if mask.sum(): out[mask]=multipletests(arr[mask], method='fdr_bh')[1]
    return out

def write(df,path): df.to_csv(path, sep='\t', index=False, na_rep='')

def parse_pairs():
    pairs=[]
    with PAIR_LIST.open() as fh:
        next(fh)
        for line in fh:
            if line.strip(): pairs.append(line.strip().split('\t')[0])
    if TASK<1 or TASK>len(pairs): raise SystemExit(f'TASK {TASK} out of range {len(pairs)}')
    return pairs[TASK-1]

def collect_pair_variants(pair_id):
    rows=[]
    high=read(ROOT/'results/phase0_extension/step3_susie_coloc/high_confidence_shared_variants.tsv')
    if not high.empty:
        high=high[high['pair_id']==pair_id]
        for _,r in high.iterrows():
            p=min([x for x in [fnum(r.get('PP.H4')), fnum(r.get('PIP_trait1')), fnum(r.get('PIP_trait2'))] if np.isfinite(x)] or [np.nan])
            pval=max(1e-300, 1-p) if np.isfinite(p) else 1e-6
            rows.append(dict(SNP=r.get('SNP'),CHR=r.get('CHR'),BP=r.get('BP'),pair_id=pair_id,P=pval,source='SuSiE_high_confidence',SuSiE_support='TRUE',PLACO_support='FALSE',CPASSOC_support='FALSE',locus_id=r.get('locus_id','')))
    sig=read(ROOT/'results/phase0_extension/step5_placo_cpassoc/cross_trait_significant_snps.tsv')
    if not sig.empty:
        sig=sig[sig['pair_id']==pair_id]
        for _,r in sig.iterrows():
            method=str(r.get('method',''))
            rows.append(dict(SNP=r.get('SNP'),CHR=r.get('CHR'),BP=r.get('BP'),pair_id=pair_id,P=r.get('p_value'),source='PLACO_significant' if method=='PLACO' else 'CPASSOC_significant',SuSiE_support='FALSE',PLACO_support=str(method=='PLACO').upper(),CPASSOC_support=str(method!='PLACO').upper(),locus_id=''))
    cand=read(ROOT/'results/phase0_extension/step5_placo_cpassoc/candidate_snps.tsv')
    if not cand.empty:
        cand=cand[cand['pair_id']==pair_id]
        for _,r in cand.iterrows():
            src=str(r.get('source',''))
            if src in {'PLACO_significant_SNP','CPASSOC_significant_SNP'}:
                ps=[fnum(r.get('P_trait1')), fnum(r.get('P_trait2'))]
                ps=[x for x in ps if np.isfinite(x) and x>0]
                rows.append(dict(SNP=r.get('SNP'),CHR=r.get('CHR'),BP=r.get('BP'),pair_id=pair_id,P=min(ps) if ps else 1e-6,source=src.replace('_SNP',''),SuSiE_support='FALSE',PLACO_support=str(src.startswith('PLACO')).upper(),CPASSOC_support=str(src.startswith('CPASSOC')).upper(),locus_id=r.get('locus_id','')))
    df=pd.DataFrame(rows)
    if df.empty: return df
    df=df[df['SNP'].notna() & df['CHR'].notna() & df['BP'].notna()]
    df['CHR']=df['CHR'].astype(str).str.replace('chr','',regex=False)
    df['BP']=pd.to_numeric(df['BP'], errors='coerce')
    df=df[df['BP'].notna()]
    agg=[]
    for (snp,chr_,bp),g in df.groupby(['SNP','CHR','BP']):
        ps=[fnum(x) for x in g['P']]; ps=[x for x in ps if np.isfinite(x) and x>0]
        agg.append(dict(SNP=snp,CHR=chr_,BP=int(bp),pair_id=pair_id,P=min(ps) if ps else 1e-6,source=';'.join(sorted(set(g['source']))),SuSiE_support=str((g['SuSiE_support']=='TRUE').any()).upper(),PLACO_support=str((g['PLACO_support']=='TRUE').any()).upper(),CPASSOC_support=str((g['CPASSOC_support']=='TRUE').any()).upper(),locus_id=';'.join(sorted(set(str(x) for x in g.get('locus_id',[]).dropna() if str(x))))))
    return pd.DataFrame(agg)

def load_genes():
    rows=[]
    with GENELOC.open() as fh:
        for line in fh:
            p=line.strip().split()
            if len(p)<6: continue
            rows.append(dict(gene_id=p[0], CHR=str(p[1]), START=int(p[2]), END=int(p[3]), strand=p[4], gene=p[5]))
    return pd.DataFrame(rows)

def map_variants(vars, genes):
    out=[]; bychr={c:g for c,g in genes.groupby('CHR')}
    for _,r in vars.iterrows():
        gs=bychr.get(str(r.CHR)); bp=int(r.BP)
        if gs is None or gs.empty: continue
        inside=gs[(gs.START<=bp)&(gs.END>=bp)]
        if not inside.empty:
            for _,g in inside.iterrows(): out.append({**r.to_dict(),'nearest_gene':g.gene,'distance_to_gene':0,'gene_chr':g.CHR,'gene_start':g.START,'gene_end':g.END})
        else:
            dist=np.minimum((gs.START-bp).abs(), (gs.END-bp).abs()); idx=dist.idxmin(); g=gs.loc[idx]
            out.append({**r.to_dict(),'nearest_gene':g.gene,'distance_to_gene':int(dist.loc[idx]),'gene_chr':g.CHR,'gene_start':g.START,'gene_end':g.END})
    return pd.DataFrame(out)

def run_magma(pair_id, vars):
    qc=[]
    out_prefix=OUTDIR/f'{pair_id}.magma'
    annot_prefix=OUTDIR/f'{pair_id}.annot'
    pfile=OUTDIR/f'{pair_id}.snps_p.tsv'
    snploc=OUTDIR/f'{pair_id}.snps_loc.tsv'
    bfile = str(RES / "plink_merged/1000G.EUR.QC")
    vars[['SNP','P']].drop_duplicates('SNP').to_csv(pfile, sep='\t', index=False)
    vars[['SNP','CHR','BP']].drop_duplicates('SNP').to_csv(snploc, sep='\t', index=False, header=False)
    if not Path(str(MAGMA)).exists() or not os.access(MAGMA, os.X_OK):
        qc.append(('MAGMA_binary','FAIL','MAGMA executable missing or not executable'))
        return pd.DataFrame(), qc
    if not Path(bfile+'.bed').exists():
        qc.append(('MAGMA_reference','FAIL','merged whole-genome PLINK reference missing'))
        return pd.DataFrame(), qc
    try:
        ann_cmd=[str(MAGMA),'--annotate','--snp-loc',str(snploc),'--gene-loc',str(GENELOC),'--out',str(annot_prefix)]
        ann=subprocess.run(ann_cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, timeout=1800)
        qc.append(('MAGMA_annotate','PASS' if ann.returncode==0 else 'FAIL', ann.stdout[-500:].replace('\n',' ')))
        if ann.returncode!=0:
            return pd.DataFrame(), qc
        gene_annot=Path(str(annot_prefix)+'.genes.annot')
        gene_cmd=[str(MAGMA),'--bfile',bfile,'--pval',str(pfile),'N=50000','--gene-annot',str(gene_annot),'--out',str(out_prefix)]
        r=subprocess.run(gene_cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, timeout=1800)
        qc.append(('MAGMA_gene_analysis','PASS' if r.returncode==0 else 'FAIL', r.stdout[-500:].replace('\n',' ')))
    except Exception as e:
        qc.append(('MAGMA_gene_analysis','FAIL',str(e)))
    genes_out=Path(str(out_prefix)+'.genes.out')
    if genes_out.exists():
        try:
            df=pd.read_csv(genes_out, delim_whitespace=True, dtype=str)
            return df, qc
        except Exception as e:
            qc.append(('MAGMA_gene_parse','FAIL',str(e)))
    return pd.DataFrame(), qc

def fallback_gene(mapped):
    cols=['gene','CHR','START','END','N_SNP','ZSTAT','P','FDR']
    if mapped is None or mapped.empty or 'nearest_gene' not in mapped.columns:
        return pd.DataFrame(columns=cols)
    rows=[]
    for gene,g in mapped.groupby('nearest_gene'):
        ps=[fnum(x) for x in g.P]; ps=[max(min(x,1),1e-300) for x in ps if np.isfinite(x) and x>0]
        if not ps: continue
        stat,p=combine_pvalues(ps, method='fisher') if len(ps)>1 else (-2*math.log(ps[0]), ps[0])
        rows.append(dict(gene=gene, CHR=g.gene_chr.iloc[0], START=int(g.gene_start.min()), END=int(g.gene_end.max()), N_SNP=g.SNP.nunique(), ZSTAT=norm.isf(max(min(p,1-1e-16),1e-300)), P=p))
    df=pd.DataFrame(rows).sort_values('P') if rows else pd.DataFrame(columns=cols[:-1])
    df['FDR']=bh(df.P) if len(df) else []
    return df[cols]

def parse_gmt(path):
    sets=[]
    with open(path, errors='ignore') as fh:
        for line in fh:
            p=line.rstrip('\n').split('\t')
            if len(p)>=3: sets.append((p[0],set(x for x in p[2:] if x)))
    return sets

def enrich(pair_id, gene_df):
    bg=set(gene_df.gene) if len(gene_df) else set(); sig=set(gene_df.loc[pd.to_numeric(gene_df.FDR,errors='coerce')<0.05,'gene'])
    if not sig: sig=set(gene_df.head(min(50,len(gene_df))).gene)
    rows=[]
    for db,name in [('GO','GO_Biological_Process_2023.gmt'),('KEGG','KEGG_2021_Human.gmt'),('Reactome','Reactome_2022.gmt')]:
        path = RES / "gmt" / name
        if not path.exists(): continue
        for pathway,genes in parse_gmt(path):
            genes=genes&bg; M=len(bg); N=len(sig); K=len(genes); x=len(genes&sig)
            if M and N and K and x:
                rows.append(dict(pair_id=pair_id,database=db,pathway=pathway,n_gene=x,beta=(x/N)/(K/M),P=hypergeom.sf(x-1,M,K,N)))
    df=pd.DataFrame(rows)
    if df.empty: return pd.DataFrame(columns=['pair_id','database','pathway','n_gene','beta','P','FDR'])
    df['FDR']=bh(df.P)
    return df.sort_values(['pair_id','FDR','P'])[['pair_id','database','pathway','n_gene','beta','P','FDR']]

def main():
    pair=parse_pairs(); qc=[]
    vars=collect_pair_variants(pair)
    if vars.empty:
        qc.append((pair,'input','FAIL','no variants for pair'))
    genes=load_genes()
    mapped=map_variants(vars, genes) if not vars.empty else pd.DataFrame()
    if mapped.empty:
        qc.append((pair,'variant_mapping','WARN','no SNP mapped to NCBI37 gene interval/nearest gene'))
    mg, qcm=run_magma(pair, vars) if not vars.empty else (pd.DataFrame(), [])
    qc += [(pair,item,status,reason) for item,status,reason in qcm]
    # Use MAGMA output if it parsed; otherwise fallback but annotate QC.
    if not mg.empty and all(c in mg.columns for c in ['GENE','CHR','START','STOP','NSNPS','ZSTAT','P']):
        id_to_symbol=dict(zip(genes['gene_id'].astype(str), genes['gene'].astype(str)))
        gene_symbol=mg.GENE.astype(str).map(id_to_symbol).fillna(mg.GENE.astype(str))
        out=pd.DataFrame({'pair_id':pair,'gene':gene_symbol,'CHR':mg.CHR,'START':mg.START,'END':mg.STOP,'N_SNP':mg.NSNPS,'ZSTAT':mg.ZSTAT,'P':mg.P})
        out['FDR']=bh(out.P)
    else:
        out=fallback_gene(mapped); out.insert(0,'pair_id',pair)
        qc.append((pair,'MAGMA_fallback','WARN','official MAGMA not run for this pair; fallback SNP-to-gene aggregation emitted with same fields'))
    vgm=mapped[['SNP','CHR','BP','nearest_gene','distance_to_gene','source','SuSiE_support','PLACO_support','CPASSOC_support']].drop_duplicates() if not mapped.empty else pd.DataFrame(columns=['SNP','CHR','BP','nearest_gene','distance_to_gene','source','SuSiE_support','PLACO_support','CPASSOC_support'])
    vgm.insert(0,'pair_id',pair)
    gp=[]
    for gene,g in mapped.groupby('nearest_gene') if not mapped.empty else []:
        ev=[]
        if (g.SuSiE_support=='TRUE').any(): ev.append('SuSiE')
        if (g.PLACO_support=='TRUE').any(): ev.append('PLACO')
        if (g.CPASSOC_support=='TRUE').any(): ev.append('CPASSOC')
        score=3*('SuSiE'in ev)+2*('PLACO'in ev)+2*('CPASSOC'in ev)+math.log1p(g.SNP.nunique())
        gp.append(dict(pair_id=pair,gene=gene,evidence_type=';'.join(ev),n_supporting_variants=g.SNP.nunique(),n_supporting_pairs=1,priority_score=score))
    gp=pd.DataFrame(gp) if gp else pd.DataFrame(columns=['pair_id','gene','evidence_type','n_supporting_variants','n_supporting_pairs','priority_score'])
    kg=[]
    for _,r in gp.iterrows():
        ev=set(str(r.evidence_type).split(';'))
        tier='Tier1' if ('SuSiE'in ev and ('PLACO'in ev or 'CPASSOC'in ev)) else ('Tier2' if 'SuSiE'in ev else 'Tier3')
        g=mapped[mapped.nearest_gene==r.gene]
        kg.append(dict(pair_id=pair,gene=r.gene,n_pairs=1,n_loci=len(set(x for x in g.locus_id.astype(str) if x)),n_variants=g.SNP.nunique(),support_evidence=r.evidence_type,priority_level=tier))
    kg=pd.DataFrame(kg) if kg else pd.DataFrame(columns=['pair_id','gene','n_pairs','n_loci','n_variants','support_evidence','priority_level'])
    enr=enrich(pair,out)
    fuma=vars[['SNP','CHR','BP','P','source']].drop_duplicates('SNP') if not vars.empty else pd.DataFrame(columns=['SNP','CHR','BP','P','source'])
    for df,suffix in [(out,'magma_gene_results'),(vgm,'variant_gene_mapping'),(gp,'gene_priority'),(kg,'key_shared_genes'),(enr,'gene_set_enrichment'),(fuma,'fuma_input_snps')]:
        write(df, OUTDIR/f'{pair}.{suffix}.tsv')
    pd.DataFrame(qc, columns=['pair_id','item','status','reason']).to_csv(OUTDIR/f'{pair}.magma_qc_report.tsv',sep='\t',index=False)
    print(pair, 'genes', len(out), 'snps', len(fuma), 'pathways', len(enr), 'qc', len(qc))
if __name__=='__main__': main()
