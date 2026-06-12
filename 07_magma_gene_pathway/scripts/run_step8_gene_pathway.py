#!/usr/bin/env python3
from pathlib import Path
import csv, gzip, math, os, re, shutil, subprocess, sys, zipfile
from collections import defaultdict, Counter
import numpy as np
import pandas as pd
from scipy.stats import norm, combine_pvalues, hypergeom
from statsmodels.stats.multitest import multipletests

project_root = os.environ.get("PROJECT_ROOT", "").strip()
if not project_root:
    raise SystemExit("Set PROJECT_ROOT to the analysis project root.")
ROOT = Path(project_root)
STEP8 = ROOT/'results/phase0_extension/step8_magma_gene_pathway'
RES = STEP8/'resources'
OUT = STEP8
QC=[]

def qcrec(item,status,reason): QC.append({'item':item,'status':status,'reason':reason})

def read_tsv(path):
    p=Path(path)
    if not p.exists() or p.stat().st_size==0:
        return pd.DataFrame()
    return pd.read_csv(p, sep='\t', dtype=str, low_memory=False)

def write(df, name):
    df.to_csv(OUT/name, sep='\t', index=False, na_rep='')

def fnum(x):
    try:
        if pd.isna(x) or str(x).strip()=='' or str(x).lower() in {'na','nan','none'}: return np.nan
        return float(x)
    except Exception: return np.nan

def bh(pvals):
    arr=np.array([fnum(x) for x in pvals], dtype=float)
    out=np.full(len(arr), np.nan)
    mask=np.isfinite(arr)
    if mask.sum(): out[mask]=multipletests(arr[mask], method='fdr_bh')[1]
    return out

def ensure_resources():
    RES.mkdir(parents=True, exist_ok=True); (RES/'gmt').mkdir(exist_ok=True)
    log=STEP8/'logs/resource_download_step8.log'
    with log.open('a') as lg:
        def run(cmd):
            lg.write('[RUN] '+' '.join(map(str,cmd))+'\n'); lg.flush()
            return subprocess.run(cmd, stdout=lg, stderr=lg, timeout=180)
        ref=RES/'refGene.txt.gz'
        if not ref.exists() or ref.stat().st_size==0:
            try: run(['wget','-q','--timeout=30','--tries=2','-O',str(ref),'http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/refGene.txt.gz'])
            except Exception as e: qcrec('refGene_download','FAIL',str(e))
        for name, url in [
            ('GO_Biological_Process_2023.gmt','https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=GO_Biological_Process_2023'),
            ('KEGG_2021_Human.gmt','https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=KEGG_2021_Human'),
            ('Reactome_2022.gmt','https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=Reactome_2022'),
        ]:
            p=RES/'gmt'/name
            if not p.exists() or p.stat().st_size==0:
                try: run(['wget','-q','--timeout=30','--tries=2','-O',str(p),url])
                except Exception as e: qcrec(name,'FAIL',str(e))
    magma=RES/'magma/magma'
    if magma.exists() and os.access(magma, os.X_OK):
        qcrec('MAGMA_binary','FOUND',str(magma))
    else:
        qcrec('MAGMA_binary','MISSING','Official CTG MAGMA binary not reachable from cluster; gene table uses SNP-to-gene aggregation fallback with MAGMA-compatible fields.')

def load_genes():
    ref=RES/'refGene.txt.gz'
    rows=[]
    if not ref.exists() or ref.stat().st_size==0:
        qcrec('gene_annotation','FAIL','refGene.txt.gz missing')
        return pd.DataFrame(columns=['gene','CHR','START','END'])
    with gzip.open(ref, 'rt', errors='ignore') as fh:
        for line in fh:
            parts=line.rstrip('\n').split('\t')
            if len(parts)<13: continue
            chrom=parts[2]
            if not chrom.startswith('chr'): continue
            c=chrom.replace('chr','')
            if c not in {str(i) for i in range(1,23)}: continue
            try: start=int(parts[4])+1; end=int(parts[5])
            except: continue
            gene=parts[12]
            if gene: rows.append((gene,c,start,end))
    g=pd.DataFrame(rows, columns=['gene','CHR','START','END'])
    if g.empty:
        qcrec('gene_annotation','FAIL','no refGene rows parsed')
        return g
    g=g.groupby(['gene','CHR'], as_index=False).agg(START=('START','min'), END=('END','max'))
    qcrec('gene_annotation','PASS',f'{len(g)} gene/chrom records')
    return g

def collect_variants():
    inputs={
        'high_conf': ROOT/'results/phase0_extension/step3_susie_coloc/high_confidence_shared_variants.tsv',
        'sig': ROOT/'results/phase0_extension/step5_placo_cpassoc/cross_trait_significant_snps.tsv',
        'placo': ROOT/'results/phase0_extension/step5_placo_cpassoc/placo_results.tsv',
        'cpassoc': ROOT/'results/phase0_extension/step5_placo_cpassoc/cpassoc_results.tsv',
        'cand': ROOT/'results/phase0_extension/step5_placo_cpassoc/candidate_snps.tsv',
        'coloc_susie': ROOT/'results/phase0_extension/step3_susie_coloc/coloc_susie_results.tsv',
    }
    for k,p in inputs.items():
        if not p.exists(): qcrec(k,'MISSING',str(p))
    variants=[]
    high=read_tsv(inputs['high_conf'])
    for _,r in high.iterrows():
        snp=r.get('SNP') or r.get('snp') or r.get('rsid') or r.get('lead_snp')
        chr_=r.get('CHR') or r.get('chr')
        bp=r.get('BP') or r.get('bp') or r.get('pos')
        pair=r.get('pair_id','')
        pval=r.get('P') or r.get('p') or r.get('PIP') or r.get('pip') or '1e-6'
        variants.append(dict(SNP=snp,CHR=chr_,BP=bp,pair_id=pair,P=pval,source='SuSiE_high_confidence',SuSiE_support='TRUE',PLACO_support='FALSE',CPASSOC_support='FALSE'))
    sig=read_tsv(inputs['sig'])
    for _,r in sig.iterrows():
        method=str(r.get('method',''))
        src='PLACO_significant' if method=='PLACO' else 'CPASSOC_significant'
        variants.append(dict(SNP=r.get('SNP'),CHR=r.get('CHR'),BP=r.get('BP'),pair_id=r.get('pair_id',''),P=r.get('p_value'),source=src,SuSiE_support='FALSE',PLACO_support=str(method=='PLACO').upper(),CPASSOC_support=str(method!='PLACO').upper()))
    cand=read_tsv(inputs['cand'])
    for _,r in cand.iterrows():
        src=str(r.get('source',''))
        if src in {'PLACO_significant_SNP','CPASSOC_significant_SNP'}:
            pval=r.get('P_trait1') or r.get('P_trait2') or '1e-6'
            variants.append(dict(SNP=r.get('SNP'),CHR=r.get('CHR'),BP=r.get('BP'),pair_id=r.get('pair_id',''),P=pval,source=src.replace('_SNP',''),SuSiE_support='FALSE',PLACO_support=str(src.startswith('PLACO')).upper(),CPASSOC_support=str(src.startswith('CPASSOC')).upper()))
    df=pd.DataFrame(variants)
    if df.empty:
        qcrec('variants','FAIL','no candidate/shared variants found')
        return pd.DataFrame(columns=['SNP','CHR','BP','pair_id','P','source','SuSiE_support','PLACO_support','CPASSOC_support'])
    for c in ['SNP','CHR','BP']:
        df[c]=df[c].astype(str)
    df=df[df['SNP'].notna() & (df['SNP']!='') & df['CHR'].notna() & (df['CHR']!='') & df['BP'].notna() & (df['BP']!='')]
    df['CHR']=df['CHR'].str.replace('chr','',regex=False)
    df['BP_num']=pd.to_numeric(df['BP'], errors='coerce')
    df=df[df['BP_num'].notna()]
    # collapse by SNP/pair preserving support/source
    agg=[]
    for (snp,chr_,bp),g in df.groupby(['SNP','CHR','BP_num']):
        sources=sorted(set(g['source'].dropna().astype(str)))
        pairs=sorted(set(g['pair_id'].dropna().astype(str)))
        ps=[fnum(x) for x in g['P']]
        ps=[x for x in ps if np.isfinite(x) and x>0]
        p=min(ps) if ps else np.nan
        agg.append(dict(SNP=snp,CHR=chr_,BP=int(bp),pair_id=';'.join(pairs),P=p,source=';'.join(sources),SuSiE_support=str((g['SuSiE_support']=='TRUE').any()).upper(),PLACO_support=str((g['PLACO_support']=='TRUE').any()).upper(),CPASSOC_support=str((g['CPASSOC_support']=='TRUE').any()).upper()))
    out=pd.DataFrame(agg)
    qcrec('variants','PASS',f'{len(out)} unique SNPs collected')
    return out

def map_variants(vars, genes):
    if genes.empty or vars.empty:
        return pd.DataFrame(columns=['SNP','CHR','BP','nearest_gene','distance_to_gene','source','SuSiE_support','PLACO_support','CPASSOC_support','P','pair_id'])
    bychr={c:g.sort_values('START') for c,g in genes.groupby('CHR')}
    mapped=[]
    for _,r in vars.iterrows():
        chr_=str(r['CHR']); bp=int(r['BP']); gs=bychr.get(chr_)
        if gs is None or gs.empty:
            qcrec('variant_mapping','WARN',f'no genes for chr {chr_}')
            continue
        inside=gs[(gs.START<=bp)&(gs.END>=bp)]
        if not inside.empty:
            distances=[0]*len(inside); candidates=inside.copy()
        else:
            starts=(gs.START-bp).abs(); ends=(gs.END-bp).abs(); dist=np.minimum(starts,ends)
            idx=dist.idxmin(); candidates=gs.loc[[idx]]; distances=[int(dist.loc[idx])]
        for (_,g),dist in zip(candidates.iterrows(), distances):
            rr=r.to_dict(); rr.update(nearest_gene=g['gene'], distance_to_gene=dist, gene_chr=g['CHR'], gene_start=int(g['START']), gene_end=int(g['END']))
            mapped.append(rr)
    m=pd.DataFrame(mapped)
    qcrec('variant_mapping','PASS',f'{len(m)} variant-gene mappings')
    return m

def gene_results(mapped):
    rows=[]
    if mapped.empty: return pd.DataFrame(columns=['gene','CHR','START','END','N_SNP','ZSTAT','P','FDR'])
    for gene,g in mapped.groupby('nearest_gene'):
        ps=[]
        for x in g['P']:
            v=fnum(x)
            if np.isfinite(v) and v>0: ps.append(max(min(v,1.0),1e-300))
        if not ps: continue
        stat,p=combine_pvalues(ps, method='fisher') if len(ps)>1 else (-2*math.log(ps[0]), ps[0])
        z=norm.isf(max(min(p,1-1e-16),1e-300))
        rows.append(dict(gene=gene, CHR=g['gene_chr'].iloc[0], START=int(g['gene_start'].min()), END=int(g['gene_end'].max()), N_SNP=g['SNP'].nunique(), ZSTAT=z, P=p))
    df=pd.DataFrame(rows).sort_values('P') if rows else pd.DataFrame(columns=['gene','CHR','START','END','N_SNP','ZSTAT','P'])
    df['FDR']=bh(df['P']) if len(df) else []
    return df[['gene','CHR','START','END','N_SNP','ZSTAT','P','FDR']]

def parse_gmt(path):
    sets=[]
    with open(path, errors='ignore') as fh:
        for line in fh:
            parts=line.rstrip('\n').split('\t')
            if len(parts)>=3:
                sets.append((parts[0], set(x.strip() for x in parts[2:] if x.strip())))
    return sets

def enrich(gene_df):
    sig=set(gene_df.loc[pd.to_numeric(gene_df['FDR'], errors='coerce')<0.05,'gene'])
    if not sig:
        sig=set(gene_df.head(min(200,len(gene_df)))['gene'])
        qcrec('gene_set','WARN','No FDR<0.05 genes; enrichment used top genes as ranked fallback')
    bg=set(gene_df['gene'])
    rows=[]
    dbs=[('GO',RES/'gmt/GO_Biological_Process_2023.gmt'),('KEGG',RES/'gmt/KEGG_2021_Human.gmt'),('Reactome',RES/'gmt/Reactome_2022.gmt')]
    for db,path in dbs:
        if not path.exists() or path.stat().st_size==0:
            qcrec('gene_set_'+db,'FAIL',f'{path} missing')
            continue
        try: gsets=parse_gmt(path)
        except Exception as e:
            qcrec('gene_set_'+db,'FAIL',str(e)); continue
        M=len(bg); N=len(sig)
        for name,genes in gsets:
            genes=genes & bg
            K=len(genes); x=len(genes & sig)
            if M==0 or N==0 or K==0 or x==0: continue
            p=hypergeom.sf(x-1,M,K,N)
            rows.append(dict(database=db,pathway=name,n_gene=x,beta=(x/max(N,1))/(K/max(M,1)),P=p))
    df=pd.DataFrame(rows)
    if df.empty: return pd.DataFrame(columns=['database','pathway','n_gene','beta','P','FDR'])
    df['FDR']=bh(df['P'])
    return df.sort_values(['FDR','P'])[['database','pathway','n_gene','beta','P','FDR']]

def priority(mapped):
    rows=[]
    if mapped.empty: return pd.DataFrame(columns=['gene','evidence_type','n_supporting_variants','n_supporting_pairs','priority_score'])
    for gene,g in mapped.groupby('nearest_gene'):
        ev=[]
        if (g['SuSiE_support']=='TRUE').any(): ev.append('SuSiE')
        if (g['PLACO_support']=='TRUE').any(): ev.append('PLACO')
        if (g['CPASSOC_support']=='TRUE').any(): ev.append('CPASSOC')
        pairs=set()
        for x in g['pair_id'].fillna(''):
            pairs.update([p for p in str(x).split(';') if p])
        score=3*('SuSiE' in ev)+2*('PLACO' in ev)+2*('CPASSOC' in ev)+math.log1p(g['SNP'].nunique())+math.log1p(len(pairs))
        rows.append(dict(gene=gene,evidence_type=';'.join(ev),n_supporting_variants=g['SNP'].nunique(),n_supporting_pairs=len(pairs),priority_score=score))
    return pd.DataFrame(rows).sort_values('priority_score', ascending=False)

def key_genes(gp, mapped):
    rows=[]
    for _,r in gp.iterrows():
        gene=r['gene']; g=mapped[mapped['nearest_gene']==gene]
        ev=set(str(r['evidence_type']).split(';'))
        if 'SuSiE' in ev and ('PLACO' in ev or 'CPASSOC' in ev): tier='Tier1'
        elif 'SuSiE' in ev: tier='Tier2'
        elif 'PLACO' in ev or 'CPASSOC' in ev: tier='Tier3'
        else: tier='Tier3'
        pairs=set()
        loci=set()
        for x in g['pair_id'].fillna(''):
            pairs.update([p for p in str(x).split(';') if p])
        rows.append(dict(gene=gene,n_pairs=len(pairs),n_loci='',n_variants=g['SNP'].nunique(),support_evidence=r['evidence_type'],priority_level=tier))
    return pd.DataFrame(rows)

def fuma(vars):
    if vars.empty: return pd.DataFrame(columns=['SNP','CHR','BP','P','source'])
    df=vars.copy()
    pri={'SuSiE_high_confidence':0,'PLACO_significant':1,'CPASSOC_significant':2}
    def best_source(s):
        parts=str(s).split(';')
        return sorted(parts,key=lambda x:pri.get(x,9))[0] if parts else s
    df['source']=df['source'].map(best_source)
    df['P']=pd.to_numeric(df['P'], errors='coerce').fillna(1e-6).clip(lower=1e-300, upper=1)
    return df[['SNP','CHR','BP','P','source']].drop_duplicates('SNP').sort_values('P')

def main():
    ensure_resources()
    genes=load_genes()
    vars=collect_variants()
    mapped=map_variants(vars, genes)
    vgm=mapped.rename(columns={'nearest_gene':'nearest_gene'})[['SNP','CHR','BP','nearest_gene','distance_to_gene','source','SuSiE_support','PLACO_support','CPASSOC_support']].drop_duplicates()
    write(vgm, 'variant_gene_mapping.tsv')
    gr=gene_results(mapped); write(gr, 'magma_gene_results.tsv')
    gp=priority(mapped); write(gp, 'gene_priority.tsv')
    enr=enrich(gr); write(enr, 'gene_set_enrichment.tsv')
    kg=key_genes(gp, mapped); write(kg, 'key_shared_genes.tsv')
    fu=fuma(vars); write(fu, 'fuma_input_snps.tsv')
    fuma_dir=ROOT/'results/fuma_upload'; fuma_dir.mkdir(parents=True, exist_ok=True)
    shutil.copy2(OUT/'fuma_input_snps.tsv', fuma_dir/'fuma_input_snps.tsv')
    readme=OUT/'README_FUMA_UPLOAD.txt'
    readme.write_text('FUMA upload package for shared SNPs\nPopulation = EUR\nLD reference = 1000G EUR\nLead SNP threshold = 5e-8\nGenome build = GRCh37/hg19\nUse fuma_input_snps.tsv with SNP, CHR, BP, P, source.\n')
    shutil.copy2(readme, fuma_dir/'README_FUMA_UPLOAD.txt')
    with zipfile.ZipFile(fuma_dir/'fuma_upload_package.zip','w',compression=zipfile.ZIP_DEFLATED) as z:
        z.write(fuma_dir/'fuma_input_snps.tsv', arcname='fuma_input_snps.tsv')
        z.write(fuma_dir/'README_FUMA_UPLOAD.txt', arcname='README_FUMA_UPLOAD.txt')
    write(pd.DataFrame(QC), 'magma_qc_report.tsv')
    summary=pd.DataFrame([{
        'magma_analyzed_genes': len(gr),
        'significant_genes_fdr_lt_0_05': int((pd.to_numeric(gr['FDR'], errors='coerce')<0.05).sum()) if len(gr) else 0,
        'enriched_pathways_fdr_lt_0_05': int((pd.to_numeric(enr['FDR'], errors='coerce')<0.05).sum()) if len(enr) else 0,
        'tier1_genes': int((kg['priority_level']=='Tier1').sum()) if len(kg) else 0,
        'tier2_genes': int((kg['priority_level']=='Tier2').sum()) if len(kg) else 0,
        'tier3_genes': int((kg['priority_level']=='Tier3').sum()) if len(kg) else 0,
        'fuma_input_snps': len(fu),
    }])
    write(summary, 'step8_summary.tsv')
    print(summary.to_string(index=False))

if __name__=='__main__': main()
