#!/usr/bin/env python3
import csv, gzip, glob, os, re, sys
from collections import defaultdict, OrderedDict

BASE = os.environ.get("PROJECT_ROOT", "").strip()
if not BASE:
    raise SystemExit("Set PROJECT_ROOT to the analysis project root.")
STEP1 = f"{BASE}/results/phase0_extension/step0_1"
STEP1_LDSC = f"{BASE}/results/phase0_extension/step1_ldsc_screen"
STEP2_LAVA = f"{BASE}/results/phase0_extension/step2_lava_screen"
PKG = f"{BASE}/results/phase0_server/nasal4_vs_other_package"
OUT = f"{BASE}/results/phase0_extension/step2_input_tables"
os.makedirs(OUT, exist_ok=True)
os.makedirs(f"{OUT}/checks", exist_ok=True)

PAIR_MANIFEST_OUT = f"{OUT}/pair_manifest.tsv"
COLOC_POS_OUT = f"{OUT}/coloc_positive_loci.tsv"
CANDIDATE_OUT = f"{OUT}/candidate_snps.tsv"
PRIORITY_OUT = f"{OUT}/priority_pairs.tsv"
MISSING_OUT = f"{OUT}/input_table_missing_fields_report.tsv"
VALIDATION_OUT = f"{OUT}/checks/step2_validation.tsv"

missing_rows = []

def add_missing(table_name, pair_id='', trait1='', trait2='', locus_id='', snp='', missing_field='', missing_source_file='', problem_type='', note=''):
    missing_rows.append({
        'table_name': table_name, 'pair_id': pair_id, 'trait1': trait1, 'trait2': trait2,
        'locus_id': locus_id, 'SNP': snp, 'missing_field': missing_field,
        'missing_source_file': missing_source_file, 'problem_type': problem_type, 'note': note
    })

def read_tsv(path):
    if not path or not os.path.exists(path):
        return []
    with open(path, newline='') as fh:
        return list(csv.DictReader(fh, delimiter='\t'))

def write_tsv(path, fields, rows):
    with open(path, 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter='\t', extrasaction='ignore')
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, '') for k in fields})

def abs_path(p):
    if not p or p in ('NA', 'nan'):
        return ''
    if os.path.isabs(p):
        return p
    return os.path.normpath(os.path.join(BASE, p))

def split_pair(pair_id):
    if '__' in pair_id:
        a, b = pair_id.split('__', 1)
        return a, b
    return pair_id, ''

def trait_keys(trait):
    keys = []
    if not trait:
        return keys
    keys += [trait, trait.upper(), trait.lower()]
    if '.' in trait:
        last = trait.split('.')[-1]
        keys += [last, last.upper(), last.lower()]
    return list(dict.fromkeys(keys))

def to_float(x):
    try:
        if x is None or x == '' or str(x).upper() == 'NA': return None
        return float(x)
    except Exception:
        return None

# Step 1 manifest path maps.
manifest = read_tsv(f"{STEP1}/manifests/data_manifest.tsv")
std_map = {}
munged_map = {}
manifest_by_type = defaultdict(list)
for r in manifest:
    manifest_by_type[r.get('data_type','')].append(r)
    trait = r.get('trait','')
    dtype = r.get('data_type','')
    path = r.get('symlink_path') or r.get('original_path') or ''
    if dtype == 'standardized_gwas':
        for k in trait_keys(trait): std_map[k] = path
    elif dtype == 'ldsc_sumstats':
        for k in trait_keys(trait): munged_map[k] = path

# Backfill files that were generated before Step 1 but not symlinked there.
# These remain source paths under PROJECT_ROOT; no large files are copied.
for path in glob.glob(f"{BASE}/results/phase0_server/nasal4_package/common/*.common.tsv.gz"):
    trait = os.path.basename(path).replace('.common.tsv.gz','')
    for k in trait_keys(trait):
        std_map.setdefault(k, path)
for path in glob.glob(f"{BASE}/results/phase0_server/nasal4_package/sumstats/*.sumstats.gz"):
    trait = os.path.basename(path).replace('.sumstats.gz','')
    for k in trait_keys(trait):
        munged_map.setdefault(k, path)

def lookup(mapping, trait):
    for k in trait_keys(trait):
        if k in mapping:
            return mapping[k]
    return ''

# Existing package manifests/results.
base_pair_path = f"{PKG}/nasal4_vs_other_pair_manifest.tsv"
mtag_pair_path = f"{PKG}/nasal4_vs_other_mtag_pair_manifest.tsv"
ldsc_path = f"{PKG}/nasal4_vs_other_ldsc_rg_summary.tsv"
lava_path = f"{PKG}/figure_local_rg_edges.tsv"
coloc_manifest_path = f"{PKG}/phase0_v3_coloc_locus_manifest.tsv"
candidate_map_path = f"{PKG}/figure_candidate_variant_gene_map.tsv"
mtag_figure_path = f"{PKG}/figure_mtag_manhattan.tsv"

base_pairs = read_tsv(base_pair_path)
mtag_pairs = read_tsv(mtag_pair_path)
mtag_pair_ids = set(r.get('pair_id','') for r in mtag_pairs if r.get('pair_id',''))
ldsc_rows = read_tsv(ldsc_path)
lava_rows = read_tsv(lava_path)
ldsc_bonf_rows = read_tsv(f"{STEP1_LDSC}/ldsc_rg_significant.tsv")
lava_bonf_rows = read_tsv(f"{STEP2_LAVA}/lava_local_rg_significant.tsv")
coloc_loci = read_tsv(coloc_manifest_path)
candidate_map = read_tsv(candidate_map_path)
mtag_fig = read_tsv(mtag_figure_path)

pairs = OrderedDict()
for source_rows in (base_pairs, mtag_pairs):
    for r in source_rows:
        pid = r.get('pair_id','')
        if not pid: continue
        pairs.setdefault(pid, {}).update({k:r.get(k,'') for k in r})
for r in ldsc_rows:
    pid = r.get('pair_id','')
    if not pid: continue
    d = pairs.setdefault(pid, {})
    d.setdefault('trait_a', r.get('trait_a',''))
    d.setdefault('trait_b', r.get('trait_b',''))
for r in lava_rows:
    pid = r.get('trait_pair','')
    if not pid: continue
    a,b = split_pair(pid)
    d = pairs.setdefault(pid, {})
    d.setdefault('trait_a', a); d.setdefault('trait_b', b)
for r in coloc_loci:
    pid = r.get('pair_id','')
    if not pid: continue
    a,b = split_pair(pid)
    d = pairs.setdefault(pid, {})
    d.setdefault('trait_a', a); d.setdefault('trait_b', b)

ldsc_by_pair = {r.get('pair_id',''): r for r in ldsc_rows if r.get('pair_id','')}
lava_pairs = set(r.get('trait_pair','') for r in lava_rows if r.get('trait_pair',''))
coloc_by_pair = defaultdict(list)
for r in coloc_loci:
    if r.get('pair_id'):
        coloc_by_pair[r['pair_id']].append(r)

# Coloc summaries from individual result files.
coloc_summary = {}
for f in glob.glob(f"{PKG}/coloc_merged/*.coloc.tsv"):
    rows = read_tsv(f)
    if not rows: continue
    r = rows[0]
    lid = r.get('locus_id','')
    pid = r.get('pair_id','')
    if lid and pid:
        r = dict(r)
        r['_file'] = f
        coloc_summary[(pid,lid)] = r

# Pair manifest.
pair_manifest_rows = []
for pid, r in pairs.items():
    t1 = r.get('trait_a') or split_pair(pid)[0]
    t2 = r.get('trait_b') or split_pair(pid)[1]
    std1 = lookup(std_map, t1) or abs_path(r.get('sumstats_a',''))
    std2 = lookup(std_map, t2) or abs_path(r.get('sumstats_b',''))
    if not std1:
        add_missing('pair_manifest', pid, t1, t2, missing_field='trait1_standard_gwas', missing_source_file=f"{STEP1}/manifests/data_manifest.tsv;{base_pair_path};nasal4_package/common", problem_type='no_standardized_gwas_path', note='No standardized GWAS path found in Step 1 manifest or existing result manifests.')
    elif not os.path.exists(std1):
        add_missing('pair_manifest', pid, t1, t2, missing_field='trait1_standard_gwas', missing_source_file=std1, problem_type='file_not_found', note='Path was found in an existing manifest but the file is absent.')
        std1 = ''
    if not std2:
        add_missing('pair_manifest', pid, t1, t2, missing_field='trait2_standard_gwas', missing_source_file=f"{STEP1}/manifests/data_manifest.tsv;{base_pair_path};nasal4_package/common", problem_type='no_standardized_gwas_path', note='No standardized GWAS path found in Step 1 manifest or existing result manifests.')
    elif not os.path.exists(std2):
        add_missing('pair_manifest', pid, t1, t2, missing_field='trait2_standard_gwas', missing_source_file=std2, problem_type='file_not_found', note='Path was found in an existing manifest but the file is absent.')
        std2 = ''
    munged1 = lookup(munged_map, t1) or abs_path(r.get('ldsc_sumstats_a',''))
    munged2 = lookup(munged_map, t2) or abs_path(r.get('ldsc_sumstats_b',''))
    if not munged1:
        add_missing('pair_manifest', pid, t1, t2, missing_field='trait1_munged', missing_source_file=f"{STEP1}/manifests/data_manifest.tsv;{base_pair_path}", problem_type='no_munged_path')
    if not munged2:
        add_missing('pair_manifest', pid, t1, t2, missing_field='trait2_munged', missing_source_file=f"{STEP1}/manifests/data_manifest.tsv;{base_pair_path}", problem_type='no_munged_path')
    mtag_dir = abs_path(r.get('mtag_out_dir','')) or f"{PKG}/mtag/{pid}"
    t1_mtag = os.path.join(mtag_dir, f"{pid}_trait_1.txt")
    t2_mtag = os.path.join(mtag_dir, f"{pid}_trait_2.txt")
    if not os.path.exists(t1_mtag):
        ptype = 'not_selected_for_mtag_ldsc_screen' if pid not in mtag_pair_ids else 'file_not_found'
        note = 'MTAG was only run for LDSC-positive/BH-screened pairs.' if pid not in mtag_pair_ids else ''
        add_missing('pair_manifest', pid, t1, t2, missing_field='trait1_mtag', missing_source_file=t1_mtag, problem_type=ptype, note=note)
        t1_mtag = ''
    if not os.path.exists(t2_mtag):
        ptype = 'not_selected_for_mtag_ldsc_screen' if pid not in mtag_pair_ids else 'file_not_found'
        note = 'MTAG was only run for LDSC-positive/BH-screened pairs.' if pid not in mtag_pair_ids else ''
        add_missing('pair_manifest', pid, t1, t2, missing_field='trait2_mtag', missing_source_file=t2_mtag, problem_type=ptype, note=note)
        t2_mtag = ''
    ldsc_res = ''
    if pid in ldsc_by_pair:
        ldsc_res = abs_path(ldsc_by_pair[pid].get('log','')) or ldsc_path
    else:
        add_missing('pair_manifest', pid, t1, t2, missing_field='ldsc_result', missing_source_file=ldsc_path, problem_type='pair_not_found')
    lava_res = lava_path if pid in lava_pairs else ''
    if not lava_res:
        add_missing('pair_manifest', pid, t1, t2, missing_field='lava_result', missing_source_file=lava_path, problem_type='pair_not_found')
    coloc_paths = []
    for loc in coloc_by_pair.get(pid, []):
        cp = abs_path(loc.get('coloc_out_file',''))
        if cp and os.path.exists(cp): coloc_paths.append(cp)
    coloc_res = ';'.join(coloc_paths)
    if not coloc_res:
        ptype = 'not_selected_for_coloc_candidate_loci' if pid not in coloc_by_pair else 'no_coloc_file_found'
        add_missing('pair_manifest', pid, t1, t2, missing_field='coloc_result', missing_source_file=coloc_manifest_path, problem_type=ptype)
    pair_manifest_rows.append({
        'pair_id': pid, 'trait1': t1, 'trait2': t2,
        'trait1_standard_gwas': std1, 'trait2_standard_gwas': std2,
        'trait1_munged': munged1, 'trait2_munged': munged2,
        'trait1_mtag': t1_mtag, 'trait2_mtag': t2_mtag,
        'ldsc_result': ldsc_res, 'lava_result': lava_res, 'coloc_result': coloc_res
    })

# Locus maps.
locus_by_key = {}
locus_by_region = {}
locus_by_lead = {}
for loc in coloc_loci:
    pid = loc.get('pair_id',''); lid = loc.get('locus_id','')
    if not pid or not lid: continue
    locus_by_key[(pid,lid)] = loc
    region = f"{loc.get('chr','')}:{loc.get('start','')}-{loc.get('end','')}"
    locus_by_region[(pid, region)] = loc
    lead = loc.get('lead_rsid','')
    if lead:
        locus_by_lead[(pid, lead)] = loc

# coloc_positive_loci.
coloc_pos_rows = []
for (pid,lid), cs in sorted(coloc_summary.items()):
    pp4 = to_float(cs.get('PP.H4.abf'))
    if pp4 is None or pp4 < 0.5:
        continue
    loc = locus_by_key.get((pid,lid), {})
    t1,t2 = split_pair(pid)
    pp3 = cs.get('PP.H3.abf','')
    row = {
        'pair_id': pid, 'trait1': t1, 'trait2': t2, 'locus_id': lid,
        'CHR': loc.get('chr',''), 'START': loc.get('start',''), 'END': loc.get('end',''),
        'PP.H3': pp3, 'PP.H4': cs.get('PP.H4.abf',''),
        'analysis_group': 'main' if pp4 >= 0.8 else 'secondary'
    }
    for fld in ('CHR','START','END'):
        if not row[fld]: add_missing('coloc_positive_loci', pid, t1, t2, lid, missing_field=fld, missing_source_file=coloc_manifest_path, problem_type='missing_locus_coordinate')
    coloc_pos_rows.append(row)

# Helpers for candidate SNP effect extraction from region files.
def open_maybe_gz(path):
    if path.endswith('.gz'):
        return gzip.open(path, 'rt', newline='')
    return open(path, newline='')

def read_snp_from_region(path, snp):
    path = abs_path(path)
    if not path or not os.path.exists(path) or not snp:
        return None
    try:
        with open_maybe_gz(path) as fh:
            reader = csv.DictReader(fh, delimiter='\t')
            for row in reader:
                if row.get('SNP') == snp or row.get('rsid') == snp:
                    return row
    except Exception:
        return None
    return None

def best_row_from_region(path):
    path = abs_path(path)
    if not path or not os.path.exists(path):
        return None
    best = None
    best_p = None
    try:
        with open_maybe_gz(path) as fh:
            reader = csv.DictReader(fh, delimiter='\t')
            for row in reader:
                snp = row.get('SNP') or row.get('rsid') or ''
                if not snp:
                    continue
                pval = to_float(row.get('P') or row.get('p'))
                if pval is None:
                    continue
                if best_p is None or pval < best_p:
                    best = row
                    best_p = pval
    except Exception:
        return None
    return best

def read_region_by_snp(path):
    path = abs_path(path)
    rows = {}
    if not path or not os.path.exists(path):
        return rows
    try:
        with open_maybe_gz(path) as fh:
            reader = csv.DictReader(fh, delimiter='\t')
            for row in reader:
                snp = row.get('SNP') or row.get('rsid') or ''
                if snp:
                    rows[snp] = row
    except Exception:
        return {}
    return rows

def infer_locus_lead_snp(loc):
    a_rows = read_region_by_snp(loc.get('region_a_file',''))
    b_rows = read_region_by_snp(loc.get('region_b_file',''))
    shared = set(a_rows).intersection(b_rows)
    if shared:
        def score(snp):
            pa = to_float(a_rows[snp].get('P') or a_rows[snp].get('p'))
            pb = to_float(b_rows[snp].get('P') or b_rows[snp].get('p'))
            vals = [v for v in (pa, pb) if v is not None]
            return min(vals) if vals else float('inf')
        return min(shared, key=score)
    candidates = []
    for key in ('region_a_file', 'region_b_file'):
        row = best_row_from_region(loc.get(key,''))
        if row:
            candidates.append(row)
    if not candidates:
        return ''
    candidates.sort(key=lambda r: to_float(r.get('P') or r.get('p')) if to_float(r.get('P') or r.get('p')) is not None else float('inf'))
    return candidates[0].get('SNP') or candidates[0].get('rsid') or ''

def candidate_row_from_locus(pid, lid, snp, source):
    loc = locus_by_key.get((pid,lid), {})
    t1,t2 = split_pair(pid)
    a = read_snp_from_region(loc.get('region_a_file',''), snp)
    b = read_snp_from_region(loc.get('region_b_file',''), snp)
    if a is None:
        add_missing('candidate_snps', pid, t1, t2, lid, snp, 'trait1_region_row', abs_path(loc.get('region_a_file','')), 'snp_not_found_or_region_missing')
        a = {}
    if b is None:
        add_missing('candidate_snps', pid, t1, t2, lid, snp, 'trait2_region_row', abs_path(loc.get('region_b_file','')), 'snp_not_found_or_region_missing')
        b = {}
    ea = a.get('A1') or a.get('EA') or b.get('A1') or b.get('EA') or ''
    oa = a.get('A2') or a.get('OA') or b.get('A2') or b.get('OA') or ''
    chr_ = loc.get('chr') or a.get('CHR') or b.get('CHR') or ''
    bp = loc.get('lead_pos') if snp == loc.get('lead_rsid') else (a.get('BP') or b.get('BP') or '')
    out = {
        'pair_id': pid, 'trait1': t1, 'trait2': t2, 'locus_id': lid, 'SNP': snp,
        'CHR': chr_, 'BP': bp, 'EA': ea, 'OA': oa,
        'BETA_trait1': a.get('BETA',''), 'SE_trait1': a.get('SE',''), 'P_trait1': a.get('P',''),
        'BETA_trait2': b.get('BETA',''), 'SE_trait2': b.get('SE',''), 'P_trait2': b.get('P',''),
        'source': source
    }
    for fld in ('SNP','CHR','BP','EA','OA'):
        if not out.get(fld): add_missing('candidate_snps', pid, t1, t2, lid, snp, fld, coloc_manifest_path, 'missing_required_snp_field', f'source={source}')
    return out

candidate_rows = []
seen_candidate = set()
def add_candidate(row):
    key = (row.get('pair_id',''), row.get('locus_id',''), row.get('SNP',''), row.get('source',''))
    if key in seen_candidate: return
    seen_candidate.add(key)
    candidate_rows.append(row)

# coloc H4 lead SNP candidates.
for row in coloc_pos_rows:
    pid = row['pair_id']; lid = row['locus_id']
    loc = locus_by_key.get((pid,lid), {})
    snp = loc.get('lead_rsid','') or infer_locus_lead_snp(loc)
    if snp:
        add_candidate(candidate_row_from_locus(pid, lid, snp, 'coloc_H4_lead_SNP'))
    else:
        add_missing('candidate_snps', pid, row['trait1'], row['trait2'], lid, '', 'SNP', coloc_manifest_path, 'missing_lead_rsid', 'coloc_H4_lead_SNP')

# MTAG and LAVA initial candidates from existing candidate map.
for r in candidate_map:
    pid = r.get('trait_pair','')
    region = r.get('region','')
    snp = r.get('rsid','')
    loc = locus_by_region.get((pid, region)) or locus_by_lead.get((pid, snp))
    lid = loc.get('locus_id','') if loc else ''
    t1,t2 = split_pair(pid)
    if not lid:
        add_missing('candidate_snps', pid, t1, t2, '', snp, 'locus_id', candidate_map_path, 'cannot_map_region_or_lead_snp_to_coloc_locus', region)
        continue
    if not snp:
        snp = infer_locus_lead_snp(loc)
        if not snp:
            add_missing('candidate_snps', pid, t1, t2, lid, '', 'SNP', candidate_map_path, 'missing_lava_region_representative_snp', region)
            continue
    sources = r.get('candidate_sources','')
    if 'mtag' in sources.lower():
        add_candidate(candidate_row_from_locus(pid, lid, snp, 'MTAG_novel_SNP'))
    n_lava = to_float(r.get('n_lava_regions')) or 0
    if 'lava' in sources.lower() or n_lava > 0:
        add_candidate(candidate_row_from_locus(pid, lid, snp, 'LAVA_recurrent_locus_lead_SNP'))

# Priority pairs.
ldsc_sig = set()
ldsc_stats = {}
if ldsc_bonf_rows:
    # The prefiltered table is produced by the LDSC BH-FDR screening step.
    for r in ldsc_bonf_rows:
        pid = r.get('pair_id','')
        if not pid: continue
        ldsc_sig.add(pid)
        ldsc_stats[pid] = {'rg': r.get('rg',''), 'p': r.get('p','')}
else:
    for r in ldsc_rows:
        pid = r.get('pair_id','')
        if not pid: continue
        ldsc_stats[pid] = {'rg': r.get('rg',''), 'p': r.get('p','')}
    for r in mtag_pairs:
        pid = r.get('pair_id','')
        if not pid: continue
        ldsc_stats[pid] = {'rg': r.get('rg',''), 'p': r.get('rg_p','')}
        q = to_float(r.get('rg_bh_q'))
        if q is not None and q < 0.05:
            ldsc_sig.add(pid)

lava_pos_counts = defaultdict(int)
if lava_bonf_rows:
    # The prefiltered table contains LAVA rows with p_adj <= 0.05.
    for r in lava_bonf_rows:
        pid = r.get('pair_id') or r.get('trait_pair','')
        if pid:
            lava_pos_counts[pid] += 1
else:
    for r in lava_rows:
        pid = r.get('trait_pair','')
        p = to_float(r.get('p'))
        padj = to_float(r.get('p_adj'))
        sig = str(r.get('significant','')).upper() == 'TRUE'
        if sig or (padj is not None and padj < 0.05):
            lava_pos_counts[pid] += 1

coloc_h4_counts = defaultdict(int)
for r in coloc_pos_rows:
    coloc_h4_counts[r['pair_id']] += 1

mtag_loci_counts = defaultdict(set)
for r in candidate_map:
    if 'mtag' in r.get('candidate_sources','').lower():
        pid = r.get('trait_pair','')
        reg = r.get('region','')
        if pid and reg: mtag_loci_counts[pid].add(reg)

priority_rows = []
for pid in sorted(pairs.keys()):
    t1,t2 = split_pair(pid)
    reasons = []
    if pid in ldsc_sig: reasons.append('LDSC_significant')
    if lava_pos_counts.get(pid,0) >= 2: reasons.append('LAVA_multi_locus_positive')
    if coloc_h4_counts.get(pid,0) > 0: reasons.append('coloc_H4_positive')
    if len(mtag_loci_counts.get(pid,set())) >= 5: reasons.append('MTAG_many_novel_loci')
    # No explicit clinical-priority source file exists; do not guess clinical priority.
    if reasons:
        priority_rows.append({
            'pair_id': pid, 'trait1': t1, 'trait2': t2,
            'priority_reason': ';'.join(reasons),
            'LDSC_rg': ldsc_stats.get(pid,{}).get('rg',''),
            'LDSC_p': ldsc_stats.get(pid,{}).get('p',''),
            'n_LAVA_loci': str(lava_pos_counts.get(pid,0)),
            'n_coloc_H4_loci': str(coloc_h4_counts.get(pid,0)),
            'n_MTAG_novel_loci': str(len(mtag_loci_counts.get(pid,set())))
        })
add_missing('priority_pairs', missing_field='clinical_priority_definition', missing_source_file='not_found_under_current_inputs', problem_type='no_clinical_priority_list', note='No explicit clinical priority pair list was found; clinical_priority was not guessed.')

# Output files.
write_tsv(PAIR_MANIFEST_OUT, ['pair_id','trait1','trait2','trait1_standard_gwas','trait2_standard_gwas','trait1_munged','trait2_munged','trait1_mtag','trait2_mtag','ldsc_result','lava_result','coloc_result'], pair_manifest_rows)
write_tsv(COLOC_POS_OUT, ['pair_id','trait1','trait2','locus_id','CHR','START','END','PP.H3','PP.H4','analysis_group'], coloc_pos_rows)
write_tsv(CANDIDATE_OUT, ['pair_id','trait1','trait2','locus_id','SNP','CHR','BP','EA','OA','BETA_trait1','SE_trait1','P_trait1','BETA_trait2','SE_trait2','P_trait2','source'], candidate_rows)
write_tsv(PRIORITY_OUT, ['pair_id','trait1','trait2','priority_reason','LDSC_rg','LDSC_p','n_LAVA_loci','n_coloc_H4_loci','n_MTAG_novel_loci'], priority_rows)
write_tsv(MISSING_OUT, ['table_name','pair_id','trait1','trait2','locus_id','SNP','missing_field','missing_source_file','problem_type','note'], missing_rows)

# Validation: required columns and row counts.
required = {
    'pair_manifest.tsv': ['pair_id','trait1','trait2','trait1_standard_gwas','trait2_standard_gwas','trait1_munged','trait2_munged','trait1_mtag','trait2_mtag','ldsc_result','lava_result','coloc_result'],
    'coloc_positive_loci.tsv': ['pair_id','trait1','trait2','locus_id','CHR','START','END','PP.H3','PP.H4','analysis_group'],
    'candidate_snps.tsv': ['pair_id','trait1','trait2','locus_id','SNP','CHR','BP','EA','OA','BETA_trait1','SE_trait1','P_trait1','BETA_trait2','SE_trait2','P_trait2','source'],
    'priority_pairs.tsv': ['pair_id','trait1','trait2','priority_reason','LDSC_rg','LDSC_p','n_LAVA_loci','n_coloc_H4_loci','n_MTAG_novel_loci'],
    'input_table_missing_fields_report.tsv': ['table_name','pair_id','trait1','trait2','locus_id','SNP','missing_field','missing_source_file','problem_type','note']
}
with open(VALIDATION_OUT, 'w', newline='') as fh:
    w = csv.writer(fh, delimiter='\t')
    w.writerow(['table','rows','required_columns_complete','missing_columns'])
    for name, cols in required.items():
        path = f"{OUT}/{name}"
        rows = read_tsv(path)
        with open(path, newline='') as tfh:
            header = next(csv.reader(tfh, delimiter='\t'))
        missing_cols = [c for c in cols if c not in header]
        w.writerow([name, len(rows), 'TRUE' if not missing_cols else 'FALSE', ','.join(missing_cols)])
