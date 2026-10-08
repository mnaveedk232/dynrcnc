#!/usr/bin/env python3
"""
ML BASELINE for DynRCNC — CLEAN 265-pair benchmark (logic identical to ml_baseline_271.py).
=================================================================
Generates the 19 numerical features directly from dynrcnc.py's own pipeline
(no external features.csv needed), then trains RandomForest / GradientBoosting /
LogisticRegression and reports in-sample + nested LOPOCV (protein-weighted) MCC.

HOW TO RUN:
  1. Put this file in the SAME folder as dynrcnc.py
  2. python3 ml_baseline_271.py
     (or:  DYNRCNC_BASE=/path/to/data python3 ml_baseline_271.py)
  3. It writes  ml_baseline_271_results.csv   <-- give this file to Claude.

Notes:
  - Uses fixed seed (42) everywhere for reproducibility.
  - Features are computed by replicating dynrcnc.py's predict() feature block
    exactly, so they match the model's own criteria inputs.
"""
import warnings; warnings.filterwarnings('ignore')
import types, os, sys
import numpy as np, pandas as pd
from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler

SEED = 42
np.random.seed(SEED)

# ── 1. Load dynrcnc.py as a module (functions only, skip its __main__) ──
HERE = os.path.dirname(os.path.abspath(__file__))
CLEAN = os.path.expanduser(os.environ.get('CLEAN_DIR', '~/Desktop/revision_analysis/00_clean_benchmark'))
sys.argv = ['dynrcnc.py', '--base', CLEAN, '--ddg', os.path.join(CLEAN, 'benchmark_271pairs.ddg')]
MODEL_PATH = os.path.join(CLEAN, 'dynrcnc.py')
src = open(MODEL_PATH).read().split("if __name__")[0]
mod = types.ModuleType("dyn"); mod.__dict__['__name__'] = 'dyn'
exec(compile(src, 'dynrcnc.py', 'exec'), mod.__dict__)

PROTEINS = mod.PROTEINS
DDG_FILE = mod.DDG_FILE
DDDG_CUT = mod.DDDG_CUT
get_ss = mod.get_ss
is_charged_orig = mod.is_charged_orig
is_charged_mut  = mod.is_charged_mut
has_major_vol_change = mod.has_major_vol_change

ddg = pd.read_csv(DDG_FILE, sep='\t', on_bad_lines='skip')
ddg.columns = ddg.columns.str.strip()

# ── 2. Build the 19-feature table for all 271 pairs ──
FEATS = ['seq_dist','same_community','adj_comm','comm_linked','in_direct',
         'z_dccm','raw_dccm','bc_max','bc_min','rmsf_z_max','rmsf_z_min',
         'n_common_neigh','both_coil','both_helix','both_struct',
         'both_isolated','one_isolated','both_charged_kept','vol_change_any']

print("Building features (1QJP absent, so it is skipped) (this computes MD features per protein)...")
rows = []
for pdb, cfg in PROTEINS.items():
    miss = [k for k in ['node_file','edge_file','gro','xtc'] if not os.path.exists(cfg[k])]
    if miss:
        print(f"  SKIP {pdb}: missing {miss}"); continue
    comm, r2c, direct, G_nb, bc, linked_comm, last = mod.build_network(cfg['node_file'], cfg['edge_file'])
    resids, rmap, dccm, z, rmsf, ss = mod.compute_md_features(cfg['gro'], cfg['xtc'], cfg['fps'])
    rmsf_mean = rmsf.mean(); rmsf_std = max(rmsf.std(), 0.001)
    for _, row in ddg[ddg.PDB == pdb].iterrows():
        parts = str(row['Mutation_Double']).split(',')
        nums = [int(''.join(c for c in p if c.isdigit())) for p in parts if any(c.isdigit() for c in p)]
        if len(nums) != 2: continue
        s1, s2 = nums
        seq = abs(s1 - s2)
        c1 = r2c.get(s1, -1); c2 = r2c.get(s2, -1)
        same = (c1 != -1 and c2 != -1 and c1 == c2)
        ss1 = get_ss(s1, rmap, ss); ss2 = get_ss(s2, rmap, ss)
        pair = (min(s1, s2), max(s1, s2))
        i = rmap.get(s1, -1); j = rmap.get(s2, -1)
        zv = abs(z[i, j])    if i >= 0 and j >= 0 else 0.0
        dc = abs(dccm[i, j]) if i >= 0 and j >= 0 else 0.0
        rf1 = rmsf[i] if i >= 0 else 0.0
        rf2 = rmsf[j] if j >= 0 else 0.0
        bc1 = bc.get(s1, 0); bc2 = bc.get(s2, 0)
        n1 = set(G_nb.neighbors(s1)) if s1 in G_nb else set()
        n2 = set(G_nb.neighbors(s2)) if s2 in G_nb else set()
        common = n1 & n2
        adj_comm = (c1 != -1 and c2 != -1 and c1 != c2)
        both_coil = (ss1 == 'C' and ss2 == 'C')
        both_helix = (ss1 == 'H' and ss2 == 'H')
        both_struc = (ss1 in ('H', 'E') and ss2 in ('H', 'E'))
        ch1 = is_charged_orig(parts[0]) and is_charged_mut(parts[0]) if parts else False
        ch2 = is_charged_orig(parts[1]) and is_charged_mut(parts[1]) if len(parts) > 1 else False
        vol1 = has_major_vol_change(parts[0]) if parts else False
        vol2 = has_major_vol_change(parts[1]) if len(parts) > 1 else False
        vol_sig = vol1 or vol2
        rf1_z = (rf1 - rmsf_mean) / rmsf_std; rf2_z = (rf2 - rmsf_mean) / rmsf_std
        both_iso = (c1 == -1 and c2 == -1)
        one_iso  = (c1 == -1 or c2 == -1)
        comm_linked = (c1 != -1 and c2 != -1 and c1 != c2 and c2 in linked_comm.get(c1, set()))
        rows.append(dict(
            PDB=pdb, mutation=row['Mutation_Double'],
            label=1 if abs(float(row['dddG'])) >= DDDG_CUT else 0,
            seq_dist=seq, same_community=int(same), adj_comm=int(adj_comm),
            comm_linked=int(comm_linked), in_direct=int(pair in direct),
            z_dccm=zv, raw_dccm=dc, bc_max=max(bc1, bc2), bc_min=min(bc1, bc2),
            rmsf_z_max=max(rf1_z, rf2_z), rmsf_z_min=min(rf1_z, rf2_z),
            n_common_neigh=len(common), both_coil=int(both_coil),
            both_helix=int(both_helix), both_struct=int(both_struc),
            both_isolated=int(both_iso), one_isolated=int(one_iso),
            both_charged_kept=int(ch1 and ch2), vol_change_any=int(vol_sig),
        ))

df = pd.DataFrame(rows)
df.to_csv(os.path.join(HERE, 'features_265.csv'), index=False)
print(f"  features_265.csv written: {len(df)} pairs  (NA={df.label.sum()} ADD={(df.label==0).sum()})")

X = df[FEATS].values.astype(float)
y = df['label'].values.astype(int)
groups = df['PDB'].values

# ── 3. Metrics + models ──
def metrics(yt, yp):
    yt = np.asarray(yt); yp = np.asarray(yp)
    tp = int(((yt==1)&(yp==1)).sum()); tn = int(((yt==0)&(yp==0)).sum())
    fp = int(((yt==0)&(yp==1)).sum()); fn = int(((yt==1)&(yp==0)).sum())
    sens = tp/(tp+fn) if tp+fn else 0.0
    spec = tn/(tn+fp) if tn+fp else 0.0
    den = np.sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn))
    mcc = (tp*tn-fp*fn)/den if den > 0 else 0.0
    return mcc, sens, spec

def make(name):
    if name == 'RandomForest':
        return RandomForestClassifier(n_estimators=300, class_weight='balanced',
                                      random_state=SEED, n_jobs=-1)
    if name == 'GradientBoosting':
        return GradientBoostingClassifier(n_estimators=200, max_depth=3, random_state=SEED)
    if name == 'LogisticRegression':
        return LogisticRegression(class_weight='balanced', max_iter=2000, random_state=SEED)

print('\n' + '='*70)
print(f'  ML BASELINE on DynRCNC features  ({len(df)} pairs, {len(FEATS)} features, seed={SEED})')
print(f'  NA={y.sum()}  ADD={(y==0).sum()}')
print('='*70)

out = []
for name in ['RandomForest', 'GradientBoosting', 'LogisticRegression']:
    # in-sample (resubstitution)
    sc = StandardScaler().fit(X)
    clf = make(name).fit(sc.transform(X), y)
    mcc_in, se_in, sp_in = metrics(y, clf.predict(sc.transform(X)))
    # nested LOPOCV (protein-weighted)
    yp_cv = np.zeros_like(y); per = {}
    for held in np.unique(groups):
        tr = groups != held; te = groups == held
        s = StandardScaler().fit(X[tr])
        c = make(name).fit(s.transform(X[tr]), y[tr])
        yp_cv[te] = c.predict(s.transform(X[te]))
        m, _, _ = metrics(y[te], yp_cv[te])
        per[held] = (m, int(te.sum()), int(y[te].sum()))
    w_mcc = sum(per[h][0]*per[h][1] for h in per)/len(y)
    mcc_pool, se_cv, sp_cv = metrics(y, yp_cv)
    gap = mcc_in - w_mcc
    print(f'\n{name}')
    print(f'  In-sample : MCC={mcc_in:.3f}  Sens={se_in*100:.1f}%  Spec={sp_in*100:.1f}%')
    print(f'  LOPOCV    : weighted MCC={w_mcc:.3f}  pooled MCC={mcc_pool:.3f}  Sens={se_cv*100:.1f}%  Spec={sp_cv*100:.1f}%')
    print(f'  Gap       : {gap:.3f}')
    for h in sorted(per):
        m, n, na = per[h]; print(f'    {h:6} N={n:3d} NA={na:2d}  MCC={m:.3f}')
    out.append(dict(method=name, in_sample_MCC=round(mcc_in,4),
                    LOPOCV_weighted_MCC=round(w_mcc,4), LOPOCV_pooled_MCC=round(mcc_pool,4),
                    gap=round(gap,4), LOPOCV_sens=round(se_cv,4), LOPOCV_spec=round(sp_cv,4)))

res = pd.DataFrame(out)
res.to_csv(os.path.join(HERE, 'ml_baseline_265_results.csv'), index=False)
print('\n' + '='*70)
print('  Saved -> ml_baseline_265_results.csv   AND   features_265.csv')
print('  (Give ml_baseline_271_results.csv to Claude for verification.)')
print('='*70)
