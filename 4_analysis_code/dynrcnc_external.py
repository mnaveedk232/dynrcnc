#!/usr/bin/env python3
"""
DynRCNC EXTERNAL VALIDATION
===========================
Applies the published DynRCNC model, unchanged, to proteins that were not used at any stage
of rule design or threshold calibration, and prints the results in the same layout as the
benchmark run (dynrcnc.py / run_clean_benchmark.py).

Nothing is refitted here. The model code, the six criteria and all eighteen numerical
thresholds are imported from the benchmark copy of dynrcnc.py; only the input paths and the
pair table differ. MD features are computed with the same functions, from 200 ns trajectories
sampled at 100 ps and corrected for periodic boundaries, exactly as for the benchmark.

PATHS (edit these four if the folders are elsewhere)
  MODEL_DIR  ./model
               dynrcnc.py                     the published model, imported unchanged
               benchmark_265pairs.ddg         required by dynrcnc.py at import time
  EXT_DIR    ./data
               <pdb>_md/<pdb>.N               RING 4.0 nodes, default settings
               <pdb>_md/<pdb>.E               RING 4.0 edges
               <pdb>_md/<pdb>_md_first_frame.gro
               <pdb>_md/<pdb>_md_reduced.xtc  200 ns, 2001 frames, PBC corrected
  DDG_FILE   ./data/external_final.ddg
               columns: PDB, Mutation_Double, dddG, pH, Method
  FEAT_FILE  ./data/features_265.csv
               benchmark feature table, used only to TRAIN the ML baselines

OUTPUT (written to OUT_DIR = ./results/)
  external_all_predictions.csv   one row per pair, same columns as clean_all_predictions.csv
  external_summary.csv           per-protein and combined metrics
  external_baselines.csv         DynRCNC against static RCNC, all-additive and three ML models

RUN
  cd external_validation_4proteins
  python3 dynrcnc_external.py
Requires: pandas, numpy, mdtraj, networkx, scipy, scikit-learn
"""
import os, sys, time, warnings, importlib.util
warnings.filterwarnings('ignore')
import numpy as np, pandas as pd, mdtraj as md

# ───────────────────────────────── paths
MODEL_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'model')
EXT_DIR   = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data')
DDG_FILE  = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data', 'external_final.ddg')
FEAT_FILE = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data', 'features_265.csv')
OUT_DIR   = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'results')

# external proteins, none of which appear in the benchmark
EXTERNAL = {
    '1JQZ': 'FGF-1',
    '1C9O': 'Cold shock protein',
    '1MJC': 'Cold shock protein A (E. coli)',
    '1OIA': 'U1A RNP domain',
}

def cfg(pdb):
    p = pdb.lower(); d = os.path.join(EXT_DIR, f'{p}_md')
    return {'node_file': os.path.join(d, f'{p}.N'),
            'edge_file': os.path.join(d, f'{p}.E'),
            'gro':       os.path.join(d, f'{p}_md_first_frame.gro'),
            'xtc':       os.path.join(d, f'{p}_md_reduced.xtc'),
            'fps':       100}

# ───────────────────────────────── import the published model
sys.argv = ['dynrcnc.py', '--base', MODEL_DIR,
            '--ddg', os.path.join(MODEL_DIR, 'benchmark_265pairs.ddg')]
spec = importlib.util.spec_from_file_location('dynrcnc', os.path.join(MODEL_DIR, 'dynrcnc.py'))
dm = importlib.util.module_from_spec(spec)
spec.loader.exec_module(dm)

def metrics(df):
    """same definition as compute_metrics() in dynrcnc.py"""
    TP = int(((df.exp == 'Non-additive') & (df.pred == 'Non-additive')).sum())
    TN = int(((df.exp == 'Additive')     & (df.pred == 'Additive')).sum())
    FP = int(((df.exp == 'Additive')     & (df.pred == 'Non-additive')).sum())
    FN = int(((df.exp == 'Non-additive') & (df.pred == 'Additive')).sum())
    N, NA, ADD = TP+TN+FP+FN, TP+FN, TN+FP
    se = TP/NA if NA else 0.0; sp = TN/ADD if ADD else 0.0
    pr = TP/(TP+FP) if TP+FP else 0.0
    f1 = 2*pr*se/(pr+se) if pr+se else 0.0
    den = np.sqrt(float((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN)))
    return dict(N=N, NA=NA, ADD=ADD, TP=TP, TN=TN, FP=FP, FN=FN,
                Sensitivity=se, Specificity=sp, Precision=pr, F1=f1,
                MCC=(TP*TN-FP*FN)/den if den else 0.0)

def bootstrap_ci(df, key, n=1000, seed=42):
    rng = np.random.default_rng(seed); v = []
    for _ in range(n):
        s = df.sample(len(df), replace=True, random_state=int(rng.integers(1e9)))
        v.append(metrics(s)[key])
    return float(np.percentile(v, 2.5)), float(np.percentile(v, 97.5))

def run_protein(pdb, c, ddg_df):
    """same layout as run_protein() in dynrcnc.py"""
    print(f"\n{'─'*55}\n  {pdb}  [{EXTERNAL[pdb]}]\n{'─'*55}")
    communities, r2c, direct, G_nb, bc, linked_comm, last = dm.build_network(c['node_file'], c['edge_file'])
    print(f"  Communities: {len(communities)}")
    resids, rmap, dccm, z, rmsf, ss = dm.compute_md_features(c['gro'], c['xtc'], c['fps'])

    t = md.load(c['xtc'], top=c['gro'])
    t = t.atom_slice(t.topology.select('name CA')); t.superpose(t, 0)
    r = md.rmsd(t, t, 0)*10.0; ns = (t.time-t.time[0])/1000.0
    print(f"  Trajectory: {t.n_frames} frames, {ns[-1]:.0f} ns, "
          f"Ca-RMSD last 100 ns = {r[ns>=100].mean():.2f} A, max = {r.max():.2f} A")

    sub = ddg_df[ddg_df['PDB'] == pdb].copy()
    if len(sub) == 0:
        print("  No DDG data."); return None
    rows = []
    for _, row in sub.iterrows():
        parts = str(row['Mutation_Double']).split(',')
        nums = [int(''.join(ch for ch in p if ch.isdigit())) for p in parts if any(ch.isdigit() for ch in p)]
        if len(nums) != 2:
            print(f"  [WARN] {pdb}: skipping malformed mutation string '{row['Mutation_Double']}'")
            continue
        s1, s2 = nums
        if s1 not in rmap or s2 not in rmap:
            print(f"  [WARN] {pdb}: residue(s) {[x for x in (s1,s2) if x not in rmap]} for "
                  f"'{row['Mutation_Double']}' not in MD structure; Z-DCCM treated as 0")
        dddg = float(row['dddG'])
        exp = 'Non-additive' if abs(dddg) >= dm.DDDG_CUT else 'Additive'
        pred, mech = dm.predict(s1, s2, row['Mutation_Double'], communities, r2c, G_nb, bc,
                                direct, linked_comm, last, rmap, ss, dccm, z, rmsf)
        i_, j_ = rmap.get(s1, -1), rmap.get(s2, -1)
        rows.append({'PDB': pdb, 'mutation': row['Mutation_Double'],
                     'pH': row.get('pH', ''), 'Method': row.get('Method', ''),
                     'dddG': dddg, 'exp': exp, 'pred': pred, 'correct': exp == pred,
                     'mechanism': mech, 'seq_dist': abs(s1-s2), 'res_i': s1, 'res_j': s2,
                     'comm_i': r2c.get(s1, -1), 'comm_j': r2c.get(s2, -1),
                     'zdccm': round(abs(z[i_, j_]), 4) if i_ >= 0 and j_ >= 0 else 0,
                     'ss_i': dm.get_ss(s1, rmap, ss), 'ss_j': dm.get_ss(s2, rmap, ss)})
    df = pd.DataFrame(rows); m = metrics(df)
    print(f"\n  RESULTS: {pdb}")
    print(f"  TP={m['TP']} TN={m['TN']} FP={m['FP']} FN={m['FN']} N={m['N']}")
    print(f"  Sensitivity : {m['Sensitivity']*100:.1f}%  ({m['TP']}/{m['NA']})")
    print(f"  Specificity : {m['Specificity']*100:.1f}%  ({m['TN']}/{m['ADD']})")
    print(f"  Precision   : {m['Precision']*100:.1f}%")
    print(f"  F1          : {m['F1']:.3f}")
    print(f"  MCC         : {m['MCC']:.3f}")
    cols = ['mutation','pH','Method','dddG','exp','pred','correct',
            'mechanism','seq_dist','comm_i','comm_j','zdccm','ss_i','ss_j']
    print(f"\n  ── Non-additive pairs ({m['NA']}) ──")
    print(df[df['exp'] == 'Non-additive'][cols].to_string(index=False))
    print(f"\n  ── Additive pairs ({m['ADD']}) ──")
    print(df[df['exp'] == 'Additive'][cols].to_string(index=False))
    print(f"\n  Mechanism breakdown:")
    print(df[df['pred'] == 'Non-additive']['mechanism'].value_counts().to_string())
    return df, m

def baselines(all_df):
    """Table 6 baselines on the same external pairs. The ML models are trained on the
       265-pair benchmark and only applied here: no external data is used for fitting."""
    from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier
    from sklearn.linear_model import LogisticRegression
    from sklearn.preprocessing import StandardScaler
    FEAT = ['seq_dist','same_community','adj_comm','comm_linked','in_direct','z_dccm','raw_dccm',
            'bc_max','bc_min','rmsf_z_max','rmsf_z_min','n_common_neigh','both_coil','both_helix',
            'both_struct','both_isolated','one_isolated','both_charged_kept','vol_change_any']
    charged = lambda m: dm.is_charged_orig(m) and dm.is_charged_mut(m)
    feats, static = [], []
    for pdb in EXTERNAL:
        c = cfg(pdb)
        communities, r2c, direct, G_nb, bc, linked, last = dm.build_network(c['node_file'], c['edge_file'])
        res, rmap, dccm, z, rmsf, ss = dm.compute_md_features(c['gro'], c['xtc'], c['fps'])
        rz = (rmsf - rmsf.mean())/(rmsf.std() if rmsf.std() else 1)
        for x in all_df[all_df.PDB == pdb].itertuples():
            a, b = int(x.res_i), int(x.res_j); parts = x.mutation.split(',')
            i, j = rmap.get(a, -1), rmap.get(b, -1); ca, cb = r2c.get(a, -1), r2c.get(b, -1)
            s1, s2 = dm.get_ss(a, rmap, ss), dm.get_ss(b, rmap, ss)
            n1 = set(G_nb.neighbors(a)) if a in G_nb else set()
            n2 = set(G_nb.neighbors(b)) if b in G_nb else set()
            same = int(ca != -1 and ca == cb); hel = int(s1 == 'H' and s2 == 'H')
            static.append(bool(same and not hel))
            feats.append(dict(seq_dist=abs(a-b), same_community=same,
                adj_comm=int(ca != -1 and cb != -1 and ca != cb),
                comm_linked=int(ca != -1 and cb != -1 and ca != cb and cb in linked.get(ca, set())),
                in_direct=int((min(a,b), max(a,b)) in direct),
                z_dccm=abs(z[i,j]) if i >= 0 and j >= 0 else 0.0,
                raw_dccm=abs(dccm[i,j]) if i >= 0 and j >= 0 else 0.0,
                bc_max=max(bc.get(a,0), bc.get(b,0)), bc_min=min(bc.get(a,0), bc.get(b,0)),
                rmsf_z_max=max(rz[i] if i >= 0 else 0, rz[j] if j >= 0 else 0),
                rmsf_z_min=min(rz[i] if i >= 0 else 0, rz[j] if j >= 0 else 0),
                n_common_neigh=len(n1 & n2), both_coil=int(s1 == 'C' and s2 == 'C'),
                both_helix=hel, both_struct=int(s1 in 'HE' and s2 in 'HE'),
                both_isolated=int(ca == -1 and cb == -1), one_isolated=int(ca == -1 or cb == -1),
                both_charged_kept=int(charged(parts[0]) and charged(parts[1])),
                vol_change_any=int(dm.has_major_vol_change(parts[0]) or dm.has_major_vol_change(parts[1]))))
    tr = pd.read_csv(FEAT_FILE)
    sc = StandardScaler().fit(tr[FEAT].values.astype(float))
    Xtr, ytr = sc.transform(tr[FEAT].values.astype(float)), tr.label.values
    Xe = sc.transform(pd.DataFrame(feats)[FEAT].values.astype(float))
    def row(name, pred):
        d = all_df.copy(); d['pred'] = np.where(pred, 'Non-additive', 'Additive')
        m = metrics(d)
        pw = float(np.mean([metrics(d[d.PDB == k])['MCC'] for k in EXTERNAL]))
        return dict(model=name, TP=m['TP'], TN=m['TN'], FP=m['FP'], FN=m['FN'],
                    Sens=round(m['Sensitivity']*100,1), Spec=round(m['Specificity']*100,1),
                    MCC=round(m['MCC'],3), protein_weighted_MCC=round(pw,3))
    R = [row('DynRCNC', (all_df.pred == 'Non-additive').values),
         row('Static RCNC', np.array(static)),
         row('All-additive reference', np.zeros(len(all_df), bool))]
    for nm, mdl in [('Random forest', RandomForestClassifier(n_estimators=300, class_weight='balanced',
                                                             random_state=42, n_jobs=-1)),
                    ('Gradient boosting', GradientBoostingClassifier(n_estimators=200, max_depth=3, random_state=42)),
                    ('Logistic regression', LogisticRegression(class_weight='balanced', max_iter=2000, random_state=42))]:
        mdl.fit(Xtr, ytr)
        R.append(row(nm + ' (trained on benchmark)', mdl.predict(Xe).astype(bool)))
    return pd.DataFrame(R)

def main():
    t0 = time.time()
    os.makedirs(OUT_DIR, exist_ok=True)
    print(f"\n{'='*62}")
    print(f"  DynRCNC  EXTERNAL VALIDATION")
    print(f"  14 active rules | 3 tiers | {len(EXTERNAL)} external proteins")
    print(f"  model: {os.path.join(MODEL_DIR,'dynrcnc.py')} (unchanged, nothing refitted)")
    print(f"  class boundary: |dddG| >= {dm.DDDG_CUT} kcal/mol")
    print(f"{'='*62}")

    for pdb in EXTERNAL:
        for k, v in cfg(pdb).items():
            if k != 'fps' and not os.path.exists(v): sys.exit(f"\n  missing input: {v}")
    if not os.path.exists(DDG_FILE): sys.exit(f"\n  missing pair table: {DDG_FILE}")
    ddg_df = pd.read_csv(DDG_FILE, sep='\t', on_bad_lines='skip')
    ddg_df.columns = ddg_df.columns.str.strip()
    print(f"\nDDG loaded: {len(ddg_df)} pairs  [{DDG_FILE}]")

    protein_dfs, summary = {}, []
    for pdb in EXTERNAL:
        out = run_protein(pdb, cfg(pdb), ddg_df)
        if out: protein_dfs[pdb], m = out[0], out[1]; summary.append({'PDB': pdb, **m})
    if not protein_dfs: print("\nNo proteins processed."); return

    all_df = pd.concat(protein_dfs.values(), ignore_index=True)
    comb = metrics(all_df); summary.append({'PDB': 'COMBINED', **comb})
    print(f"\n  Computing bootstrap CIs...")
    ci = {k: bootstrap_ci(all_df, k) for k in ['MCC','F1','Sensitivity','Specificity']}
    yt = (all_df.exp == 'Non-additive').astype(int).values
    yp = (all_df.pred == 'Non-additive').astype(int).values
    p_perm = dm.permutation_test(yt, yp)
    y_static = np.array([1 if (a != -1 and b != -1 and a == b) else 0
                         for a, b in zip(all_df.comm_i.values, all_df.comm_j.values)])
    p_mcn, b_mcn, c_mcn = dm.mcnemar_test(yt, yp, y_static)

    print(f"\n{'='*62}")
    print(f"  COMBINED: {len(protein_dfs)} external proteins, {len(all_df)} pairs")
    print(f"{'='*62}")
    print(f"  TP={comb['TP']} TN={comb['TN']} FP={comb['FP']} FN={comb['FN']}")
    for met in ['Sensitivity','Specificity','Precision','F1','MCC']:
        lo, hi = ci.get(met, (None, None))
        s = f"  95%CI [{lo:.3f}–{hi:.3f}]" if lo is not None else ""
        print(f"  {met:14}: {comb[met]*100:.2f}%{s}")
    print(f"\n  Statistical Tests:")
    print(f"  Permutation p    : {p_perm:.4f} {'✓' if p_perm < 0.05 else '✗'}")
    print(f"  McNemar p        : {p_mcn:.4f}  [vs static RCNC clique baseline; b={b_mcn}, c={c_mcn}]")

    sdf = pd.DataFrame(summary)
    print(f"\n{'─'*78}")
    print(f"  {'PDB':8} {'N':>4} {'NA':>4} {'ADD':>4} {'Sens%':>7} {'Spec%':>7} {'F1':>6} {'MCC':>6}")
    print(f"{'─'*78}")
    for _, r in sdf.iterrows():
        print(f"  {r['PDB']:8} {r['N']:>4} {r['NA']:>4} {r['ADD']:>4} "
              f"{r['Sensitivity']*100:>7.1f} {r['Specificity']*100:>7.1f} {r['F1']:>6.3f} {r['MCC']:>6.3f}")
    pw = sdf[sdf.PDB != 'COMBINED'].MCC.mean()
    print(f"\n  Protein-weighted mean MCC: {pw:.3f}")

    print(f"\n  BASELINES ON THE SAME EXTERNAL PAIRS")
    print(f"  (ML models trained on the 265-pair benchmark; no external data used for fitting)")
    b = baselines(all_df)
    print(f"{'─'*98}")
    print(f"  {'Model':42} {'TP':>4} {'TN':>4} {'FP':>4} {'FN':>4} {'Sens%':>7} {'Spec%':>7} {'MCC':>7} {'wMCC':>7}")
    print(f"{'─'*98}")
    for _, r in b.iterrows():
        print(f"  {r['model']:42} {r['TP']:>4} {r['TN']:>4} {r['FP']:>4} {r['FN']:>4} "
              f"{r['Sens']:>7.1f} {r['Spec']:>7.1f} {r['MCC']:>7.3f} {r['protein_weighted_MCC']:>7.3f}")

    all_df.to_csv(os.path.join(OUT_DIR, 'external_all_predictions.csv'), index=False)
    sdf.to_csv(os.path.join(OUT_DIR, 'external_summary.csv'), index=False)
    b.to_csv(os.path.join(OUT_DIR, 'external_baselines.csv'), index=False)
    print(f"\n  Results saved → {os.path.abspath(OUT_DIR)}/")
    print(f"  Runtime: {time.time()-t0:.1f}s")
    print(f"{'='*62}\n")
    print("Wrote: external_all_predictions.csv, external_summary.csv, external_baselines.csv")

if __name__ == '__main__':
    main()
