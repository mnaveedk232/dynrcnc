#!/usr/bin/env python3
"""
Dynamic-gate ablation on the CLEAN 265-pair benchmark (Table S2 + gates-only row
of Table 7 / Section 3.8).

"Opening" a dynamic condition = making that |Z-DCCM| or DCCM test always true,
while every structural pre-condition and every other criterion stays unchanged.
Ten dynamic conditions in total:
  C1 : zv > C1_COUPLING                                  (1)
  C2 : zv > C2_COUPLING                                  (1)
  C3 : coil, helix-one-iso, helix-both-iso, short-helix  (4)
  C4 : dccm > C4_BACKBONE_CORR                           (1)
  C5 : common-neighbour + linked-community cutoffs       (2, "forward cutoffs")
  C5 : hub exclusion gate (zv < 0.15 -> additive)        (1)
  C6 : no dynamic condition

MD features are computed ONCE per protein with dynrcnc.py's own functions,
so the Z-DCCM method is identical to the paper.

Run inside 00_clean_benchmark/ (after run_clean_benchmark.py has made the
<pdb>_md_first_frame.gro / <pdb>_md_reduced.xtc names):
    python3 ablation_dynamic_gates.py
Writes: ablation_dynamic_gates.csv
"""
import os, sys, importlib.util
import numpy as np, pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.argv = ['dynrcnc.py', '--base', HERE,
            '--ddg', os.path.join(HERE, 'benchmark_271pairs.ddg'),
            '--out', os.path.join(HERE, 'clean_out')]
spec = importlib.util.spec_from_file_location('dynrcnc', os.path.join(HERE, 'dynrcnc.py'))
dm = importlib.util.module_from_spec(spec); spec.loader.exec_module(dm)

KEEP = ['1STN', '1BNI', '2LZM', '1PGA', '1CSP', '2RN2', '2CI2']   # 1QJP dropped


def predict_open(s1, s2, mut_str, communities, r2c, G_nb, bc, direct,
                 linked_comm, rmap, ss, dccm, z, op):
    """Copy of dynrcnc.predict() with switchable dynamic conditions.
    op = set of opened gates. Empty set reproduces predict() exactly."""
    seq = abs(s1-s2)
    c1 = r2c.get(s1, -1); c2 = r2c.get(s2, -1)
    same = (c1 != -1 and c2 != -1 and c1 == c2)
    ss1 = dm.get_ss(s1, rmap, ss); ss2 = dm.get_ss(s2, rmap, ss)
    pair = (min(s1, s2), max(s1, s2))
    i = rmap.get(s1, -1); j = rmap.get(s2, -1)
    zv = abs(z[i, j]) if i >= 0 and j >= 0 else 0.0
    dc = abs(dccm[i, j]) if i >= 0 and j >= 0 else 0.0
    bc1 = bc.get(s1, 0); bc2 = bc.get(s2, 0)
    n1 = set(G_nb.neighbors(s1)) if s1 in G_nb else set()
    n2 = set(G_nb.neighbors(s2)) if s2 in G_nb else set()
    common = n1 & n2
    adj_comm = (c1 != -1 and c2 != -1 and c1 != c2)
    both_coil = (ss1 == 'C' and ss2 == 'C')
    both_helix = (ss1 == 'H' and ss2 == 'H')
    parts = str(mut_str).split(',')
    ch1 = dm.is_charged_orig(parts[0]) and dm.is_charged_mut(parts[0]) if parts else False
    ch2 = dm.is_charged_orig(parts[1]) and dm.is_charged_mut(parts[1]) if len(parts) > 1 else False
    vol = (dm.has_major_vol_change(parts[0]) if parts else False) or \
          (dm.has_major_vol_change(parts[1]) if len(parts) > 1 else False)
    both_iso = (c1 == -1 and c2 == -1)
    one_iso = (c1 == -1 or c2 == -1)
    linked = (c1 != -1 and c2 != -1 and c1 != c2 and c2 in linked_comm.get(c1, set()))
    g = lambda name, cond: True if name in op else cond

    if same:
        hf = sum(1 for r in communities[c1] if dm.get_ss(r, rmap, ss) == 'H') / max(len(communities[c1]), 1)
        if hf >= dm.C1_HELIX_FRAC: return 'Additive', 'C1'
        if (ss1 == 'E' and ss2 == 'H') or (ss1 == 'H' and ss2 == 'E'): return 'Additive', 'C1'
        if g('C1', zv > dm.C1_COUPLING): return 'Non-additive', 'C1'
        return 'Additive', 'C1'
    if pair in direct and seq > 2 and g('C2', zv > dm.C2_COUPLING):
        if not (seq <= 3 and both_helix): return 'Non-additive', 'C2'
    if both_coil and seq > 2 and g('C3a', zv > dm.C3_COIL_COUPLING):
        return 'Non-additive', 'C3'
    if both_helix and one_iso and g('C3b', zv > dm.C3_HELIX_ONE_ISO_COUPLING) and seq <= dm.C3_HELIX_ISO_SEP_MAX:
        return 'Non-additive', 'C3'
    if (both_helix and both_iso and g('C3c', zv > dm.C3_HELIX_BOTH_ISO_COUPLING) and
            dm.C3_HELIX_BOTH_ISO_SEP_MIN < seq <= dm.C3_HELIX_BOTH_ISO_SEP_MAX and vol):
        return 'Non-additive', 'C3'
    if seq <= dm.C3_SHORT_HELIX_SEP and both_helix and adj_comm and g('C3d', zv > dm.C3_SHORT_HELIX_COUPLING):
        return 'Non-additive', 'C3'
    if seq <= 2 and g('C4', dc > dm.C4_BACKBONE_CORR) and both_iso:
        return 'Non-additive', 'C4'
    if len(common) > 0 and seq > 5 and not adj_comm and g('C5a', zv > dm.C5_COMMON_NEIGHBOR_COUPLING):
        return 'Non-additive', 'C5'
    if (bc1 > dm.C5_HUB_BETWEENNESS or bc2 > dm.C5_HUB_BETWEENNESS) and seq > dm.C5_HUB_SEP_MIN:
        if 'C5hub' not in op and zv < 0.15: return 'Additive', 'C5'
        return 'Non-additive', 'C5'
    if linked and g('C5b', zv > dm.C5_LINKED_COMM_COUPLING) and seq > dm.C5_LINKED_COMM_SEP_MIN:
        return 'Non-additive', 'C5'
    if ch1 and ch2 and both_iso and seq >= 8:
        return 'Non-additive', 'C6'
    return 'Additive', 'none'


VARIANTS = [
    ('None (full DynRCNC)',          set()),
    ('C1 clique community',          {'C1'}),
    ('C2 direct contact',            {'C2'}),
    ('C3 dynamic coupling',          {'C3a', 'C3b', 'C3c', 'C3d'}),
    ('C4 sequential backbone',       {'C4'}),
    ('C5 forward coupling cutoffs',  {'C5a', 'C5b'}),
    ('C5 hub exclusion gate',        {'C5hub'}),
    ('C6 electrostatic coupling',    set()),
    ('All dynamic conditions',       {'C1', 'C2', 'C3a', 'C3b', 'C3c', 'C3d', 'C4', 'C5a', 'C5b', 'C5hub'}),
]


def metrics(y, p):
    tp = int(((y == 1) & (p == 1)).sum()); tn = int(((y == 0) & (p == 0)).sum())
    fp = int(((y == 0) & (p == 1)).sum()); fn = int(((y == 1) & (p == 0)).sum())
    d = np.sqrt(float((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn)))
    return tp, tn, fp, fn, 100*tp/(tp+fn), 100*tn/(tn+fp), ((tp*tn-fp*fn)/d if d else 0.0)


def main():
    ddg = pd.read_csv(os.path.join(HERE, 'benchmark_271pairs.ddg'), sep='\t', on_bad_lines='skip')
    ddg.columns = ddg.columns.str.strip()
    rows = []   # one entry per pair: (y, {variant: pred})
    for pdb in KEEP:
        cfg = dm.PROTEINS[pdb]
        print(f'  {pdb} ...', flush=True)
        comm, r2c, direct, G_nb, bc, linked, last = dm.build_network(cfg['node_file'], cfg['edge_file'])
        resids, rmap, dccm, z, rmsf, ss = dm.compute_md_features(cfg['gro'], cfg['xtc'], cfg['fps'])
        for _, r in ddg[ddg['PDB'] == pdb].iterrows():
            parts = str(r['Mutation_Double']).split(',')
            nums = [int(''.join(c for c in p if c.isdigit())) for p in parts if any(c.isdigit() for c in p)]
            if len(nums) != 2: continue
            y = int(abs(float(r['dddG'])) >= dm.DDDG_CUT)
            pr = {name: int(predict_open(nums[0], nums[1], r['Mutation_Double'], comm, r2c, G_nb, bc,
                                         direct, linked, rmap, ss, dccm, z, op)[0] == 'Non-additive')
                  for name, op in VARIANTS}
            rows.append((y, pr))

    y = np.array([r[0] for r in rows])
    base = np.array([r[1][VARIANTS[0][0]] for r in rows])
    bm = metrics(y, base)[-1]
    out = []
    for name, op in VARIANTS:
        p = np.array([r[1][name] for r in rows])
        tp, tn, fp, fn, se, sp, m = metrics(y, p)
        out.append({'Dynamic condition(s) opened': name, 'n': len(op) if name != 'C6 electrostatic coupling' else '-',
                    'Flips': int((p != base).sum()), 'TP': tp, 'TN': tn, 'FP': fp, 'FN': fn,
                    'Sens': round(se, 1), 'Spec': round(sp, 1), 'MCC': round(m, 3), 'dMCC': round(m-bm, 3)})
    res = pd.DataFrame(out)
    print(f'\nPairs evaluated: {len(y)}')
    print(res.to_string(index=False))
    res.to_csv(os.path.join(HERE, 'ablation_dynamic_gates.csv'), index=False)
    print('\nWrote ablation_dynamic_gates.csv')


if __name__ == '__main__':
    main()
