#!/usr/bin/env python3
"""Verify and export the RING contact files and residue-level DCCM / Z-DCCM matrices.

For each of the 7 benchmark proteins this script
  1. checks that the RING node/edge files belong to that protein (residues and residue
     names equal those of the MD topology) and cover every benchmark mutation site,
  2. recomputes DCCM and Z-DCCM from the wild-type trajectory with the same function
     that dynrcnc.py uses (compute_md_features), so the matrices are the ones behind
     the manuscript,
  3. compares |Z-DCCM| at every benchmark pair with the zdccm column of
     results/clean_all_predictions.csv (the 0.568 run),
  4. copies the files to the upload folder ONLY for proteins that pass every check.

  python3 export_ring_dccm.py --base BASE --dynrcnc ~/Desktop/github2/dynrcnc.py \
      --pred ~/Desktop/github2/results/clean_all_predictions.csv --out ~/Desktop/github_upload
BASE holds one <pdb>_md/ folder per protein (.N, .E, *_first_frame.gro, *_reduced.xtc).
"""
import argparse, os, sys, shutil, importlib.util
import numpy as np, pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument('--base', required=True)
ap.add_argument('--dynrcnc', default=os.path.expanduser('~/Desktop/github2/dynrcnc.py'))
ap.add_argument('--pred', default=os.path.expanduser('~/Desktop/github2/results/clean_all_predictions.csv'))
ap.add_argument('--xtcdir', default=None, help='folder holding <pdb>_md_reduced.xtc (overrides the xtc path under --base)')
ap.add_argument('--out', default=os.path.expanduser('~/Desktop/github_upload'))
a = ap.parse_args()
base, out = [os.path.expanduser(p) for p in (a.base, a.out)]

sys.argv = [sys.argv[0], '--base', base]            # dynrcnc.py reads --base at import
spec = importlib.util.spec_from_file_location('dynrcnc', os.path.expanduser(a.dynrcnc))
dyn = importlib.util.module_from_spec(spec); spec.loader.exec_module(dyn)

if a.xtcdir:
    for _p, _c in dyn.PROTEINS.items():
        _c['xtc'] = os.path.join(os.path.expanduser(a.xtcdir), _p.lower() + '_md_reduced.xtc')

pred = pd.read_csv(os.path.expanduser(a.pred))
ring_dir = os.path.join(out, '2_RING_contact_files')
dccm_dir = os.path.join(out, '3_residue_level_DCCM_matrices')
os.makedirs(ring_dir, exist_ok=True); os.makedirs(dccm_dir, exist_ok=True)

def gro_resnames(gro):
    names = {}
    for l in open(gro).read().splitlines()[2:]:
        if len(l) < 20: continue
        try: names[int(l[0:5])] = l[5:10].strip()
        except ValueError: pass
    return names

report = []; all_ok = True
for pdb, cfg in dyn.PROTEINS.items():
    prob = []
    for k in ('node_file', 'edge_file', 'gro', 'xtc'):
        if not os.path.exists(cfg[k]): prob.append('missing ' + cfg[k])
    if prob:
        report.append((pdb, 'FAIL', '; '.join(prob))); all_ok = False; continue
    # 1. RING belongs to this protein
    n = pd.read_table(cfg['node_file'], encoding='latin-1'); n.columns = n.columns.str.strip()
    ring = {int(r.split(':')[1]): r.split(':')[3] for r in n['NodeId']}
    gro = gro_resnames(cfg['gro'])
    only_ring = sorted(set(ring) - set(gro)); only_gro = sorted(set(gro) - set(ring))
    name_bad = [r for r in ring if r in gro and ring[r].upper()[:3] != gro[r].upper()[:3]
                and not (ring[r].upper().startswith('HI') and gro[r].upper().startswith('HI'))]
    sub = pred[pred.PDB == pdb]
    sites = set(sub.res_i) | set(sub.res_j)
    miss_sites = sorted(sites - set(ring))
    if only_ring or only_gro or name_bad or miss_sites:
        prob.append('RING vs topology: only in RING %s, only in topology %s, residue-name mismatches %s, mutation sites not in RING %s'
                    % (only_ring[:5], only_gro[:5], name_bad[:5], miss_sites[:5]))
    # 2. matrices from the trajectory
    resids, rmap, dccm, z, rmsf, ss = dyn.compute_md_features(cfg['gro'], cfg['xtc'], cfg['fps'])
    # 3. match with the manuscript run
    bad = 0
    for _, r in sub.iterrows():
        i, j = rmap.get(int(r.res_i), -1), rmap.get(int(r.res_j), -1)
        v = round(abs(z[i, j]), 4) if i >= 0 and j >= 0 else 0
        if abs(v - float(r.zdccm)) > 1e-3: bad += 1
    if bad: prob.append('%d of %d pairs differ from clean_all_predictions.csv (wrong trajectory?)' % (bad, len(sub)))
    status = 'FAIL' if prob else 'PASS'
    report.append((pdb, status, '; '.join(prob) if prob else
                   '%d residues, %d pairs: |Z-DCCM| identical to the 0.568 run' % (len(resids), len(sub))))
    if prob: all_ok = False; continue
    # 4. export only verified files
    d = os.path.join(ring_dir, pdb); os.makedirs(d, exist_ok=True)
    shutil.copy(cfg['node_file'], os.path.join(d, pdb.lower() + '.N'))
    shutil.copy(cfg['edge_file'], os.path.join(d, pdb.lower() + '.E'))
    lab = ['%s%d' % (gro[r], r) for r in resids]
    pd.DataFrame(dccm, index=lab, columns=lab).round(5).to_csv(os.path.join(dccm_dir, pdb + '_DCCM.csv'))
    pd.DataFrame(z, index=lab, columns=lab).round(5).to_csv(os.path.join(dccm_dir, pdb + '_Z-DCCM.csv'))

print('\n' + '=' * 70)
for pdb, st, msg in report: print('%-5s %-5s %s' % (st, pdb, msg))
print('=' * 70)
print('ALL 7 PROTEINS PASS. Output folder: ' + out if all_ok else
      'AT LEAST ONE PROTEIN FAILED. Upload only the folders listed as PASS. Output folder: ' + out)
for root, _, fs in os.walk(out):
    for f in fs:
        p = os.path.join(root, f)
        if os.path.getsize(p) > 20e6: print('WARNING above 20 MB:', p)
