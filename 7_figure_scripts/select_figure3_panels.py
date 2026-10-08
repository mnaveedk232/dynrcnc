#!/usr/bin/env python3
"""
Figure 3 panel selection: one representative pair per criterion (C1 to C6).
For each criterion it lists every correctly predicted non-additive pair decided by that criterion,
and marks the suggested panel: the one with the largest |dddG|, i.e. the clearest experimental
non-additivity among the correct calls. Everything needed for the PyMOL panels is printed:
PDB, the two residue numbers, the mutation string, |dddG|, |Z-DCCM|, sequence separation,
community membership and secondary structure of both sites.

Input : 5_results/clean_all_predictions.csv
Output: 7_figure_scripts/figure3_panels.csv  (suggested panels only)
        7_figure_scripts/figure3_candidates.csv  (all correct calls per criterion)
Run [local]:
  cd 7_figure_scripts
  python3 select_figure3_panels.py
"""
import os, pandas as pd
H = os.path.dirname(os.path.abspath(__file__))
CSV = os.path.join(H, '..', '5_results', 'clean_all_predictions.csv')
NAME = {'C1':'C1  clique community','C2':'C2  direct contact','C3':'C3  dynamic packing',
        'C4':'C4  sequential backbone','C5':'C5  network topology','C6':'C6  electrostatic coupling'}
d = pd.read_csv(CSV)
d['criterion'] = d.mechanism.str.split('_').str[1]
ok = d[(d.exp == 'Non-additive') & (d.pred == 'Non-additive')].copy()
ok['absd'] = ok.dddG.abs()
cols = ['criterion','PDB','mutation','res_i','res_j','dddG','zdccm','seq_dist','comm_i','comm_j','ss_i','ss_j','mechanism','pH','Method']
cand, panels = [], []
for c in ['C1','C2','C3','C4','C5','C6']:
    g = ok[ok.criterion == c].sort_values('absd', ascending=False)
    print(f"\n===== {NAME[c]}: {len(g)} correctly predicted non-additive pairs")
    if len(g) == 0:
        print("   none"); continue
    print(g[cols].to_string(index=False))
    best = g.iloc[0]
    print(f"   -> suggested panel: {best.PDB}  {best.mutation}  "
          f"(residues {int(best.res_i)} and {int(best.res_j)}, |dddG| = {abs(best.dddG):.2f} kcal/mol, "
          f"|Z-DCCM| = {best.zdccm:.2f}, rule {best.mechanism})")
    cand.append(g[cols]); panels.append(best[cols])
pd.concat(cand).to_csv(os.path.join(H,'figure3_candidates.csv'), index=False)
p = pd.DataFrame(panels)[cols]
p.insert(0,'panel',list('ABCDEF')[:len(p)])
p.to_csv(os.path.join(H,'figure3_panels.csv'), index=False)
print("\n===== SUGGESTED FIGURE 3 PANELS")
print(p.to_string(index=False))
print(f"\nSaved {os.path.join(H,'figure3_panels.csv')}")
print(f"Saved {os.path.join(H,'figure3_candidates.csv')}")
