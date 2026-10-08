#!/usr/bin/env python3
"""Figure 4 data (clean 265 pairs).
A: all 135 virtual-edge configurations (sens, spec, MCC) + DynRCNC, static RCNC, best VE markers
B: best MCC across distance cut-offs at each |Z-DCCM| threshold
C: per-pair |Z-DCCM| by class + Mann-Whitney P
Outputs: fig4_A_grid.csv, fig4_A_markers.csv, fig4_B.csv, fig4_C.csv
Run inside 7_figure_scripts/:  python3 fig4_data.py"""
import os, numpy as np, pandas as pd
from scipy.stats import mannwhitneyu
HERE=os.path.dirname(os.path.abspath(__file__))
R=os.path.expanduser(os.environ.get('RESULTS_DIR',os.path.join(os.path.dirname(os.path.abspath(__file__)),'..','5_results')))
g=pd.read_csv(os.path.join(R,'virtual_edge_grid_265.csv'))
g.to_csv(os.path.join(HERE,'fig4_A_grid.csv'),index=False)
b=g.loc[g.mcc.idxmax()]
mk=pd.DataFrame([dict(model='DynRCNC',sens=40/54,spec=184/211,mcc=0.568),
                 dict(model='Static RCNC',sens=4/54,spec=204/211,mcc=0.083),
                 dict(model='Best virtual edge',sens=b.sensitivity,spec=b.specificity,mcc=b.mcc)])
mk.to_csv(os.path.join(HERE,'fig4_A_markers.csv'),index=False)
B=g.groupby('zthr').mcc.max().reset_index().rename(columns={'mcc':'best_mcc'})
B.to_csv(os.path.join(HERE,'fig4_B.csv'),index=False)
d=pd.read_csv(os.path.join(R,'clean_all_predictions.csv'))
C=d[['PDB','mutation','exp','zdccm']]; C.to_csv(os.path.join(HERE,'fig4_C.csv'),index=False)
na=C[C.exp=='Non-additive'].zdccm; ad=C[C.exp=='Additive'].zdccm
print(f"A: {len(g)} configs, MCC range {g.mcc.min():.3f} to {g.mcc.max():.3f}; best Z>{b.zthr}, d<={b.dist:.0f} A: MCC {b.mcc:.3f}, sens {b.sensitivity:.3f}, spec {b.specificity:.3f}")
print("B:",', '.join(f"{z}:{m:.3f}" for z,m in zip(B.zthr,B.best_mcc)))
for n,s in [('Non-additive',na),('Additive',ad)]:
    print(f"C: {n:12s} n={len(s)} median={s.median():.3f} Q1={s.quantile(.25):.3f} Q3={s.quantile(.75):.3f} max={s.max():.3f}")
print(f"C: Mann-Whitney P = {mannwhitneyu(na,ad).pvalue:.3f}")
