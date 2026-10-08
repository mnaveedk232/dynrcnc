#!/usr/bin/env python3
"""Static RCNC baseline (Ming 2018 / Zhang 2024 rule, as in rcnc_dyn.py predict_rcnc),
Mann-Whitney on per-pair |Z-DCCM|, and FoldX, all on the CLEAN 265 pairs.
Run inside 04_baselines_265/:  python3 stats_265.py"""
import os, numpy as np, pandas as pd
from scipy.stats import mannwhitneyu
CLEAN=os.path.expanduser('~/Desktop/revision_analysis/00_clean_benchmark')
FOLDX=os.path.join(os.path.dirname(os.path.abspath(__file__)),'..','6_supporting_data','FoldX_vs_Benchmark_vs_DynRCNC_265pairs.xlsx')
def show(name,y,q):
    tp=(y&q).sum();tn=(~y&~q).sum();fp=(~y&q).sum();fn=(y&~q).sum()
    m=(tp*tn-fp*fn)/np.sqrt(float((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn)))
    print(f"{name:12s} N={len(y)} TP={tp} TN={tn} FP={fp} FN={fn} Sens={100*tp/(tp+fn):.1f}% Spec={100*tn/(tn+fp):.1f}% MCC={m:.3f}")
d=pd.read_csv(f'{CLEAN}/clean_all_predictions.csv'); y=(d.exp=='Non-additive').values
same=((d.comm_i!=-1)&(d.comm_i==d.comm_j)).values; helix=((d.ss_i=='H')&(d.ss_j=='H')).values
show('Static RCNC',y,same&~helix)
print(f"Mann-Whitney |Z-DCCM| NA vs ADD: P={mannwhitneyu(d.zdccm[y],d.zdccm[~y]).pvalue:.3f}  (median NA={d.zdccm[y].median():.3f}, ADD={d.zdccm[~y].median():.3f})")
f=pd.read_excel(FOLDX,sheet_name='Per-pair').rename(columns={'Exp class':'exp','FoldX class':'pred'}); f=f[f.PDB.isin(['1STN','1BNI','2LZM','1PGA','1CSP','2RN2','2CI2'])]
show('FoldX',(f.exp=='Non-additive').values,(f.pred=='Non-additive').values)
