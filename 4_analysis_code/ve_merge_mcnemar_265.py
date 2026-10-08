#!/usr/bin/env python3
"""Merge ve_out/ve2_*.csv into the 135-cell virtual-edge grid
and run McNemar DynRCNC vs best virtual edge (logic of mcnemar_ve.py), CLEAN 265 pairs.
Run inside 04_baselines_265/:  python3 ve_merge_mcnemar_265.py"""
import os, glob, math, numpy as np, pandas as pd
from scipy.stats import chi2
HERE=os.path.dirname(os.path.abspath(__file__))
CLEAN=os.path.expanduser('~/Desktop/revision_analysis/00_clean_benchmark')
big=pd.concat([pd.read_csv(f,keep_default_na=False) for f in sorted(glob.glob(f'{HERE}/ve_out/ve2_*.csv'))],ignore_index=True)
def mcc(tp,tn,fp,fn):
    d=math.sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn)); return (tp*tn-fp*fn)/d if d>0 else 0.0
rows=[]
for (zt,dc),g in big.groupby(['zthr','dist']):
    TP=((g.pred=='NA')&(g.exp=='Non-additive')).sum(); TN=((g.pred=='ADD')&(g.exp=='Additive')).sum()
    FP=((g.pred=='NA')&(g.exp=='Additive')).sum(); FN=((g.pred=='ADD')&(g.exp=='Non-additive')).sum()
    rows.append(dict(zthr=zt,dist=dc,TP=int(TP),TN=int(TN),FP=int(FP),FN=int(FN),
                     sensitivity=TP/(TP+FN),specificity=TN/(TN+FP),mcc=mcc(TP,TN,FP,FN)))
grid=pd.DataFrame(rows); grid.to_csv(f'{HERE}/virtual_edge_grid_265.csv',index=False)
b=grid.loc[grid.mcc.idxmax()]
print(f"Pairs in grid: {len(big)//135}")
print(f"BEST: Z>{b.zthr}, d<={b.dist} A -> MCC={b.mcc:.3f} Sens={b.sensitivity*100:.1f}% Spec={b.specificity*100:.1f}% (TP={int(b.TP)} TN={int(b.TN)} FP={int(b.FP)} FN={int(b.FN)})")
print(f"MCC range: {grid.mcc.min():.3f} to {grid.mcc.max():.3f}")
p=pd.read_csv(f'{CLEAN}/clean_all_predictions.csv')
ve={}
for f in glob.glob(f'{HERE}/ve_out/ve2_*.csv'):
    d=pd.read_csv(f,keep_default_na=False); pdb=os.path.basename(f)[4:8]
    for _,r in d[(d.zthr==b.zthr)&(d.dist==b.dist)].iterrows():
        ve[(pdb,)+tuple(sorted((int(r.s1),int(r.s2))))]=int(r.pred=='NA')
yt=(p.exp=='Non-additive').astype(int).values; yp=(p.pred=='Non-additive').astype(int).values
yv=np.array([ve[(r.PDB,)+tuple(sorted((int(r.res_i),int(r.res_j))))] for r in p.itertuples()])
bb=int(((yp==yt)&(yv!=yt)).sum()); cc=int(((yp!=yt)&(yv==yt)).sum())
print(f"McNemar DynRCNC vs best VE ({len(yt)} rows): b={bb} c={cc} P={chi2.sf((abs(bb-cc)-1)**2/(bb+cc),1):.2e}")
