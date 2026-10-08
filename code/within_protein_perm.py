#!/usr/bin/env python3
"""Within-protein label permutation P and protein-clustered bootstrap CI of MCC (clean 265).
Run inside 04_baselines_265/:  python3 within_protein_perm.py"""
import os, numpy as np, pandas as pd
d=pd.read_csv(os.path.expanduser('~/Desktop/revision_analysis/00_clean_benchmark/clean_all_predictions.csv'))
y=(d.exp=='Non-additive').astype(int).values; p=(d.pred=='Non-additive').astype(int).values; g=d.PDB.values
def mcc(y,p):
    tp=((y==1)&(p==1)).sum();tn=((y==0)&(p==0)).sum();fp=((y==0)&(p==1)).sum();fn=((y==1)&(p==0)).sum()
    den=np.sqrt(float((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn))); return (tp*tn-fp*fn)/den if den else 0
obs=mcc(y,p); rng=np.random.default_rng(42); c=0
for _ in range(1000):
    yy=y.copy()
    for k in np.unique(g):
        i=np.where(g==k)[0]; yy[i]=rng.permutation(yy[i])
    c+=mcc(yy,p)>=obs
G={k:np.where(g==k)[0] for k in np.unique(g)}; ks=list(G); b=[]
for _ in range(1000):
    idx=np.concatenate([G[k] for k in rng.choice(ks,len(ks))]); b.append(mcc(y[idx],p[idx]))
print(f"MCC={obs:.3f}  within-protein permutation P={(c+1)/1001:.4f}  protein-clustered 95% CI {np.percentile(b,2.5):.3f} to {np.percentile(b,97.5):.3f}")
