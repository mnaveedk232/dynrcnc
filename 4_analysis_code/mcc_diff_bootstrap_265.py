#!/usr/bin/env python3
"""Protein-clustered bootstrap of the MCC difference DynRCNC - static RCNC (clean 265 pairs).
Run inside 04_baselines_265/:  python3 mcc_diff_bootstrap_265.py"""
import os, numpy as np, pandas as pd
d=pd.read_csv(os.path.expanduser('~/Desktop/revision_analysis/00_clean_benchmark/clean_all_predictions.csv'))
d['y']=(d.exp=='Non-additive').astype(int); d['p']=(d.pred=='Non-additive').astype(int)
d['s']=((d.comm_i!=-1)&(d.comm_i==d.comm_j)&~((d.ss_i=='H')&(d.ss_j=='H'))).astype(int)
def mcc(y,p):
    tp=((y==1)&(p==1)).sum();tn=((y==0)&(p==0)).sum();fp=((y==0)&(p==1)).sum();fn=((y==1)&(p==0)).sum()
    den=np.sqrt(float((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn))); return (tp*tn-fp*fn)/den if den else 0
rng=np.random.default_rng(42); pro=d.PDB.unique(); G={k:g for k,g in d.groupby('PDB')}
diff=[]
for _ in range(1000):
    s=pd.concat([G[k] for k in rng.choice(pro,len(pro))])
    diff.append(mcc(s.y.values,s.p.values)-mcc(s.y.values,s.s.values))
    rng.integers(1e9)
diff=np.array(diff); obs=mcc(d.y.values,d.p.values)-mcc(d.y.values,d.s.values)
print(f"Observed dMCC={obs:.3f}  protein-clustered 95% CI {np.percentile(diff,2.5):.3f} to {np.percentile(diff,97.5):.3f}  fraction<=0: {(diff<=0).mean():.3f}")
