#!/usr/bin/env python3
"""Figure 2 data (clean 265-pair benchmark): per-protein MCC for DynRCNC and static RCNC,
pair-level bootstrap 95% CI (1000 resamples within each protein, seed 42), TP/NA labels,
and the improvement DynRCNC - static. Output: fig2_data.csv
Run inside 7_figure_scripts/:  python3 fig2_data.py"""
import os, numpy as np, pandas as pd
CLEAN=os.path.expanduser(os.environ.get('RESULTS_DIR',os.path.join(os.path.dirname(os.path.abspath(__file__)),'..','5_results')))
HERE=os.path.dirname(os.path.abspath(__file__))
d=pd.read_csv(os.path.join(CLEAN,'clean_all_predictions.csv'))
d['y']=(d.exp=='Non-additive').astype(int); d['p']=(d.pred=='Non-additive').astype(int)
d['s']=((d.comm_i!=-1)&(d.comm_i==d.comm_j)&~((d.ss_i=='H')&(d.ss_j=='H'))).astype(int)
def mcc(y,p):
    tp=((y==1)&(p==1)).sum();tn=((y==0)&(p==0)).sum();fp=((y==0)&(p==1)).sum();fn=((y==1)&(p==0)).sum()
    den=np.sqrt(float((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn))); return (tp*tn-fp*fn)/den if den else 0.0
rng=np.random.default_rng(42); rows=[]
for pdb,g in d.groupby('PDB'):
    y=g.y.values; p=g.p.values; s=g.s.values; n=len(g); bd=[];bs=[]
    for _ in range(1000):
        i=rng.integers(0,n,n); bd.append(mcc(y[i],p[i])); bs.append(mcc(y[i],s[i]))
    rows.append(dict(PDB=pdb,N=n,NA=int(y.sum()),TP_dyn=int(((y==1)&(p==1)).sum()),TP_static=int(((y==1)&(s==1)).sum()),
        MCC_dyn=round(mcc(y,p),3),CI_dyn_lo=round(np.percentile(bd,2.5),3),CI_dyn_hi=round(np.percentile(bd,97.5),3),
        MCC_static=round(mcc(y,s),3),CI_static_lo=round(np.percentile(bs,2.5),3),CI_static_hi=round(np.percentile(bs,97.5),3)))
r=pd.DataFrame(rows); r['Improvement']=(r.MCC_dyn-r.MCC_static).round(3)
r=r.sort_values('MCC_dyn',ascending=False)
r.to_csv(os.path.join(HERE,'fig2_data.csv'),index=False)
pd.set_option('display.width',200); print(r.to_string(index=False))
