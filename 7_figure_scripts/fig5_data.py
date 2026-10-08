#!/usr/bin/env python3
"""Figure 5 data (clean 265 pairs): in-sample and LOPOCV MCC per method.
DynRCNC: in-sample from clean_all_predictions.csv, LOPOCV = N-weighted mean of clean_lopocv.csv
Static RCNC: same rule as static_rcnc_265.py (no fitted thresholds, so LOPOCV = in-sample)
FoldX: untrained physics baseline (single value, shown as LOPOCV only, as in the submitted figure)
ML (default settings, as in the submitted figure): recomputed here at full precision from
  5_results/features_265.csv with the exact models of ml_baseline_265.py
Output: fig5_data.csv.  Run inside 7_figure_scripts/:  python3 fig5_data.py"""
import os, numpy as np, pandas as pd
H=os.path.dirname(os.path.abspath(__file__)); R=os.path.expanduser(os.environ.get('RESULTS_DIR',os.path.join(os.path.dirname(os.path.abspath(__file__)),'..','5_results')))
C=R
def mcc(y,p):
    tp=((y==1)&(p==1)).sum();tn=((y==0)&(p==0)).sum();fp=((y==0)&(p==1)).sum();fn=((y==1)&(p==0)).sum()
    return (tp*tn-fp*fn)/np.sqrt(float((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn)))
d=pd.read_csv(os.path.join(C,'clean_all_predictions.csv')); y=(d.exp=='Non-additive').astype(int).values
dyn_in=mcc(y,(d.pred=='Non-additive').astype(int).values)
l=pd.read_csv(os.path.join(C,'clean_lopocv.csv')); dyn_lo=(l.MCC*l.N).sum()/l.N.sum()
st=((d.comm_i!=-1)&(d.comm_i==d.comm_j)&~((d.ss_i=='H')&(d.ss_j=='H'))).astype(int).values; st_m=mcc(y,st)
f=pd.read_excel(os.path.join(os.path.dirname(os.path.abspath(__file__)),'..','6_supporting_data','FoldX_vs_Benchmark_vs_DynRCNC_265pairs.xlsx'),sheet_name='Per-pair'); f=f[f.PDB.isin(['1STN','1BNI','2LZM','1PGA','1CSP','2RN2','2CI2'])]
fx=mcc((f['Exp class']=='Non-additive').astype(int).values,(f['FoldX class']=='Non-additive').astype(int).values)
import warnings; warnings.filterwarnings('ignore')
from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
FEATS=['seq_dist','same_community','adj_comm','comm_linked','in_direct','z_dccm','raw_dccm','bc_max','bc_min','rmsf_z_max','rmsf_z_min',
       'n_common_neigh','both_coil','both_helix','both_struct','both_isolated','one_isolated','both_charged_kept','vol_change_any']
F=pd.read_csv(os.path.join(R,'features_265.csv')); X=F[FEATS].values.astype(float); yy=F.label.values.astype(int); G=F.PDB.values
def make(n):
    if n=='RandomForest': return RandomForestClassifier(n_estimators=300,class_weight='balanced',random_state=42,n_jobs=-1)
    if n=='GradientBoosting': return GradientBoostingClassifier(n_estimators=200,max_depth=3,random_state=42)
    return LogisticRegression(class_weight='balanced',max_iter=2000,random_state=42)
def mcc0(a,b):
    den=np.sqrt(float((((a==1)&(b==1)).sum()+((a==0)&(b==1)).sum())*(((a==1)&(b==1)).sum()+((a==1)&(b==0)).sum())*(((a==0)&(b==0)).sum()+((a==0)&(b==1)).sum())*(((a==0)&(b==0)).sum()+((a==1)&(b==0)).sum())))
    return mcc(a,b) if den>0 else 0.0
mlv={}
for n in ['RandomForest','GradientBoosting','LogisticRegression']:
    sc=StandardScaler().fit(X); m_in=mcc0(yy,make(n).fit(sc.transform(X),yy).predict(sc.transform(X)))
    w=0.0
    for h in np.unique(G):
        tr=G!=h; te=G==h; s2=StandardScaler().fit(X[tr])
        p=make(n).fit(s2.transform(X[tr]),yy[tr]).predict(s2.transform(X[te])); w+=mcc0(yy[te],p)*te.sum()
    mlv[n]=(m_in,w/len(yy))
rows=[('Static\nRCNC',st_m,st_m),('FoldX',np.nan,fx),
      ('Logistic\nreg.',*mlv['LogisticRegression']),
      ('Random\nforest',*mlv['RandomForest']),
      ('Gradient\nboost.',*mlv['GradientBoosting']),
      ('DynRCNC',dyn_in,dyn_lo)]
out=pd.DataFrame(rows,columns=['method','in_sample','lopocv'])
out.to_csv(os.path.join(H,'fig5_data.csv'),index=False)
print(out.assign(method=out.method.str.replace('\n',' ')).to_string(index=False,float_format=lambda v: f'{v:.3f}'))
