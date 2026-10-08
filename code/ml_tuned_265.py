#!/usr/bin/env python3
"""Nested-tuned ML baselines on the CLEAN 265-pair benchmark (Table 6 'tuned' rows).
Outer loop: leave-one-protein-out. Inner loop: hyperparameter grid selected by
grouped (leave-one-protein-out) CV on the training proteins only, scored by MCC.
Reads features_265.csv (written by ml_baseline_265.py). Seed 42.
Run:  python3 ml_tuned_265.py      ->  ml_tuned_265_results.csv"""
import warnings; warnings.filterwarnings('ignore')
import os, numpy as np, pandas as pd
from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline
from sklearn.model_selection import GridSearchCV, LeaveOneGroupOut
from sklearn.metrics import make_scorer, matthews_corrcoef
SEED=42; HERE=os.path.dirname(os.path.abspath(__file__))
FEATS=['seq_dist','same_community','adj_comm','comm_linked','in_direct','z_dccm','raw_dccm','bc_max','bc_min',
       'rmsf_z_max','rmsf_z_min','n_common_neigh','both_coil','both_helix','both_struct','both_isolated',
       'one_isolated','both_charged_kept','vol_change_any']
df=pd.read_csv(os.path.join(HERE,'features_265.csv'))
X=df[FEATS].values.astype(float); y=df.label.values.astype(int); g=df.PDB.values
MODELS={
 'RandomForest (tuned)':(RandomForestClassifier(class_weight='balanced',random_state=SEED,n_jobs=-1),
    {'m__max_depth':[3,6,None],'m__min_samples_leaf':[1,5],'m__n_estimators':[100,300]}),
 'GradientBoosting (tuned)':(GradientBoostingClassifier(random_state=SEED),
    {'m__max_depth':[2,3],'m__min_samples_leaf':[1,5],'m__n_estimators':[100,200]}),
 'LogisticRegression (tuned)':(LogisticRegression(class_weight='balanced',max_iter=2000,random_state=SEED),
    {'m__C':[0.01,0.1,1,10]})}
def met(t,p):
    tp=((t==1)&(p==1)).sum();tn=((t==0)&(p==0)).sum();fp=((t==0)&(p==1)).sum();fn=((t==1)&(p==0)).sum()
    d=np.sqrt(float((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn)))
    return ((tp*tn-fp*fn)/d if d else 0.0), tp/(tp+fn), tn/(tn+fp)
def search(est,grid):
    return GridSearchCV(Pipeline([('s',StandardScaler()),('m',est)]),grid,
                        scoring=make_scorer(matthews_corrcoef),cv=LeaveOneGroupOut(),n_jobs=-1)
out=[]
for name,(est,grid) in MODELS.items():
    full=search(est,grid).fit(X,y,groups=g); m_in=met(y,full.predict(X))[0]
    yp=np.zeros_like(y); w=0.0
    for h in np.unique(g):
        tr=g!=h; te=g==h
        b=search(est,grid).fit(X[tr],y[tr],groups=g[tr]); yp[te]=b.predict(X[te])
        w+=met(y[te],yp[te])[0]*te.sum()
    w/=len(y); _,se,sp=met(y,yp)
    print(f'{name:28s} in-sample={m_in:.3f}  LOPOCV weighted={w:.3f}  gap={m_in-w:.3f}  Sens={se*100:.1f}  Spec={sp*100:.1f}',flush=True)
    out.append(dict(method=name,in_sample_MCC=round(m_in,4),LOPOCV_weighted_MCC=round(w,4),gap=round(m_in-w,4),
                    LOPOCV_sens=round(se,4),LOPOCV_spec=round(sp,4)))
pd.DataFrame(out).to_csv(os.path.join(HERE,'ml_tuned_265_results.csv'),index=False)
print('Saved ml_tuned_265_results.csv')
