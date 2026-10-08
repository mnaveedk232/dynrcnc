#!/usr/bin/env python3
"""Static RCNC baseline run directly through the model pipeline (clean 265 pairs).
Rule = rcnc_dyn.py predict_rcnc (Ming 2018 / Zhang 2024): same k=3 clique community
-> Non-additive, except both sites in helix -> Additive; otherwise Additive.
Run inside 04_baselines_265/:  python3 static_rcnc_265.py"""
import os, sys, importlib.util, numpy as np, pandas as pd
CLEAN=os.path.expanduser(os.environ.get('CLEAN_DIR','~/Desktop/revision_analysis/00_clean_benchmark'))
sys.argv=['dynrcnc.py','--base',CLEAN,'--ddg',os.path.join(CLEAN,'benchmark_265pairs.ddg')]
spec=importlib.util.spec_from_file_location('dm',os.path.join(CLEAN,'dynrcnc.py'))
dm=importlib.util.module_from_spec(spec); spec.loader.exec_module(dm)
ddg=pd.read_csv(os.path.join(CLEAN,'benchmark_265pairs.ddg'),sep='\t',on_bad_lines='skip'); ddg.columns=ddg.columns.str.strip()
yt=[];yp=[]
for pdb in ['1STN','1BNI','2LZM','1PGA','1CSP','2RN2','2CI2']:
    c=dm.PROTEINS[pdb]
    comm,r2c,direct,G_nb,bc,linked,last=dm.build_network(c['node_file'],c['edge_file'])
    res,rmap,dccm,z,rmsf,ss=dm.compute_md_features(c['gro'],c['xtc'],c['fps'])
    for _,r in ddg[ddg.PDB==pdb].iterrows():
        parts=str(r['Mutation_Double']).split(',')
        n=[int(''.join(ch for ch in p if ch.isdigit())) for p in parts if any(ch.isdigit() for ch in p)]
        if len(n)!=2: continue
        c1=r2c.get(n[0],-1); c2=r2c.get(n[1],-1)
        same=(c1!=-1 and c1==c2)
        helix=(dm.get_ss(n[0],rmap,ss)=='H' and dm.get_ss(n[1],rmap,ss)=='H')
        yp.append(int(same and not helix)); yt.append(int(abs(float(r['dddG']))>=dm.DDDG_CUT))
y=np.array(yt,bool); q=np.array(yp,bool)
tp=(y&q).sum();tn=(~y&~q).sum();fp=(~y&q).sum();fn=(y&~q).sum()
m=(tp*tn-fp*fn)/np.sqrt(float((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn)))
print(f"Static RCNC (model pipeline) N={len(y)} TP={tp} TN={tn} FP={fp} FN={fn} Sens={100*tp/(tp+fn):.1f}% Spec={100*tn/(tn+fp):.1f}% MCC={m:.3f}")
