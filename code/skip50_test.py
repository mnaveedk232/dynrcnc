#!/usr/bin/env python3
"""Equilibration-window check: recompute DynRCNC on the clean 265 pairs after discarding the
first 50 ns of each 200 ns trajectory (DCCM/Z-DCCM from 50-200 ns only). Thresholds unchanged.
Run inside 04_baselines_265/:  python3 skip50_test.py"""
import os, sys, importlib.util, numpy as np, pandas as pd
CLEAN=os.path.expanduser(os.environ.get('CLEAN_DIR','~/Desktop/revision_analysis/00_clean_benchmark'))
sys.argv=['dynrcnc.py','--base',CLEAN,'--ddg',os.path.join(CLEAN,'benchmark_271pairs.ddg')]
spec=importlib.util.spec_from_file_location('dm',os.path.join(CLEAN,'dynrcnc.py'))
dm=importlib.util.module_from_spec(spec); spec.loader.exec_module(dm)
PS_PER_FRAME={'1STN':100,'1BNI':100,'2LZM':100,'1PGA':40,'1CSP':40,'2RN2':100,'2CI2':100}
ddg=pd.read_csv(os.path.join(CLEAN,'benchmark_271pairs.ddg'),sep='\t',on_bad_lines='skip'); ddg.columns=ddg.columns.str.strip()
def run(skip):
    dm.SKIP_NS=skip; yt=[];yp=[];per={}
    for pdb,fr in PS_PER_FRAME.items():
        c=dm.PROTEINS[pdb]
        comm,r2c,direct,G_nb,bc,linked,last=dm.build_network(c['node_file'],c['edge_file'])
        res,rmap,dccm,z,rmsf,ss=dm.compute_md_features(c['gro'],c['xtc'],fr)
        for _,r in ddg[ddg.PDB==pdb].iterrows():
            parts=str(r['Mutation_Double']).split(',')
            n=[int(''.join(ch for ch in p if ch.isdigit())) for p in parts if any(ch.isdigit() for ch in p)]
            if len(n)!=2: continue
            pr,_=dm.predict(n[0],n[1],r['Mutation_Double'],comm,r2c,G_nb,bc,direct,linked,last,rmap,ss,dccm,z,rmsf)
            yt.append(int(abs(float(r['dddG']))>=dm.DDDG_CUT)); yp.append(int(pr=='Non-additive'))
    return np.array(yt),np.array(yp)
def met(y,p):
    tp=((y==1)&(p==1)).sum();tn=((y==0)&(p==0)).sum();fp=((y==0)&(p==1)).sum();fn=((y==1)&(p==0)).sum()
    return tp,tn,fp,fn,(tp*tn-fp*fn)/np.sqrt(float((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn)))
y0,p0=run(0); y1,p1=run(50)
print("Full 0-200 ns : TP=%d TN=%d FP=%d FN=%d MCC=%.3f"%met(y0,p0))
print("Skip 50 ns    : TP=%d TN=%d FP=%d FN=%d MCC=%.3f"%met(y1,p1))
print("Pairs with identical prediction: %d/%d"%((p0==p1).sum(),len(p0)))
