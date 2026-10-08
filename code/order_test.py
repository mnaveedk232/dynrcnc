#!/usr/bin/env python3
"""Cascade-order test (clean 265 pairs): evaluates every criterion independently for each pair,
counts pairs on which more than one criterion fires, and compares the confusion matrix of the
published order (C1..C6) with the reverse tier order (C6, C5, C4, C3, C2, C1).
Run inside 04_baselines_265/:  python3 order_test.py"""
import os, sys, importlib.util, numpy as np, pandas as pd
C=os.path.expanduser(os.environ.get('CLEAN_DIR','~/Desktop/revision_analysis/00_clean_benchmark'))
sys.argv=['dynrcnc.py','--base',C,'--ddg',os.path.join(C,'benchmark_271pairs.ddg')]
spec=importlib.util.spec_from_file_location('dm',os.path.join(C,'dynrcnc.py')); dm=importlib.util.module_from_spec(spec); spec.loader.exec_module(dm)
def rules(s1,s2,mut,comm,r2c,G_nb,bc,direct,linked,rmap,ss,dccm,z):
    """returns list of (criterion, verdict) for every criterion whose condition holds, in published order"""
    seq=abs(s1-s2); c1=r2c.get(s1,-1); c2=r2c.get(s2,-1); same=(c1!=-1 and c2!=-1 and c1==c2)
    ss1=dm.get_ss(s1,rmap,ss); ss2=dm.get_ss(s2,rmap,ss); pair=(min(s1,s2),max(s1,s2))
    i=rmap.get(s1,-1); j=rmap.get(s2,-1); zv=abs(z[i,j]) if i>=0 and j>=0 else 0.0; dc=abs(dccm[i,j]) if i>=0 and j>=0 else 0.0
    bc1=bc.get(s1,0); bc2=bc.get(s2,0)
    n1=set(G_nb.neighbors(s1)) if s1 in G_nb else set(); n2=set(G_nb.neighbors(s2)) if s2 in G_nb else set()
    adj=(c1!=-1 and c2!=-1 and c1!=c2); coil=(ss1=='C' and ss2=='C'); hel=(ss1=='H' and ss2=='H')
    p=str(mut).split(',')
    ch1=dm.is_charged_orig(p[0]) and dm.is_charged_mut(p[0]); ch2=dm.is_charged_orig(p[1]) and dm.is_charged_mut(p[1]) if len(p)>1 else False
    vol=dm.has_major_vol_change(p[0]) or (dm.has_major_vol_change(p[1]) if len(p)>1 else False)
    biso=(c1==-1 and c2==-1); oiso=(c1==-1 or c2==-1)
    lk=(c1!=-1 and c2!=-1 and c1!=c2 and c2 in linked.get(c1,set()))
    out=[]
    if same:
        hf=sum(1 for r in comm[c1] if dm.get_ss(r,rmap,ss)=='H')/max(len(comm[c1]),1)
        if hf>=dm.C1_HELIX_FRAC or (ss1,ss2) in [('E','H'),('H','E')]: out.append(('C1','A'))
        else: out.append(('C1','N' if zv>dm.C1_COUPLING else 'A'))
    if pair in direct and seq>2 and zv>dm.C2_COUPLING and not (seq<=3 and hel): out.append(('C2','N'))
    if (coil and seq>2 and zv>dm.C3_COIL_COUPLING) or (hel and oiso and zv>dm.C3_HELIX_ONE_ISO_COUPLING and seq<=dm.C3_HELIX_ISO_SEP_MAX) or \
       (hel and biso and zv>dm.C3_HELIX_BOTH_ISO_COUPLING and dm.C3_HELIX_BOTH_ISO_SEP_MIN<seq<=dm.C3_HELIX_BOTH_ISO_SEP_MAX and vol) or \
       (seq<=dm.C3_SHORT_HELIX_SEP and hel and adj and zv>dm.C3_SHORT_HELIX_COUPLING): out.append(('C3','N'))
    if seq<=2 and dc>dm.C4_BACKBONE_CORR and biso: out.append(('C4','N'))
    c5=None
    if len(n1&n2)>0 and seq>5 and not adj and zv>dm.C5_COMMON_NEIGHBOR_COUPLING: c5='N'
    elif (bc1>dm.C5_HUB_BETWEENNESS or bc2>dm.C5_HUB_BETWEENNESS) and seq>dm.C5_HUB_SEP_MIN: c5='A' if zv<0.15 else 'N'
    elif lk and zv>dm.C5_LINKED_COMM_COUPLING and seq>dm.C5_LINKED_COMM_SEP_MIN: c5='N'
    if c5: out.append(('C5',c5))
    if ch1 and ch2 and biso and seq>=8: out.append(('C6','N'))
    return out
ddg=pd.read_csv(os.path.join(C,'benchmark_271pairs.ddg'),sep='\t',on_bad_lines='skip'); ddg.columns=ddg.columns.str.strip()
Y=[];F=[];R=[];multi=0;conflict=0
for pdb in ['1STN','1BNI','2LZM','1PGA','1CSP','2RN2','2CI2']:
    c=dm.PROTEINS[pdb]; comm,r2c,direct,G_nb,bc,linked,last=dm.build_network(c['node_file'],c['edge_file'])
    res,rmap,dccm,z,rmsf,ss=dm.compute_md_features(c['gro'],c['xtc'],c['fps'])
    for _,r in ddg[ddg.PDB==pdb].iterrows():
        pp=str(r['Mutation_Double']).split(',')
        n=[int(''.join(ch for ch in q if ch.isdigit())) for q in pp if any(ch.isdigit() for ch in q)]
        if len(n)!=2: continue
        o=rules(n[0],n[1],r['Mutation_Double'],comm,r2c,G_nb,bc,direct,linked,rmap,ss,dccm,z)
        if len(o)>1: multi+=1; conflict+=len(set(v for _,v in o))>1
        F.append(o[0][1] if o else 'A'); R.append(o[-1][1] if o else 'A'); Y.append(int(abs(float(r['dddG']))>=dm.DDDG_CUT))
Y=np.array(Y)
def cm(P):
    p=np.array([x=='N' for x in P]).astype(int); tp=((Y==1)&(p==1)).sum();tn=((Y==0)&(p==0)).sum();fp=((Y==0)&(p==1)).sum();fn=((Y==1)&(p==0)).sum()
    return tp,tn,fp,fn,(tp*tn-fp*fn)/np.sqrt(float((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn)))
print(f"Pairs: {len(Y)}; pairs with >1 criterion firing: {multi}; of these with conflicting verdicts: {conflict}")
print("Published order  C1..C6 : TP=%d TN=%d FP=%d FN=%d MCC=%.3f"%cm(F))
print("Reverse tier order C6..C1: TP=%d TN=%d FP=%d FN=%d MCC=%.3f"%cm(R))
print("Pairs whose predicted class differs:",sum(a!=b for a,b in zip(F,R)))
