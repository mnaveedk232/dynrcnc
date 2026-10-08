#!/usr/bin/env python3
"""Figure S2 data: one-at-a-time sensitivity of all 18 numerical thresholds (Table S1) on the
clean 265-pair benchmark. Each threshold is swept over a range while the other 17 stay at their
adopted values. MD features are computed once per protein with dynrcnc.py's own functions.
Outputs: figS2_sweep.csv (threshold, value, MCC) and figS2_summary.csv
Run inside 7_figure_scripts/:  python3 figS2_data.py"""
import os, sys, importlib.util, numpy as np, pandas as pd
H=os.path.dirname(os.path.abspath(__file__))
C=os.path.expanduser(os.environ.get('MD_BASE',os.path.join(os.path.dirname(os.path.abspath(__file__)),'..','data')))
DDG=os.path.join(H,'..','1_benchmark_dataset','benchmark_265pairs.ddg')
sys.argv=['dynrcnc.py','--base',C,'--ddg',DDG]
spec=importlib.util.spec_from_file_location('dynrcnc',os.path.join(H,'..','4_analysis_code','dynrcnc.py'))
dm=importlib.util.module_from_spec(spec); spec.loader.exec_module(dm)
HUB_MIN=[0.15]
def predict_open(s1, s2, mut_str, communities, r2c, G_nb, bc, direct,
                 linked_comm, rmap, ss, dccm, z, op):
    """Copy of dynrcnc.predict() with switchable dynamic conditions.
    op = set of opened gates. Empty set reproduces predict() exactly."""
    seq = abs(s1-s2)
    c1 = r2c.get(s1, -1); c2 = r2c.get(s2, -1)
    same = (c1 != -1 and c2 != -1 and c1 == c2)
    ss1 = dm.get_ss(s1, rmap, ss); ss2 = dm.get_ss(s2, rmap, ss)
    pair = (min(s1, s2), max(s1, s2))
    i = rmap.get(s1, -1); j = rmap.get(s2, -1)
    zv = abs(z[i, j]) if i >= 0 and j >= 0 else 0.0
    dc = abs(dccm[i, j]) if i >= 0 and j >= 0 else 0.0
    bc1 = bc.get(s1, 0); bc2 = bc.get(s2, 0)
    n1 = set(G_nb.neighbors(s1)) if s1 in G_nb else set()
    n2 = set(G_nb.neighbors(s2)) if s2 in G_nb else set()
    common = n1 & n2
    adj_comm = (c1 != -1 and c2 != -1 and c1 != c2)
    both_coil = (ss1 == 'C' and ss2 == 'C')
    both_helix = (ss1 == 'H' and ss2 == 'H')
    parts = str(mut_str).split(',')
    ch1 = dm.is_charged_orig(parts[0]) and dm.is_charged_mut(parts[0]) if parts else False
    ch2 = dm.is_charged_orig(parts[1]) and dm.is_charged_mut(parts[1]) if len(parts) > 1 else False
    vol = (dm.has_major_vol_change(parts[0]) if parts else False) or \
          (dm.has_major_vol_change(parts[1]) if len(parts) > 1 else False)
    both_iso = (c1 == -1 and c2 == -1)
    one_iso = (c1 == -1 or c2 == -1)
    linked = (c1 != -1 and c2 != -1 and c1 != c2 and c2 in linked_comm.get(c1, set()))
    g = lambda name, cond: True if name in op else cond

    if same:
        hf = sum(1 for r in communities[c1] if dm.get_ss(r, rmap, ss) == 'H') / max(len(communities[c1]), 1)
        if hf >= dm.C1_HELIX_FRAC: return 'Additive', 'C1'
        if (ss1 == 'E' and ss2 == 'H') or (ss1 == 'H' and ss2 == 'E'): return 'Additive', 'C1'
        if g('C1', zv > dm.C1_COUPLING): return 'Non-additive', 'C1'
        return 'Additive', 'C1'
    if pair in direct and seq > 2 and g('C2', zv > dm.C2_COUPLING):
        if not (seq <= 3 and both_helix): return 'Non-additive', 'C2'
    if both_coil and seq > 2 and g('C3a', zv > dm.C3_COIL_COUPLING):
        return 'Non-additive', 'C3'
    if both_helix and one_iso and g('C3b', zv > dm.C3_HELIX_ONE_ISO_COUPLING) and seq <= dm.C3_HELIX_ISO_SEP_MAX:
        return 'Non-additive', 'C3'
    if (both_helix and both_iso and g('C3c', zv > dm.C3_HELIX_BOTH_ISO_COUPLING) and
            dm.C3_HELIX_BOTH_ISO_SEP_MIN < seq <= dm.C3_HELIX_BOTH_ISO_SEP_MAX and vol):
        return 'Non-additive', 'C3'
    if seq <= dm.C3_SHORT_HELIX_SEP and both_helix and adj_comm and g('C3d', zv > dm.C3_SHORT_HELIX_COUPLING):
        return 'Non-additive', 'C3'
    if seq <= 2 and g('C4', dc > dm.C4_BACKBONE_CORR) and both_iso:
        return 'Non-additive', 'C4'
    if len(common) > 0 and seq > 5 and not adj_comm and g('C5a', zv > dm.C5_COMMON_NEIGHBOR_COUPLING):
        return 'Non-additive', 'C5'
    if (bc1 > dm.C5_HUB_BETWEENNESS or bc2 > dm.C5_HUB_BETWEENNESS) and seq > dm.C5_HUB_SEP_MIN:
        if 'C5hub' not in op and zv < HUB_MIN[0]: return 'Additive', 'C5'
        return 'Non-additive', 'C5'
    if linked and g('C5b', zv > dm.C5_LINKED_COMM_COUPLING) and seq > dm.C5_LINKED_COMM_SEP_MIN:
        return 'Non-additive', 'C5'
    if ch1 and ch2 and both_iso and seq >= 8:
        return 'Non-additive', 'C6'
    return 'Additive', 'none'



SWEEPS=[('C1_COUPLING',0.30,np.round(np.arange(-0.50,1.051,0.05),3)),
 ('C1_HELIX_FRAC',0.75,np.round(np.arange(0.50,1.101,0.02),3)),
 ('C2_COUPLING',0.45,np.round(np.arange(-0.50,1.001,0.05),3)),
 ('C3_SHORT_HELIX_COUPLING',0.40,np.round(np.arange(-0.25,0.951,0.05),3)),
 ('C3_HELIX_BOTH_ISO_COUPLING',1.00,np.round(np.arange(0.0,1.601,0.05),3)),
 ('C3_COIL_COUPLING',1.20,np.round(np.arange(0.80,1.751,0.05),3)),
 ('C3_HELIX_ONE_ISO_COUPLING',1.40,np.round(np.arange(0.90,2.001,0.05),3)),
 ('C3_SHORT_HELIX_SEP',4,np.arange(-3,16,1)),
 ('C3_HELIX_ISO_SEP_MAX',30,np.arange(0,121,2)),
 ('C3_HELIX_BOTH_ISO_SEP_MIN',4,np.arange(-5,36,1)),
 ('C3_HELIX_BOTH_ISO_SEP_MAX',30,np.arange(0,121,2)),
 ('C4_BACKBONE_CORR',0.55,np.round(np.arange(0.25,0.821,0.01),3)),
 ('C5_COMMON_NEIGHBOR_COUPLING',0.65,np.round(np.arange(-0.50,2.101,0.05),3)),
 ('C5_LINKED_COMM_COUPLING',1.00,np.round(np.arange(-0.50,3.301,0.05),3)),
 ('C5_HUB_BETWEENNESS',0.03,np.round(np.arange(0.020,0.0451,0.0005),4)),
 ('C5_HUB_SEP_MIN',15,np.arange(13,19,1)),
 ('C5_LINKED_COMM_SEP_MIN',30,np.arange(5,48,1)),
 ('C5_HUB_MIN_COUPLING',0.15,np.round(np.arange(0.0,0.501,0.01),3))]
ddg=pd.read_csv(DDG,sep='\t',on_bad_lines='skip'); ddg.columns=ddg.columns.str.strip()
cache=[]
for pdb in ['1STN','1BNI','2LZM','1PGA','1CSP','2RN2','2CI2']:
    c=dm.PROTEINS[pdb]; print('  loading',pdb,flush=True)
    net=dm.build_network(c['node_file'],c['edge_file']); md=dm.compute_md_features(c['gro'],c['xtc'],c['fps'])
    for _,r in ddg[ddg.PDB==pdb].iterrows():
        parts=str(r['Mutation_Double']).split(',')
        n=[int(''.join(ch for ch in p if ch.isdigit())) for p in parts if any(ch.isdigit() for ch in p)]
        if len(n)!=2: continue
        cache.append((n,r['Mutation_Double'],int(abs(float(r['dddG']))>=dm.DDDG_CUT),net,md))
def run():
    y=[];p=[]
    for n,mut,lab,(comm,r2c,direct,G_nb,bc,linked,last),(res,rmap,dccm,z,rmsf,ss) in cache:
        pr,_=predict_open(n[0],n[1],mut,comm,r2c,G_nb,bc,direct,linked,rmap,ss,dccm,z,set())
        y.append(lab); p.append(int(pr=='Non-additive'))
    y=np.array(y); p=np.array(p)
    tp=((y==1)&(p==1)).sum();tn=((y==0)&(p==0)).sum();fp=((y==0)&(p==1)).sum();fn=((y==1)&(p==0)).sum()
    den=np.sqrt(float((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn))); return (tp*tn-fp*fn)/den if den else 0.0
base=run(); print(f'Pairs: {len(cache)}  baseline MCC = {base:.3f}')
rows=[];summ=[]
for name,adopt,vals in SWEEPS:
    if name=='C5_HUB_MIN_COUPLING': old=HUB_MIN[0]
    else: old=getattr(dm,name)
    ms=[]
    for v in vals:
        if name=='C5_HUB_MIN_COUPLING': HUB_MIN[0]=float(v)
        else: setattr(dm,name,type(old)(v) if isinstance(old,int) else float(v))
        m=run(); ms.append(m); rows.append(dict(threshold=name,value=float(v),mcc=m))
    if name=='C5_HUB_MIN_COUPLING': HUB_MIN[0]=old
    else: setattr(dm,name,old)
    ms=np.array(ms); vals=np.array(vals,dtype=float); same=np.isclose(ms,base,atol=1e-9)
    i=int(np.argmin(abs(vals-adopt))); lo=hi=i
    while lo>0 and same[lo-1]: lo-=1
    while hi<len(vals)-1 and same[hi+1]: hi+=1
    brk=[float((vals[k]+vals[k+1])/2) for k in range(len(vals)-1) if not np.isclose(ms[k],ms[k+1],atol=1e-9)]
    summ.append(dict(threshold=name,adopted=adopt,plateau_lo=vals[lo],plateau_hi=vals[hi],mcc_min=ms.min(),mcc_max=ms.max(),
                     n_breaks=len(brk),breaks=';'.join(f'{b:g}' for b in brk)))
    print(f"{name:28s} adopted={adopt:<6g} plateau {vals[lo]:g} to {vals[hi]:g}   MCC range {ms.min():.3f} to {ms.max():.3f}")
pd.DataFrame(rows).to_csv(os.path.join(H,'figS2_sweep.csv'),index=False)
pd.DataFrame(summ).to_csv(os.path.join(H,'figS2_summary.csv'),index=False)
print('Saved figS2_sweep.csv, figS2_summary.csv')
