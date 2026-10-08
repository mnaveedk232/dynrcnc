import warnings; warnings.filterwarnings('ignore')
import sys, pandas as pd, numpy as np, networkx as nx
from networkx.algorithms.community import k_clique_communities
import mdtraj as md, MDAnalysis as mda

# Clean-benchmark version of ve_one_user.py (logic unchanged; only paths are arguments)
# usage: python3 ve_clean.py <code> <PDB> <fps> <data_dir> <ddg_file> <out_dir>
code=sys.argv[1]; PDB=sys.argv[2]; fps=int(sys.argv[3])
BASE=sys.argv[4]
FOLD=f'{BASE}/{code}_md'
DDG=sys.argv[5]
OUTDIR=sys.argv[6]
ZTHRS=[0.5,0.75,1.0,1.25,1.5,1.75,2.0,2.25,2.5,2.75,3.0,3.25,3.5,3.75,4.0]
DISTS=[6,8,10,12,14,16,18,20,25]
MAXADD=1200

def build_base_G(nf,ef):
    n_df=pd.read_table(nf,encoding='latin-1'); n_df.columns=n_df.columns.str.strip()
    e_df=pd.read_table(ef,encoding='latin-1'); e_df.columns=e_df.columns.str.strip()
    G=nx.Graph(); sd=set()
    for _,row in e_df.iterrows():
        s1=int(str(row['NodeId1']).split(':')[1]); s2=int(str(row['NodeId2']).split(':')[1])
        it=str(row['Interaction']).split(':')[0]
        if it=='HBOND':
            if abs(s1-s2)>4: sd.add(tuple(sorted((s1,s2))))
        else: sd.add(tuple(sorted((s1,s2))))
    for a,b in sd: G.add_edge(a,b)
    last=int(str(n_df['NodeId'].iloc[-1]).split(':')[1])
    for i in range(1,last): G.add_edge(i,i+1)
    return G,last

def zdccm_coords(gro,xtc,fps):
    u=mda.Universe(gro,xtc); resids=np.array([r.resid for r in u.residues])
    rmap={r.resid:k for k,r in enumerate(u.residues)}
    traj=md.load(xtc,top=gro); ca=traj.topology.select('name CA')
    tca=traj.atom_slice(ca); tca.superpose(tca,0); coords=tca.xyz
    delta=coords-coords.mean(axis=0)
    num=np.einsum('fic,fjc->ij',delta,delta); sq=np.sqrt(np.einsum('fic,fic->i',delta,delta))
    dccm=np.where(np.outer(sq,sq)>0,num/np.outer(sq,sq),0.0)
    z=np.zeros_like(dccm); n=len(resids)
    for d in range(int(resids.max()-resids.min())+1):
        pr=[(i,j) for i in range(n) for j in range(n) if i!=j and abs(resids[i]-resids[j])==d]
        if len(pr)<5: continue
        vals=[dccm[i,j] for i,j in pr]; mu,sg=np.mean(vals),max(np.std(vals),0.01)
        for i,j in pr: z[i,j]=(dccm[i,j]-mu)/sg
    xyz=traj.atom_slice(ca).xyz.mean(axis=0)*10.0
    res_order=[traj.topology.atom(a).residue.resSeq for a in ca]; rid2idx={r:i for i,r in enumerate(res_order)}
    return rmap,z,xyz,rid2idx

G,last=build_base_G(f'{FOLD}/{code}.N',f'{FOLD}/{code}.E')
base=list(G.edges())
rmap,z,xyz,rid2idx=zdccm_coords(f'{FOLD}/md_first_frame.gro',f'{FOLD}/md_reduced.xtc',fps)
resids=list(rmap.keys()); cand=[]
for a in range(len(resids)):
    for b in range(a+1,len(resids)):
        ra,rb=resids[a],resids[b]; ia,ib=rmap[ra],rmap[rb]
        if ia>=z.shape[0] or ib>=z.shape[0]: continue
        za=abs(z[ia,ib])
        if za<=0.5: continue
        if ra in rid2idx and rb in rid2idx:
            dist=float(np.linalg.norm(xyz[rid2idx[ra]]-xyz[rid2idx[rb]]))
            cand.append((ra,rb,za,dist))
cand=np.array(cand) if cand else np.zeros((0,4))
ddg=pd.read_csv(DDG,sep='\t',on_bad_lines='skip'); ddg.columns=ddg.columns.str.strip()
pairs=[]
for _,row in ddg[ddg['PDB']==PDB].iterrows():
    parts=str(row['Mutation_Double']).split(','); nums=[int(''.join(c for c in p if c.isdigit())) for p in parts if any(c.isdigit() for c in p)]
    if len(nums)!=2: continue
    exp='Non-additive' if abs(float(row['dddG']))>=1.0 else 'Additive'
    pairs.append((nums[0],nums[1],exp))

out=[]
for zt in ZTHRS:
    for dc in DISTS:
        G2=nx.Graph(); G2.add_edges_from(base)
        if len(cand):
            sel=cand[(cand[:,2]>zt)&(cand[:,3]<=dc)]
            if len(sel)>MAXADD:
                sel=sel[np.argsort(-sel[:,2])[:MAXADD]]
            G2.add_edges_from([(int(r[0]),int(r[1])) for r in sel])
        comm=list(k_clique_communities(G2,3))
        r2c={}
        for ci,cc in enumerate(comm):
            for r in cc: r2c.setdefault(r,set()).add(ci)
        for s1,s2,exp in pairs:
            same=bool(r2c.get(s1,set()) & r2c.get(s2,set()))
            out.append((zt,dc,s1,s2,exp,'NA' if same else 'ADD'))
pd.DataFrame(out,columns=['zthr','dist','s1','s2','exp','pred']).to_csv(f'{OUTDIR}/ve2_{PDB}.csv',index=False)
print(f"{PDB}: {len(pairs)} pairs done, saved ve2_{PDB}.csv",flush=True)
