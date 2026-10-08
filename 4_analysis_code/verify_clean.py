import mdtraj as md, numpy as np, os
HERE=os.path.dirname(os.path.abspath(__file__))
print(f"{'protein':8}{'frames':>8}{'max_jump(nm)':>14}{'RMSF(nm)':>10}  status")
for p in ['1stn','1bni','2lzm','1pga','1csp','2rn2','2ci2']:
    gro=f'{HERE}/{p}_md/md_first_frame.gro'
    xtc=f'{HERE}/{p}_md/md_reduced.xtc'
    t=md.load(xtc,top=gro)
    ca=t.topology.select('name CA')
    tca=t.atom_slice(ca); tca.superpose(tca,0)
    c=tca.xyz
    jump=np.linalg.norm(np.diff(c,axis=0),axis=2).max()
    rmsf=c.std(axis=0).mean(axis=1).mean()
    st='CLEAN' if jump<2 else 'STILL BROKEN'
    print(f"{p:8}{t.n_frames:>8}{jump:>14.3f}{rmsf:>10.4f}  {st}")
