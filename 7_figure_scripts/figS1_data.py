#!/usr/bin/env python3
"""Figure S1 data: C-alpha RMSD (A) of each 200 ns production trajectory, clean benchmark
(1PGA and 2CI2 use the PBC-corrected trajectories in 00_clean_benchmark).
Superposed on C-alpha of the first production frame. Output: figS1_rmsd.csv
Run inside 7_figure_scripts/:  python3 figS1_data.py"""
import os, numpy as np, pandas as pd, mdtraj as md
H=os.path.dirname(os.path.abspath(__file__))
C=os.path.expanduser(os.environ.get('MD_BASE',os.path.join(os.path.dirname(os.path.abspath(__file__)),'..','data')))
out=[]
print(f"{'PDB':5s} {'frames':>6s} {'dt_ps':>6s} {'end_ns':>7s} {'mean_last100':>12s} {'max':>6s} {'mean_0-20ns':>11s}")
for code in ['1stn','1bni','2lzm','1pga','1csp','2rn2','2ci2']:
    f=os.path.join(C,f'{code}_md')
    t=md.load(os.path.join(f,f'{code}_md_reduced.xtc'),top=os.path.join(f,f'{code}_md_first_frame.gro'))
    ca=t.topology.select('name CA'); t=t.atom_slice(ca); t.superpose(t,0)
    r=md.rmsd(t,t,0)*10.0; ns=(t.time-t.time[0])/1000.0
    out.append(pd.DataFrame(dict(PDB=code.upper(),time_ns=ns,rmsd_A=r)))
    last=r[ns>=ns[-1]-100]
    print(f"{code.upper():5s} {t.n_frames:6d} {t.time[1]-t.time[0]:6.0f} {ns[-1]:7.1f} {last.mean():12.3f} {r.max():6.3f} {r[ns<=20].mean():11.3f}")
pd.concat(out).to_csv(os.path.join(H,'figS1_rmsd.csv'),index=False)
print('Saved figS1_rmsd.csv')
