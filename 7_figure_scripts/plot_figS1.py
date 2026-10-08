#!/usr/bin/env python3
"""Figure S1: C-alpha RMSD per protein (7 panels), same design as the submitted Figure S1.
Reads figS1_rmsd.csv (made on this system by figS1_data.py).
Writes FigureS1.tiff / .png (300 dpi) and .pdf (editable text).
Run inside 7_figure_scripts/:  python3 plot_figS1.py"""
import os, numpy as np, pandas as pd, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams.update({'font.family':'sans-serif','font.sans-serif':['Arial','Liberation Sans','DejaVu Sans'],
                     'pdf.fonttype':42,'ps.fonttype':42,'font.size':10})
H=os.path.dirname(os.path.abspath(__file__))
d=pd.read_csv(os.path.join(H,'figS1_rmsd.csv'))
order=['1STN','1BNI','2LZM','1PGA','1CSP','2RN2','2CI2']
ymax=max(3.2,np.ceil(d.rmsd_A.max()*10)/10+0.1)
fig,axs=plt.subplots(2,4,figsize=(13.8,7.2),sharex=True,sharey=True)
for ax,p in zip(axs.flat,order):
    g=d[d.PDB==p]; last=g[g.time_ns>=g.time_ns.max()-100].rmsd_A.mean()
    ax.plot(g.time_ns,g.rmsd_A,color='#3182bd',lw=0.6)
    ax.axhline(last,color='#e6550d',ls=(0,(4,3)),lw=1.1)
    ax.set_xlim(-2,202); ax.set_ylim(0,ymax); ax.set_xticks([0,50,100,150,200]); ax.set_yticks([0,1,2,3])
    ax.grid(color='#e5e5e5',lw=0.8); ax.set_axisbelow(True)
    for s in ax.spines.values(): s.set_color('#333333'); s.set_linewidth(0.8)
    ax.tick_params(colors='#333333',labelsize=10,length=3)
    from matplotlib.patches import Rectangle
    ax.add_patch(Rectangle((0,1.0),1,0.12,transform=ax.transAxes,fc='#262626',ec='#262626',clip_on=False))
    ax.text(0.5,1.06,p,transform=ax.transAxes,ha='center',va='center',color='white',fontsize=12,fontweight='bold')
axs.flat[-1].set_visible(False)
axs.flat[3].tick_params(labelbottom=True)
fig.subplots_adjust(left=0.065,right=0.99,bottom=0.09,top=0.93,wspace=0.06,hspace=0.28)
fig.text(0.53,0.02,'Time (ns)',ha='center',fontsize=13)
fig.text(0.018,0.51,'Cα-backbone RMSD (Å)',va='center',rotation=90,fontsize=13)
for ext in ['tiff','png','pdf']:
    fig.savefig(os.path.join(H,f'FigureS1.{ext}'),dpi=300,bbox_inches='tight',**({'pil_kwargs':{'compression':'tiff_lzw'}} if ext=='tiff' else {}))
print('Saved FigureS1.tiff, FigureS1.png, FigureS1.pdf in',H)
