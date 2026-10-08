#!/usr/bin/env python3
"""Figure S2: one-at-a-time threshold sensitivity (18 panels), same design as the submitted Figure S2.
Reads figS2_sweep.csv and figS2_summary.csv (made on this system by figS2_data.py).
Writes FigureS2.tiff / .png (300 dpi) and .pdf (editable text).
Run inside 7_figure_scripts/:  python3 plot_figS2.py"""
import os, string, numpy as np, pandas as pd, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
plt.rcParams.update({'font.family':'sans-serif','font.sans-serif':['Arial','Liberation Sans','DejaVu Sans'],
                     'pdf.fonttype':42,'ps.fonttype':42,'font.size':9})
H=os.path.dirname(os.path.abspath(__file__))
sw=pd.read_csv(os.path.join(H,'figS2_sweep.csv')); sm=pd.read_csv(os.path.join(H,'figS2_summary.csv'))
BASE=sw[(sw.threshold=='C1_COUPLING')&np.isclose(sw.value,0.30)].mcc.iloc[0]
PHYS={'C1_HELIX_FRAC':1.0}
fig,axs=plt.subplots(6,3,figsize=(11,14.5))
for k,(ax,(_,s)) in enumerate(zip(axs.flat,sm.iterrows())):
    g=sw[sw.threshold==s.threshold].sort_values('value'); x=g.value.values; y=g.mcc.values
    ax.axvspan(s.plateau_lo,s.plateau_hi,color='#cfe3c4',zorder=0)
    brk=[float(b) for b in str(s.breaks).split(';') if b not in ('','nan')]
    left=[b for b in brk if b<s.plateau_lo]; right=[b for b in brk if b>s.plateau_hi]
    for b in ([max(left)] if left else [])+([min(right)] if right else []):
        ax.axvline(b,color='#f08080',lw=1.2,zorder=1)
    if s.threshold in PHYS: ax.axvline(PHYS[s.threshold],color='#d98c1f',ls=(0,(1,1.5)),lw=1.4,zorder=2)
    ax.axhline(BASE,color='#8c8c8c',ls='--',lw=0.9,zorder=1)
    ax.plot(x,y,color='#2166ac',lw=1.6,zorder=3)
    ai=int(np.argmin(abs(x-s.adopted)))
    ax.scatter([s.adopted],[y[ai]],s=34,color='#e0402a',edgecolor='#b22222',zorder=4)
    ax.set_ylim(0.30,0.62); ax.set_yticks(np.arange(0.30,0.61,0.05))
    ax.set_yticklabels([f'{v:.2f}' for v in np.arange(0.30,0.61,0.05)],fontsize=8)
    pad=(x.max()-x.min())*0.04; ax.set_xlim(x.min()-pad,x.max()+pad); ax.tick_params(axis='x',labelsize=8)
    ax.grid(axis='y',color='#ebebeb',lw=0.6); ax.set_axisbelow(True)
    ax.set_xlabel('threshold value',fontsize=9,labelpad=1); ax.set_ylabel('MCC',fontsize=9,labelpad=1)
    ax.text(0.01,1.02,f'{k+1}. {s.threshold}',transform=ax.transAxes,fontsize=7.5,fontweight='bold',va='bottom')
    ax.text(0.99,1.02,f'adopted = {s.adopted:g}',transform=ax.transAxes,fontsize=7.5,color='#707070',ha='right',va='bottom')
    ax.text(-0.2,1.08,string.ascii_uppercase[k],transform=ax.transAxes,fontsize=13,fontweight='bold',va='bottom')
    for sp in ax.spines.values(): sp.set_linewidth(0.8)
h=[Line2D([],[],color='#2166ac',lw=1.6,label='MCC vs threshold value'),Patch(color='#cfe3c4',label='Stable plateau (MCC unchanged)'),
   Line2D([],[],color='#d98c1f',ls=(0,(1,1.5)),lw=1.4,label='Physical maximum'),
   Line2D([],[],marker='o',ls='',ms=6,mfc='#e0402a',mec='#b22222',label='Adopted value'),
   Line2D([],[],color='#f08080',lw=1.2,label='Break point (MCC changes)'),Line2D([],[],color='#8c8c8c',ls='--',lw=0.9,label=f'Baseline MCC {BASE:.3f}')]
fig.legend(handles=h,loc='lower center',ncol=3,frameon=False,fontsize=9,bbox_to_anchor=(0.5,0.0))
plt.tight_layout(rect=(0,0.04,1,1),h_pad=1.8,w_pad=1.6)
for ext in ['tiff','png','pdf']:
    fig.savefig(os.path.join(H,f'FigureS2.{ext}'),dpi=300,bbox_inches='tight',**({'pil_kwargs':{'compression':'tiff_lzw'}} if ext=='tiff' else {}))
print('Saved FigureS2.tiff, FigureS2.png, FigureS2.pdf in',H)
