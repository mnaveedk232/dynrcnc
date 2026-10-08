#!/usr/bin/env python3
"""Figure 5 (clean 265 pairs): radial bar chart of in-sample vs LOPOCV MCC, same design as the
submitted Figure 5 (legend moved to a small corner box; no legend above the plot).
Reads fig5_data.csv. Writes Figure5.tiff / .png (300 dpi) and .pdf (editable text).
Run inside 7_figure_scripts/:  python3 plot_fig5.py"""
import os, numpy as np, pandas as pd, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
plt.rcParams.update({'font.family':'sans-serif','font.sans-serif':['Arial','Liberation Sans','DejaVu Sans'],
                     'pdf.fonttype':42,'ps.fonttype':42})
H=os.path.dirname(os.path.abspath(__file__))
d=pd.read_csv(os.path.join(H,'fig5_data.csv')); d['method']=d.method.str.replace('\\n','\n',regex=False)
YEL='#f2c661'; GRN='#1f8a70'; TXY='#c8962a'; TXG='#1b6b58'
n=len(d); th=np.deg2rad(90-np.arange(n)*360/n)      # clockwise from top
HOLE=0.30; SC=1.0                                    # radius = HOLE + value*SC
r=lambda v: HOLE+max(v,0)*SC
fig=plt.figure(figsize=(7.2,8)); ax=fig.add_subplot(111,projection='polar')
wi=np.deg2rad(360/n*0.62); wl=np.deg2rad(360/n*0.32)
for t,(m,vi,vl) in zip(th,zip(d.method,d.in_sample,d.lopocv)):
    if not np.isnan(vi): ax.bar(t,max(vi,0)*SC,width=wi,bottom=HOLE,color=YEL,zorder=2,linewidth=0)
    ax.bar(t,max(vl,0.012)*SC if vl>=0 else 0.012,width=wl,bottom=HOLE,color=GRN,zorder=3,linewidth=0)
    top=r(np.nanmax([vi,vl])); lr=top+0.10
    kw=dict(fontweight='bold',fontsize=11,ha='center',va='center',textcoords='offset points')
    if np.isnan(vi): ax.annotate(f'{vl:.3f}',(t,lr),xytext=(0,0),color=TXG,**kw)
    else:
        ax.annotate(f'{vl:.3f}',(t,lr),xytext=(0,7.5),color=TXG,**kw)
        ax.annotate(f'{vi:.3f}',(t,lr),xytext=(0,-7.5),color=TXY,**kw)
    ax.text(t,1.47 if abs(np.sin(t))>0.9 else 1.48,m,fontweight='bold',fontsize=12,ha='center',va='center')
for rr in [HOLE,HOLE+0.25,HOLE+0.5,HOLE+0.75,HOLE+0.95]:
    ax.plot(np.linspace(0,2*np.pi,400),[rr]*400,color='#d0d0d0',lw=1,zorder=1)
ax.plot(np.linspace(0,2*np.pi,100),[0.03]*100,color='#a0a0a0',lw=1)
ax.set_ylim(0,1.25); ax.set_xticks([]); ax.set_yticks([]); ax.spines['polar'].set_visible(False); ax.grid(False)
fig.legend(handles=[Patch(color=YEL,label='In-sample'),Patch(color=GRN,label='LOPOCV')],loc='upper left',
           bbox_to_anchor=(0.02,0.98),frameon=False,fontsize=10,handlelength=1.2)
for ext in ['tiff','png','pdf']:
    fig.savefig(os.path.join(H,f'Figure5.{ext}'),dpi=300,bbox_inches='tight',**({'pil_kwargs':{'compression':'tiff_lzw'}} if ext=='tiff' else {}))
print('Saved Figure5.tiff, Figure5.png, Figure5.pdf in',H)
