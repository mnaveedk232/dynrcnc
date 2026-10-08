#!/usr/bin/env python3
"""Figure 4 (clean 265 pairs), same layout as the submitted Figure 4.
Reads fig4_A_grid.csv, fig4_A_markers.csv, fig4_B.csv, fig4_C.csv.
Writes Figure4.tiff / .png (300 dpi) and .pdf (editable text).
Run inside 7_figure_scripts/:  python3 plot_fig4.py"""
import os, numpy as np, pandas as pd, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
plt.rcParams.update({'font.family':'sans-serif','font.sans-serif':['Arial','Liberation Sans','DejaVu Sans'],
                     'pdf.fonttype':42,'ps.fonttype':42,'font.size':11})
H=os.path.dirname(os.path.abspath(__file__))
g=pd.read_csv(os.path.join(H,'fig4_A_grid.csv')); mk=pd.read_csv(os.path.join(H,'fig4_A_markers.csv')).set_index('model')
B=pd.read_csv(os.path.join(H,'fig4_B.csv')); C=pd.read_csv(os.path.join(H,'fig4_C.csv'))
DYN='#2e7d32'; VE='#e8772e'; STA='#8c8c8c'
fig=plt.figure(figsize=(10.5,10))
gs=fig.add_gridspec(2,2,height_ratios=[1,1],hspace=0.32,wspace=0.45)
a=fig.add_subplot(gs[0,0]); b=fig.add_subplot(gs[0,1]); c=fig.add_subplot(gs[1,:])
# A
sc=a.scatter(g.specificity,g.sensitivity,c=g.mcc,cmap='viridis',s=38,zorder=2,edgecolors='none')
cb=fig.colorbar(sc,ax=a,fraction=0.07,pad=0.12,shrink=0.55,anchor=(0,0.85)); cb.ax.set_title('MCC',fontsize=12,pad=8); cb.ax.tick_params(labelsize=10)
a.scatter(mk.loc['DynRCNC','spec'],mk.loc['DynRCNC','sens'],marker='D',s=130,color=DYN,edgecolor='black',linewidth=1,zorder=4)
a.scatter(mk.loc['Static RCNC','spec'],mk.loc['Static RCNC','sens'],marker='^',s=170,color='#a6a6a6',edgecolor='black',linewidth=1,zorder=4)
a.scatter(mk.loc['Best virtual edge','spec'],mk.loc['Best virtual edge','sens'],marker='o',s=150,color=VE,edgecolor='#b03a10',linewidth=1.5,zorder=4)
h=[Line2D([],[],marker='D',ls='',ms=9,mfc=DYN,mec='black',label='DynRCNC'),
   Line2D([],[],marker='^',ls='',ms=11,mfc='#a6a6a6',mec='black',label='Static RCNC'),
   Line2D([],[],marker='o',ls='',ms=10,mfc=VE,mec='#b03a10',label='Best virtual edge')]
a.legend(handles=h,loc='lower left',frameon=False,fontsize=10,labelspacing=1.1,borderpad=0.8)
a.set_xlim(-0.04,1.05); a.set_ylim(-0.04,1.05)
a.set_xticks([0,0.25,0.5,0.75,1]); a.set_yticks([0,0.25,0.5,0.75,1])
a.set_xticklabels(['0.00','0.25','0.50','0.75','1.00']); a.set_yticklabels(['0.00','0.25','0.50','0.75','1.00'])
a.set_xlabel('Specificity',fontsize=13); a.set_ylabel('Sensitivity',fontsize=13)
# B
b.plot(B.zthr,B.best_mcc,'-o',color=VE,lw=2,ms=7,zorder=3)
dyn=mk.loc['DynRCNC','mcc']; sta=mk.loc['Static RCNC','mcc']; i=B.best_mcc.idxmax()
b.axhline(dyn,color=DYN,ls=(0,(1,2.5)),lw=2); b.axhline(sta,color=STA,ls='--',lw=1.5)
b.text(1.0,dyn+0.018,f'DynRCNC ({dyn:.3f})',color=DYN,fontweight='bold',fontsize=11)
b.text(3.05,sta+0.012,f'Static RCNC ({sta:.3f})',color=STA,fontsize=10)
b.text(B.zthr[i],B.best_mcc[i]+0.02,f'peak {B.best_mcc[i]:.3f}',color=VE,fontweight='bold',fontsize=11,ha='center')
b.set_ylim(-0.02,0.65); b.set_xlim(0.3,4.2); b.set_xticks([1,2,3,4]); b.set_yticks([0,0.2,0.4,0.6])
b.set_yticklabels(['0.0','0.2','0.4','0.6'])
b.set_xlabel('|Z-DCCM| threshold',fontsize=13); b.set_ylabel('Best MCC (across distance cut-offs)',fontsize=13)
# C
rng=np.random.default_rng(42); grp=[('Non-additive','#c1e3bd'),('Additive','#d9d9d9')]
for k,(lab,col) in enumerate(grp):
    v=C[C.exp==lab].zdccm.values
    c.boxplot(v,positions=[k],widths=0.34,patch_artist=True,showfliers=False,
              boxprops=dict(facecolor=col,edgecolor='#404040',lw=1.3),medianprops=dict(color='#404040',lw=2.5),
              whiskerprops=dict(color='#404040',lw=1),capprops=dict(lw=0))
    c.scatter(k+rng.uniform(-0.1,0.1,len(v)),v,s=10,color='#7f7f7f',alpha=0.55,zorder=3,edgecolors='#555555',linewidths=0.3)
c.set_xticks([0,1]); c.set_xticklabels([g_[0] for g_ in grp],fontsize=12); c.set_xlim(-0.55,1.55)
c.set_ylabel('|Z-DCCM|',fontsize=13); c.set_ylim(-0.15,3.5); c.set_yticks([0,1,2,3])
for ax,lab in [(a,'A'),(b,'B'),(c,'C')]:
    ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False); ax.tick_params(labelsize=11)
    ax.text(-0.14 if ax is not c else -0.06,1.04,lab,transform=ax.transAxes,fontsize=18,fontweight='bold',va='bottom')
for ext in ['tiff','png','pdf']:
    fig.savefig(os.path.join(H,f'Figure4.{ext}'),dpi=300,bbox_inches='tight',**({'pil_kwargs':{'compression':'tiff_lzw'}} if ext=='tiff' else {}))
print('Saved Figure4.tiff, Figure4.png, Figure4.pdf in',H)
