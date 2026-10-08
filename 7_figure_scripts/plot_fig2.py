#!/usr/bin/env python3
"""Figure 2 (clean 265-pair benchmark), same layout as the submitted Figure 2.
Reads fig2_data.csv, writes Figure2.tiff / Figure2.png (300 dpi) and Figure2.pdf (editable text).
Run inside 7_figure_scripts/:  python3 plot_fig2.py"""
import os, pandas as pd, numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams.update({'font.family':'sans-serif','font.sans-serif':['Arial','Liberation Sans','DejaVu Sans'],
                     'pdf.fonttype':42,'ps.fonttype':42,'svg.fonttype':'none','font.size':11})
HERE=os.path.dirname(os.path.abspath(__file__))
d=pd.read_csv(os.path.join(HERE,'fig2_data.csv'))
DYN='#1d7a8a'; STA='#b9b1a8'; POS='#3a6ea5'; NEG='#c0504d'; ZERO='#707070'
x=np.arange(len(d)); w=0.36
fig,(a,b)=plt.subplots(1,2,figsize=(13,6))
# Panel A
a.bar(x-w/2,d.MCC_dyn,w,color=DYN,edgecolor='black',linewidth=0.6,label='DynRCNC',zorder=2)
a.bar(x+w/2,d.MCC_static,w,color=STA,edgecolor='black',linewidth=0.6,label='Static RCNC',zorder=2)
lo=d.MCC_dyn-d.CI_dyn_lo; hi=d.CI_dyn_hi-d.MCC_dyn
a.errorbar(x-w/2,d.MCC_dyn,yerr=[lo.clip(lower=0),hi.clip(lower=0)],fmt='none',ecolor='black',elinewidth=1,capsize=4,zorder=3)
for xi,(m,h,tp,na) in enumerate(zip(d.MCC_dyn,d.CI_dyn_hi,d.TP_dyn,d.NA)):
    top=max(m,h,0)
    a.text(xi-w/2,top+0.025,f'{tp}/{na}',ha='center',va='bottom',color=DYN,fontweight='bold',fontsize=11)
a.set_ylabel('Matthews correlation coefficient (MCC)',fontsize=13)
a.legend(loc='lower left',frameon=True,fancybox=False,edgecolor='#999999',fontsize=11,borderpad=0.8)
# Panel B
imp=d.Improvement.values
cols=[POS if v>0 else (NEG if v<0 else ZERO) for v in imp]
b.bar(x,imp,0.6,color=cols,edgecolor='black',linewidth=0.6,zorder=2)
for xi,v in enumerate(imp):
    lab=f'+{v:.2f}' if v>0 else f'{v:.2f}'
    col=POS if v>0 else (NEG if v<0 else ZERO)
    b.text(xi,v+0.025 if v>=0 else v-0.025,lab,ha='center',va='bottom' if v>=0 else 'top',color=col,fontweight='bold',fontsize=11)
b.set_ylabel('Improvement over Static  (ΔMCC)',fontsize=13)
for ax,lab in [(a,'A'),(b,'B')]:
    ax.set_xticks(x); ax.set_xticklabels(d.PDB,fontsize=11)
    ax.set_ylim(-0.48,1.3); ax.set_yticks(np.arange(-0.4,1.21,0.2))
    ax.set_yticklabels([f'{v:.1f}' for v in np.arange(-0.4,1.21,0.2)],fontsize=11)
    ax.axhline(0,color='black',linestyle=(0,(6,3,1,3)),linewidth=0.8,zorder=1)
    ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)
    ax.set_xlim(-0.6,len(d)-0.4)
    ax.text(-0.13,1.04,lab,transform=ax.transAxes,fontsize=18,fontweight='bold',va='bottom')
plt.tight_layout(w_pad=4)
for ext in ['tiff','png','pdf']:
    fig.savefig(os.path.join(HERE,f'Figure2.{ext}'),dpi=300,bbox_inches='tight',**({'pil_kwargs':{'compression':'tiff_lzw'}} if ext=='tiff' else {}))
print('Saved Figure2.tiff, Figure2.png, Figure2.pdf in',HERE)
