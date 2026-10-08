#!/usr/bin/env python3
"""Figure S3: MCC heatmap of the 135-cell virtual-edge grid, same design as the submitted Figure S3.
Reads fig4_A_grid.csv (copy of 5_results/virtual_edge_grid_265.csv made on this system by fig4_data.py).
Prints the grid for checking, then writes FigureS3.tiff / .png (300 dpi) and .pdf (editable text).
Run inside 7_figure_scripts/:  python3 plot_figS3.py"""
import os, numpy as np, pandas as pd, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
from matplotlib.patches import Rectangle
plt.rcParams.update({'font.family':'sans-serif','font.sans-serif':['Arial','Liberation Sans','DejaVu Sans'],
                     'pdf.fonttype':42,'ps.fonttype':42,'font.size':10})
H=os.path.dirname(os.path.abspath(__file__))
g=pd.read_csv(os.path.join(H,'fig4_A_grid.csv'))
M=g.pivot(index='dist',columns='zthr',values='mcc').sort_index(ascending=False)
pd.set_option('display.width',250); print(M.round(3).to_string()); 
b=g.loc[g.mcc.idxmax()]; print(f"\nBest cell: Z>{b.zthr}, d<={b.dist:.0f} A, MCC={b.mcc:.3f}; min MCC={g.mcc.min():.3f}")
cmap=LinearSegmentedColormap.from_list('bo',['#5b8ec4','#f2f2f2','#e6550d'])
norm=TwoSlopeNorm(vmin=min(g.mcc.min(),-0.01),vcenter=0.0,vmax=g.mcc.max())
fig,ax=plt.subplots(figsize=(13,6.2))
im=ax.imshow(M.values,cmap=cmap,norm=norm,aspect='auto')
for i in range(M.shape[0]):
    for j in range(M.shape[1]):
        ax.text(j,i,('0.00' if abs(M.values[i,j])<0.005 else f'{M.values[i,j]:.2f}'),ha='center',va='center',fontsize=7.5,color='#222222')
ax.set_xticks(range(M.shape[1])); ax.set_xticklabels([f'{v:g}' for v in M.columns],rotation=45)
ax.set_yticks(range(M.shape[0])); ax.set_yticklabels([f'{v:g}' for v in M.index])
ax.set_xticks(np.arange(-0.5,M.shape[1]),minor=True); ax.set_yticks(np.arange(-0.5,M.shape[0]),minor=True)
ax.grid(which='minor',color='white',lw=1.2); ax.tick_params(which='minor',length=0)
bi=list(M.index).index(b.dist); bj=list(M.columns).index(b.zthr)
ax.add_patch(Rectangle((bj-0.5,bi-0.5),1,1,fill=False,ec='black',lw=2.5,zorder=10,clip_on=False))
ax.set_xlabel('|Z-DCCM| threshold',fontsize=12); ax.set_ylabel('Cα distance cut-off (Å)',fontsize=12)
for s in ['top','right']: ax.spines[s].set_visible(False)
cb=fig.colorbar(im,ax=ax,fraction=0.03,pad=0.03,shrink=0.45,anchor=(0,0.8)); cb.ax.set_title('MCC',fontsize=11,pad=6)
for ext in ['tiff','png','pdf']:
    fig.savefig(os.path.join(H,f'FigureS3.{ext}'),dpi=300,bbox_inches='tight',**({'pil_kwargs':{'compression':'tiff_lzw'}} if ext=='tiff' else {}))
print('Saved FigureS3.tiff, FigureS3.png, FigureS3.pdf in',H)
