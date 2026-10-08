#!/usr/bin/env python3
"""Figure 1: the DynRCNC decision cascade (schematic, no data).
Run:  python3 make_figure1.py     -> Figure1.tiff (600 dpi, RGB, LZW), Figure1.pdf (vector), Figure1.png (600 dpi)
Every condition and threshold is taken from Section 2.4 and Table 2 of the manuscript."""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Circle, FancyArrowPatch
from matplotlib import font_manager
import os
from PIL import Image
OUT = os.path.dirname(os.path.abspath(__file__))
plt.rcParams['font.family'] = 'Liberation Sans' if any('Liberation Sans' in f.name for f in font_manager.fontManager.ttflist) else 'DejaVu Sans'
plt.rcParams['pdf.fonttype'] = 42
W = 7.2; U = 20.0                      # inches wide; drawing units per inch (the figure is printed at 6.5 in wide)
TEAL, TEAL_L, TEAL_D = '#2a8c82', '#e3f2ef', '#1d6b63'
ORG, ORG_L, ORG_D = '#d9962b', '#fdf1d6', '#8a5a0c'
NAVY, GREY, GREY_L = '#1f3a6e', '#6e6e6e', '#f1f1f1'
FS, FS_B, FS_H, FS_T = 7.0, 7.4, 8.2, 6.8
LH = 3.1
CX0, CX1 = 36.0, 143.0; TXT0 = CX0 + 8.0; TAG1 = CX1 - 2.6
# ------------------------------------------------------------------ content
T1 = [('C1', 'Clique community',
       ['Both sites in the same k = 3 community and |Z-DCCM| > 0.30',
        'Exception: helix-dominant community (helix fraction ≥ 0.75), or helix-strand boundary pair'],
       [('Non-additive', 'na', 0), ('Additive', 'ad', 2)])]
T2 = [('C2', 'Direct contact',
       ['RING contact between sites more than 2 residues apart, |Z-DCCM| > 0.45; adjacent helix pairs excluded'],
       [('Non-additive', 'na', 0)]),
      ('C3', 'Dynamic packing',
       ['Coil pair, |Z-DCCM| > 1.20; helix pair with one isolated site, > 1.40; with both sites isolated,',
        '> 1.00 and a side-chain volume change; short-range helix pair in adjacent communities, > 0.40'],
       [('Non-additive', 'na', 0)]),
      ('C4', 'Sequential backbone',
       ['Sequence separation ≤ 2, raw DCCM > 0.55, both sites network-isolated'],
       [('Non-additive', 'na', 0)]),
      ('C5', 'Network topology',
       ['Hub (betweenness > 0.03, separation > 15); shared non-backbone neighbor',
        '(|Z-DCCM| > 0.65, separation > 5); or linked communities (|Z-DCCM| > 1.00)',
        'A hub pair with |Z-DCCM| < 0.15 is returned as additive'],
       [('Non-additive', 'na', 0), ('Additive', 'ad', 3)])]
T3 = [('C6', 'Electrostatic coupling',
       ['Both sites retain a charge after mutation, both are network-isolated, separation ≥ 8'],
       [('Non-additive', 'na', 0)])]
def card_h(c): return 2.2 + LH * (1 + len(c[2])) + 0.6
def panel_h(cards): return 6.4 + sum(card_h(c) for c in cards) + 1.3 * (len(cards) - 1) + 1.8
GAP, DEF_H, MARG = 4.4, 7.0, 2.0
hT1, hT2, hT3 = panel_h(T1), panel_h(T2), panel_h(T3)
TOTAL = MARG + hT1 + GAP + hT2 + GAP + hT3 + GAP + DEF_H + MARG
H = TOTAL / U
fig = plt.figure(figsize=(W, H)); ax = fig.add_axes([0, 0, 1, 1]); ax.set_xlim(0, W * U); ax.set_ylim(0, TOTAL); ax.axis('off')
def box(x0, y0, x1, y1, fc, ec, lw=0.9, r=1.4, z=1): ax.add_patch(FancyBboxPatch((x0, y0), x1 - x0, y1 - y0, boxstyle=f'round,pad=0,rounding_size={r}', fc=fc, ec=ec, lw=lw, zorder=z))
def arrow(x0, y0, x1, y1, color='#444', lw=1.3, ms=9): ax.add_patch(FancyArrowPatch((x0, y0), (x1, y1), arrowstyle='-|>', mutation_scale=ms, color=color, lw=lw, zorder=3, shrinkA=0, shrinkB=0))
def tag(xr, yc, text, kind):
    w = 2.6 + 1.28 * len(text); fc, tc, ec = (NAVY, 'white', NAVY) if kind == 'na' else (GREY_L, '#333', GREY)
    box(xr - w, yc - 1.9, xr, yc + 1.9, fc, ec, lw=0.8, r=1.2, z=4); ax.text(xr - w / 2, yc, text, ha='center', va='center', fontsize=FS_T, color=tc, fontweight='bold', zorder=5)
def card(ytop, c, ec):
    code, title, lines, tags = c; h = card_h(c)
    box(CX0 + 1.5, ytop - h, CX1 - 1.5, ytop, 'white', ec, lw=0.9, r=1.2, z=2)
    y0 = ytop - 3.5
    ax.add_patch(Circle((CX0 + 5.2, y0), 2.35, fc=ec, ec='none', zorder=4)); ax.text(CX0 + 5.2, y0, code, ha='center', va='center', color='white', fontsize=FS_B, fontweight='bold', zorder=5)
    ax.text(TXT0, y0, title, ha='left', va='center', fontsize=FS_B, fontweight='bold', color='#222', zorder=5)
    for k, ln in enumerate(lines, 1): ax.text(TXT0, y0 - LH * k, ln, ha='left', va='center', fontsize=FS, color='#222', zorder=5)
    for text, kind, idx in tags: tag(TAG1, y0 + (0.5 if idx == 0 else 0) - LH * idx, text, kind)
    return h
def tier(ytop, h, title, cards, fc, ec, hdr, note=None):
    box(CX0, ytop - h, CX1, ytop, fc, ec, lw=1.3, r=2.0, z=0)
    ax.text(CX0 + 2.6, ytop - 3.1, title, ha='left', va='center', fontsize=FS_H, fontweight='bold', color=hdr, zorder=5)
    if note: ax.text(CX1 - 2.6, ytop - 3.1, note, ha='right', va='center', fontsize=FS, color=hdr, style='italic', zorder=5)
    yc = ytop - 6.4
    for c in cards: yc -= card(yc, c, ec) + 1.3
yT1 = TOTAL - MARG; yT2 = yT1 - hT1 - GAP; yT3 = yT2 - hT2 - GAP; yD = yT3 - hT3 - GAP
tier(yT1, hT1, 'Tier 1   Static topology', T1, TEAL_L, TEAL, TEAL_D)
tier(yT2, hT2, 'Tier 2   MD-gated conditional filters', T2, TEAL_L, TEAL, TEAL_D, note='tested in the order C2, C3, C4, C5')
tier(yT3, hT3, 'Tier 3   Electrostatics', T3, ORG_L, ORG, ORG_D)
box(CX0, yD - DEF_H, CX1, yD, GREY_L, GREY, lw=1.2, r=2.0, z=1)
ax.text(CX0 + 3.0, yD - DEF_H / 2, 'No criterion satisfied', ha='left', va='center', fontsize=FS_B, fontweight='bold', color='#222', zorder=5)
tag(TAG1 + 0.4, yD - DEF_H / 2, 'Additive', 'ad'); ax.text(TAG1 + 0.4 - 2.6 - 1.28 * 8 - 1.5, yD - DEF_H / 2, 'default prediction', ha='right', va='center', fontsize=FS, color='#444', zorder=5)
for ya, yb in [(yT1 - hT1, yT2), (yT2 - hT2, yT3), (yT3 - hT3, yD)]:
    xm = CX0 + 12.0; arrow(xm, ya, xm, yb); ax.text(xm + 1.8, (ya + yb) / 2, 'not satisfied', ha='left', va='center', fontsize=FS - 0.4, color='#444', style='italic', zorder=5)
# input panel
IX0, IX1 = 1.0, 31.0; ymid = (TOTAL - MARG + (yD - DEF_H)) / 2; ih = 62.0
box(IX0, ymid - ih / 2, IX1, ymid + ih / 2, '#f6f6f6', '#888', lw=1.1, r=2.0, z=0)
xc = (IX0 + IX1) / 2; ytop_i = ymid + ih / 2
ax.text(xc, ytop_i - 3.4, 'Input', ha='center', va='center', fontsize=FS_H, fontweight='bold', color='#222')
icy = ytop_i - 10.8
ax.plot([9.6, 22.4], [icy, icy], ls=(0, (2, 2)), color='#e07a00', lw=1.3, zorder=2)
ax.add_patch(Circle((9.6, icy), 2.7, fc='white', ec='#444', lw=1.1, zorder=3)); ax.add_patch(Circle((22.4, icy), 2.7, fc='#c0189c', ec='#444', lw=1.1, zorder=3))
ax.text(9.6, icy, 'i', ha='center', va='center', fontsize=FS_B, style='italic', zorder=4); ax.text(22.4, icy, 'j', ha='center', va='center', fontsize=FS_B, style='italic', color='white', zorder=4)
items = [('Double-mutant pair', True), ('(sites i and j)', False), ('', False), ('Native contact network:', True), ('RING 4.0, k = 3 clique', False), ('communities', False), ('', False), ('Wild-type MD:', True), ('|Z-DCCM| of the pair', False)]
for k, (ln, b) in enumerate(items): ax.text(xc, icy - 7.0 - k * LH, ln, ha='center', va='center', fontsize=FS, color='#222', fontweight='bold' if b else 'normal')
y_note = icy - 7.0 - len(items) * LH - 0.8
for k, ln in enumerate(['The first criterion', 'satisfied decides the', 'prediction and is', 'reported with it']): ax.text(xc, y_note - k * (LH - 0.3), ln, ha='center', va='center', fontsize=FS - 0.4, color='#444', style='italic')
arrow(IX1, ymid, CX0, ymid)
for ext, kw in [('png', dict(dpi=600)), ('pdf', dict())]: fig.savefig(os.path.join(OUT, f'Figure1.{ext}'), facecolor='white', **kw)
im = Image.open(os.path.join(OUT, 'Figure1.png')).convert('RGB'); im.save(os.path.join(OUT, 'Figure1.tiff'), compression='tiff_lzw', dpi=(600, 600))
print('Figure1 written:', im.size, 'pixels =', round(im.size[0] / 600, 2), 'x', round(im.size[1] / 600, 2), 'in at 600 dpi; printed at 6.5 in wide the height is', round(6.5 * H / W, 2), 'in')
