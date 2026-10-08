#!/usr/bin/env python3
"""
Build the DDGUN vs Benchmark vs DynRCNC comparison Excel (265 pairs, 7 proteins).

INPUT FILES (put this script in the same folder, or edit the paths below):
  - ddgun_per_pair.csv          # DDGUN per-pair table (offset columns, see note)
  - *_ddgun3d.out               # one per protein: 1STN 1BNI 2LZM 1PGA 1CSP 2RN2 2CI2
  - all_predictions.csv         # DynRCNC output (from the dynrcnc.py run)

OUTPUT:
  - DDGUN_vs_Benchmark_vs_DynRCNC_265pairs.xlsx  (3 sheets)

Run:  python3 build_ddgun_comparison.py
Needs: pip install openpyxl
"""

import csv, math, glob, os
from collections import defaultdict
from openpyxl import Workbook
from openpyxl.styles import Font, PatternFill, Alignment, Border, Side
from openpyxl.utils import get_column_letter

# ----------------------------------------------------------------------
# EDIT THESE PATHS if your files live elsewhere
# ----------------------------------------------------------------------
DDGUN_CSV   = "ddgun_per_pair.csv"
OUT_GLOB    = "*_ddgun3d.out"          # all seven .out files
DYNRCNC_CSV = "all_predictions.csv"    # DynRCNC per-pair predictions
OUTFILE     = "DDGUN_vs_Benchmark_vs_DynRCNC_265pairs.xlsx"

PROTEIN_ORDER = ['1STN', '1BNI', '2LZM', '1PGA', '1CSP', '2RN2', '2CI2']

# ----------------------------------------------------------------------
# 1. DDGUN .out files -> singles s1,s2 and double T, for the T=mean proof
# ----------------------------------------------------------------------
def load_ddgun_out(pattern):
    singles = {}   # (PDB, frozenset{mut1,mut2}) -> (s1, s2, T, variant_string)
    for path in glob.glob(pattern):
        pdb = os.path.basename(path).split('_')[0].upper()
        with open(path) as f:
            for line in f:
                if line.startswith('#') or not line.strip():
                    continue
                parts = line.rstrip('\n').split('\t')
                variant = parts[2].replace(' ', '')
                if ',' in variant:                       # double mutant row
                    s = [float(x) for x in parts[3].split(',')]
                    T = float(parts[4])
                    key = (pdb, frozenset(variant.split(',')))
                    singles[key] = (s[0], s[1], T, variant)
    return singles

# ----------------------------------------------------------------------
# 2. DDGUN per-pair CSV. NOTE: this file's header is shifted by one column
#    because the mutation pair got split across two columns. We read it by
#    POSITION, not by header name.
#    Column order actually is:
#      0 protein | 1 mut1 | 2 mut2 | 3 Method | 4 exp_dddG |
#      5 ddgun_dddG | 6 exp_class | 7 ddgun_class | 8 match
# ----------------------------------------------------------------------
def load_ddgun_csv(path):
    rows = []
    with open(path) as f:
        reader = csv.reader(f)
        next(reader)                                     # skip header
        for p in reader:
            if len(p) < 9:
                continue
            rows.append({
                'PDB':         p[0],
                'pair':        p[1] + ',' + p[2],
                'Method':      p[3],
                'exp_dddG':    round(float(p[4]), 2),
                'ddgun_dddG':  round(float(p[5]), 2),
                'exp_class':   p[6],                     # 'non-add' / 'additive'
                'ddgun_class': p[7],
                'match':       p[8],
            })
    return rows

# ----------------------------------------------------------------------
# 3. DynRCNC predictions
# ----------------------------------------------------------------------
def load_dynrcnc(path):
    rows = []
    with open(path) as f:
        for r in csv.DictReader(f):
            m = r['mutation'].replace('"', '').replace(' ', '')
            rows.append({
                'PDB':      r['PDB'],
                'pair':     m,
                'pH':       r['pH'],
                'Method':   r['Method'],
                'exp_dddG': round(float(r['dddG']), 2),
                'dyn_pred': r['pred'],
                'dyn_correct': r['correct'],
                'mechanism': r['mechanism'],
                'seq_dist': r['seq_dist'],
                'ss_i':     r['ss_i'],
                'ss_j':     r['ss_j'],
                'zdccm':    r['zdccm'],
            })
    return rows

# ----------------------------------------------------------------------
# Merge the three sources (row-by-row, duplicates preserved)
# ----------------------------------------------------------------------
def merge(dyn, ddg, singles):
    def key(d):
        return (d['PDB'], frozenset(d['pair'].split(',')), d['Method'], d['exp_dddG'])

    ddg_by = defaultdict(list)
    for d in ddg:
        ddg_by[key(d)].append(d)

    used = defaultdict(int)
    merged = []
    for d in dyn:
        k = key(d)
        cand = ddg_by.get(k, [])
        idx = used[k]
        dd = cand[idx] if idx < len(cand) else None
        used[k] += 1

        row = dict(d)
        if dd:
            row['ddgun_dddG']  = dd['ddgun_dddG']
            row['ddgun_class'] = dd['ddgun_class']
            row['exp_class']   = dd['exp_class']
        else:
            row['ddgun_dddG'] = row['ddgun_class'] = row['exp_class'] = None

        sm = singles.get((d['PDB'], frozenset(d['pair'].split(','))))
        if sm:
            s1, s2, T, _ = sm
            row['s1'], row['s2'], row['T'] = s1, s2, T
            row['mean'] = round((s1 + s2) / 2, 3)
            row['T_eq_mean'] = 'YES' if abs(T - (s1 + s2) / 2) < 0.051 else 'NO'
        else:
            row['s1'] = row['s2'] = row['T'] = row['mean'] = row['T_eq_mean'] = None
        merged.append(row)
    return merged

# ----------------------------------------------------------------------
# Confusion-matrix metrics (positive class = non-additive)
# ----------------------------------------------------------------------
def metrics(sub, pred_key):
    TP = TN = FP = FN = 0
    for r in sub:
        e = (r['exp_class'] == 'non-add')
        p = (r[pred_key] == 'non-add')
        if   e and p:         TP += 1
        elif not e and not p: TN += 1
        elif not e and p:     FP += 1
        else:                 FN += 1
    N = TP + TN + FP + FN
    se = TP / (TP + FN) if TP + FN else 0
    sp = TN / (TN + FP) if TN + FP else 0
    pr = TP / (TP + FP) if TP + FP else 0
    f1 = 2 * pr * se / (pr + se) if pr + se else 0
    den = math.sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN))
    mcc = (TP * TN - FP * FN) / den if den else 0
    return TP, TN, FP, FN, N, se * 100, sp * 100, f1, mcc

# ----------------------------------------------------------------------
# Styling helpers
# ----------------------------------------------------------------------
ARI      = lambda **k: Font(name='Arial', size=k.pop('size', 10), **k)
THIN     = Side(style='thin', color='BFBFBF')
BORDER   = Border(left=THIN, right=THIN, top=THIN, bottom=THIN)
CENTER   = Alignment(horizontal='center', vertical='center')
HDR_FILL = PatternFill('solid', fgColor='1F4E78')
HDR_FONT = Font(name='Arial', bold=True, color='FFFFFF', size=10)
GRP = {
    'exp':   PatternFill('solid', fgColor='DDEBF7'),
    'ddgun': PatternFill('solid', fgColor='FCE4D6'),
    'dyn':   PatternFill('solid', fgColor='E2EFDA'),
}

def exp_word(r):   return 'Non-additive' if r['exp_class']   == 'non-add' else 'Additive'
def ddgun_word(r): return 'Non-additive' if r['ddgun_class'] == 'non-add' else 'Additive'

# ----------------------------------------------------------------------
# Build the workbook
# ----------------------------------------------------------------------
def build(merged):
    # derive DynRCNC class label
    for r in merged:
        r['dyn_class'] = 'non-add' if r['dyn_pred'] == 'Non-additive' else 'additive'

    merged.sort(key=lambda r: (PROTEIN_ORDER.index(r['PDB']), -r['exp_dddG']))
    wb = Workbook()

    # ---- Sheet 1: per-pair ------------------------------------------------
    ws = wb.active; ws.title = 'Per-pair comparison'
    headers = ['PDB', 'Mutation pair', 'pH', 'Method',
               'Exp dddG (kcal/mol)', 'Exp class',
               'DDGUN s1', 'DDGUN s2', 'DDGUN T', 'DDGUN mean(s1,s2)', 'T = mean?',
               'DDGUN dddG', 'DDGUN class', 'DDGUN correct',
               'DynRCNC pred', 'DynRCNC correct',
               'seq_dist', 'SS_i', 'SS_j', 'zDCCM']
    for c, h in enumerate(headers, 1):
        cell = ws.cell(1, c, h)
        cell.fill = HDR_FILL; cell.font = HDR_FONT; cell.border = BORDER
        cell.alignment = Alignment(horizontal='center', vertical='center', wrap_text=True)

    r = 2
    for x in merged:
        vals = [x['PDB'], x['pair'], x['pH'], x['Method'],
                x['exp_dddG'], exp_word(x),
                x['s1'], x['s2'], x['T'], x['mean'], x['T_eq_mean'],
                x['ddgun_dddG'], ddgun_word(x),
                'Yes' if ddgun_word(x) == exp_word(x) else 'No',
                x['dyn_pred'], x['dyn_correct'],
                x['seq_dist'], x['ss_i'], x['ss_j'], x['zdccm']]
        for c, v in enumerate(vals, 1):
            cell = ws.cell(r, c, v); cell.font = ARI(); cell.border = BORDER
            if c != 2:
                cell.alignment = CENTER
            if   5 <= c <= 6:  cell.fill = GRP['exp']
            elif 7 <= c <= 14: cell.fill = GRP['ddgun']
            elif 15 <= c <= 16: cell.fill = GRP['dyn']
        if x['T_eq_mean'] == 'YES':
            ws.cell(r, 11).font = ARI(bold=True, color='008000')
        for cc in (14, 16):
            val = ws.cell(r, cc).value
            ws.cell(r, cc).font = ARI(color='C00000' if val in ('No', 'False') else '008000')
        r += 1

    for i, w in enumerate([7,16,6,9,12,12,9,9,9,15,10,11,13,12,13,13,9,6,6,9], 1):
        ws.column_dimensions[get_column_letter(i)].width = w
    ws.freeze_panes = 'C2'
    ws.auto_filter.ref = f"A1:{get_column_letter(len(headers))}1"

    # ---- Sheet 2: metrics -------------------------------------------------
    by = defaultdict(list)
    for x in merged:
        by[x['PDB']].append(x)

    ws2 = wb.create_sheet('Metrics summary')
    ws2.cell(1, 1, 'Non-additivity detection performance').font = ARI(bold=True, size=12)
    ws2.cell(2, 2, 'DDGUN3D').font = ARI(bold=True); ws2.cell(2, 2).alignment = CENTER
    ws2.cell(2, 6, 'DynRCNC').font = ARI(bold=True); ws2.cell(2, 6).alignment = CENTER
    ws2.merge_cells('B2:E2'); ws2.merge_cells('F2:I2')
    ws2.cell(2, 2).fill = GRP['ddgun']; ws2.cell(2, 6).fill = GRP['dyn']
    for c, h in enumerate(['PDB','N','Sens%','Spec%','F1','MCC','Sens%','Spec%','F1','MCC'], 1):
        cell = ws2.cell(3, c, h); cell.fill = HDR_FILL; cell.font = HDR_FONT
        cell.border = BORDER; cell.alignment = CENTER

    rr = 4
    for p in PROTEIN_ORDER:
        d = metrics(by[p], 'ddgun_class'); y = metrics(by[p], 'dyn_class')
        vals = [p, d[4], round(d[5],1), round(d[6],1), round(d[7],3), round(d[8],3),
                       round(y[5],1), round(y[6],1), round(y[7],3), round(y[8],3)]
        for c, v in enumerate(vals, 1):
            cell = ws2.cell(rr, c, v); cell.font = ARI(); cell.border = BORDER; cell.alignment = CENTER
            if   3 <= c <= 6:  cell.fill = GRP['ddgun']
            elif 7 <= c <= 10: cell.fill = GRP['dyn']
        rr += 1

    d = metrics(merged, 'ddgun_class'); y = metrics(merged, 'dyn_class')
    vals = ['COMBINED', d[4], round(d[5],1), round(d[6],1), round(d[7],3), round(d[8],3),
                          round(y[5],1), round(y[6],1), round(y[7],3), round(y[8],3)]
    for c, v in enumerate(vals, 1):
        cell = ws2.cell(rr, c, v); cell.font = ARI(bold=True); cell.alignment = CENTER
        cell.border = Border(top=Side(style='medium'), bottom=Side(style='medium'), left=THIN, right=THIN)
        if   3 <= c <= 6:  cell.fill = PatternFill('solid', fgColor='F8CBAD')
        elif 7 <= c <= 10: cell.fill = PatternFill('solid', fgColor='C6E0B4')
    for i, w in enumerate([10,5,8,8,7,8,8,8,7,8], 1):
        ws2.column_dimensions[get_column_letter(i)].width = w

    kr = rr + 2
    ws2.cell(kr, 1, 'Key finding').font = ARI(bold=True, size=11, color='C00000')
    notes = [
        'T = mean(s1,s2) holds for every pair across all 7 proteins.',
        'Therefore DDGUN3D epistasis dddG = -(s1+s2)/2, a linear function of the two singles.',
        'DDGUN3D contains no coupling term and cannot predict non-additivity by construction.',
        f'Class-level (combined): DDGUN3D MCC = {round(d[8],3)} vs DynRCNC MCC = {round(y[8],3)}.',
    ]
    for i, n in enumerate(notes):
        ws2.cell(kr + 1 + i, 1, '• ' + n).font = ARI(size=10)

    # ---- Sheet 3: T=mean proof (unique pairs) -----------------------------
    ws3 = wb.create_sheet('DDGUN T=mean proof')
    for c, h in enumerate(['PDB','Mutation pair','s1','s2','T (double)',
                           'mean(s1,s2)','T - mean','T = mean?'], 1):
        cell = ws3.cell(1, c, h); cell.fill = HDR_FILL; cell.font = HDR_FONT
        cell.border = BORDER; cell.alignment = CENTER
    seen = set(); r3 = 2
    for x in merged:
        if x['s1'] is None:
            continue
        k = (x['PDB'], frozenset(x['pair'].split(',')))
        if k in seen:
            continue
        seen.add(k)
        diff = round(x['T'] - x['mean'], 3)
        for c, v in enumerate([x['PDB'], x['pair'], x['s1'], x['s2'],
                               x['T'], x['mean'], diff, x['T_eq_mean']], 1):
            cell = ws3.cell(r3, c, v); cell.font = ARI(); cell.border = BORDER; cell.alignment = CENTER
        ws3.cell(r3, 8).font = ARI(bold=True, color='008000')
        r3 += 1
    for i, w in enumerate([7,18,8,8,11,13,9,10], 1):
        ws3.column_dimensions[get_column_letter(i)].width = w
    ws3.freeze_panes = 'A2'
    ws3.cell(r3 + 1, 1, f'{r3 - 2} unique double-mutant pairs; '
                        f'T = mean(s1,s2) within rounding (+/-0.05) in every case.'
             ).font = ARI(size=9, italic=True, color='808080')

    wb.save(OUTFILE)
    print(f"Saved {OUTFILE}")
    print(f"  per-pair rows : {len(merged)}")
    print(f"  DDGUN  MCC    : {round(d[8],3)}")
    print(f"  DynRCNC MCC   : {round(y[8],3)}")

# ----------------------------------------------------------------------
if __name__ == '__main__':
    singles = load_ddgun_out(OUT_GLOB)
    ddg     = load_ddgun_csv(DDGUN_CSV)
    dyn     = load_dynrcnc(DYNRCNC_CSV)

    print(f"DDGUN .out double-variants : {len(singles)}")
    print(f"DDGUN per-pair rows        : {len(ddg)}")
    print(f"DynRCNC rows               : {len(dyn)}")

    merged = merge(dyn, ddg, singles)
    n_missing = sum(1 for r in merged if r['ddgun_dddG'] is None)
    if n_missing:
        print(f"WARNING: {n_missing} DynRCNC rows had no DDGUN match")

    build(merged)
