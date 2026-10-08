#!/usr/bin/env python3
"""Table 7 (criterion / tier ablation) on the CLEAN 265-pair benchmark.
A disabled component's predictions become additive (no fall-through), as in the
manuscript. Works directly on clean_all_predictions.csv, no MD needed.
Run inside 00_clean_benchmark/:   python3 ablation_table7.py"""
import pandas as pd, numpy as np
CRIT = {'R1':'C1','R2':'C1','R3':'C1','R4':'C2','R5':'C3','R6':'C3','R7':'C3','R8':'C3',
        'R9':'C4','R10':'C5','R11':'C5','R12':'C5','R13':'C5','R14':'C6'}
df = pd.read_csv('clean_all_predictions.csv')
df['crit'] = df.mechanism.map(lambda m: CRIT.get(str(m).split('_')[0], 'none'))
y = df.exp == 'Non-additive'
def met(p):
    q = p == 'Non-additive'
    tp=(y&q).sum(); tn=(~y&~q).sum(); fp=(~y&q).sum(); fn=(y&~q).sum()
    d = np.sqrt(float((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn)))
    return (tp*tn-fp*fn)/d if d else 0, tp, tn, fp, fn, 100*tp/(tp+fn), 100*tn/(tn+fp)
base = met(df.pred)[0]
rows = [('None (full model)',[]), ('Tier 2 (all MD-derived)',['C2','C3','C4','C5']),
        ('Tier 1 (static topology)',['C1']), ('C5 network topology',['C5']),
        ('C3 dynamic packing',['C3']), ('C2 direct contact',['C2']),
        ('C1 clique community',['C1']), ('C4 sequential backbone',['C4']),
        ('C6 electrostatic coupling',['C6'])]
out = []
for n, cs in rows:
    m = met(df.pred.where(~df.crit.isin(cs), 'Additive'))
    out.append({'Component removed':n,'MCC':round(m[0],3),'dMCC':round(m[0]-base,3),'TP':m[1],
                'TN':m[2],'FP':m[3],'FN':m[4],'Sens':round(m[5],1),'Spec':round(m[6],1)})
r = pd.DataFrame(out); print(f'Pairs: {len(df)}'); print(r.to_string(index=False))
r.to_csv('ablation_table7.csv', index=False)
