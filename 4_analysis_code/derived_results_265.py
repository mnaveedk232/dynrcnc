#!/usr/bin/env python3
"""Regenerates the replicate and per-criterion result files of the 265-pair analysis from the primary run.

Replicate 1 is the primary trajectory (results/clean_all_predictions.csv); replicates 2 and 3 are read from
results/replicate_zdccm_per_pair.csv. Also writes results/per_criterion_precision.csv.
Files of the older run (replicate 1 before the periodic-boundary correction of 1PGA and 2CI2) are moved to a backup folder.
Usage:  python3 derived_results_265.py [path to the github folder, default ~/Desktop/github2]"""
import sys, os, math, shutil
import numpy as np, pandas as pd
R = os.path.expanduser(sys.argv[1] if len(sys.argv) > 1 else '~/Desktop/github2'); B = os.path.expanduser('~/Desktop/github2_stale_backup'); os.makedirs(B, exist_ok=True)
P = lambda p: os.path.join(R, p)
c = pd.read_csv(P('results/clean_all_predictions.csv')); rp = pd.read_csv(P('results/replicate_zdccm_per_pair.csv'))
assert len(c) == 265 and len(rp) == 265
def stash(rel):
    if os.path.exists(P(rel)): shutil.move(P(rel), os.path.join(B, os.path.basename(rel))); print('moved to backup:', rel)
# ---- replicates
cm = c.copy(); cm['mut'] = cm.mutation.str.replace(' ', ''); cm['k'] = cm.groupby(['PDB', 'mut', 'dddG']).cumcount()
rp = rp.copy(); rp['mut'] = rp.mut.str.replace(' ', ''); rp['k'] = rp.groupby(['PDB', 'mut', 'dddG']).cumcount()
j = rp.merge(cm[['PDB', 'mut', 'dddG', 'k', 'zdccm', 'pred']], on=['PDB', 'mut', 'dddG', 'k'], how='left'); assert j.zdccm.notna().all() and len(j) == 265
shutil.copy(P('results/replicate_zdccm_per_pair.csv'), os.path.join(B, 'replicate_zdccm_per_pair_old.csv'))
out = j[['PDB', 'mut', 'exp', 'dddG']].copy(); out['z_rep1'] = j.zdccm; out['pred_rep1'] = j.pred; out['z_rep2'] = j.z_rep2; out['pred_rep2'] = j.pred_rep2; out['z_rep3'] = j.z_rep3; out['pred_rep3'] = j.pred_rep3
Z = out[['z_rep1', 'z_rep2', 'z_rep3']]; out['z_mean'] = Z.mean(axis=1); out['z_sd'] = Z.std(axis=1, ddof=1)
out.to_csv(P('results/replicate_zdccm_per_pair.csv'), index=False)
def stats(pr):
    y = (out.exp == 'Non-additive').values; q = (pr == 'Non-additive').values; tp, tn, fp, fn = int((y & q).sum()), int((~y & ~q).sum()), int((~y & q).sum()), int((y & ~q).sum())
    d = math.sqrt(float((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))); return tp, tn, fp, fn, 100 * tp / (tp + fn), 100 * tn / (tn + fp), (tp * tn - fp * fn) / d
ms = []; lines = ['=' * 60, '  REPLICATE Z-DCCM CONVERGENCE ANALYSIS', '  265 benchmark pairs x 3 replicates (replicate 1 = primary trajectory)', '=' * 60, '', '1. Z-DCCM variability across replicates',
    f'   Mean per-pair SD of |Z-DCCM| : {out.z_sd.mean():.3f}', f'   Median per-pair SD           : {out.z_sd.median():.3f}', f'   Mean |Z-DCCM|                : {Z.values.mean():.3f}', '   Inter-replicate Pearson r:']
for a, b in [('z_rep1', 'z_rep2'), ('z_rep1', 'z_rep3'), ('z_rep2', 'z_rep3')]: lines.append(f'     {a} vs {b}: {out[a].corr(out[b]):.3f}')
same = (out.pred_rep1 == out.pred_rep2) & (out.pred_rep2 == out.pred_rep3)
lines += ['', '2. Classification stability', f'   Identical prediction in all replicates: {int(same.sum())}/265 ({100 * same.mean():.1f}%)', f'   Pairs whose prediction changes        : {int((~same).sum())}', '', '3. Combined MCC per replicate']
for k in (1, 2, 3):
    tp, tn, fp, fn, se, sp, m = stats(out[f'pred_rep{k}']); ms.append(m); lines.append(f'   rep{k}: TP={tp} TN={tn} FP={fp} FN={fn}  Sens={se:.1f}% Spec={sp:.1f}%  MCC={m:.3f}')
lines += ['', f'   MCC mean = {np.mean(ms):.3f} +/- {np.std(ms, ddof=1):.3f}  (range {min(ms):.3f}-{max(ms):.3f})']
open(P('results/replicate_convergence_summary.txt'), 'w').write('\n'.join(lines) + '\n'); print('\n'.join(lines))
# ---- per-criterion precision
q = c[c.pred == 'Non-additive'].copy(); q['crit'] = q.mechanism.str.extract(r'_(C\d)_')[0]
g = q.groupby('crit').agg(fired=('exp', 'size'), TP=('exp', lambda x: int((x == 'Non-additive').sum()))); g['FP'] = g.fired - g.TP; g['precision'] = (g.TP / g.fired).round(3)
if os.path.exists(P('results/per_criterion_precision.csv')): shutil.copy(P('results/per_criterion_precision.csv'), os.path.join(B, 'per_criterion_precision_old.csv'))
g.reset_index().rename(columns={'crit': 'criterion'})[['criterion', 'fired', 'TP', 'FP', 'precision']].to_csv(P('results/per_criterion_precision.csv'), index=False)
print('\nper-criterion precision:\n' + g.reset_index().to_string(index=False))
# ---- files of the older run
for rel in ('results/cutoff_sensitivity.csv', 'results/cascade_multi_criterion.csv', 'code/cascade_precision_cutoff.py', 'code/replicate_convergence.py'): stash(rel)
tx = P('results/cascade_precision_cutoff_summary.txt')
if os.path.exists(tx) and any(k in open(tx).read() for k in ('0.554', '0.467', 'TN=182')): stash('results/cascade_precision_cutoff_summary.txt')
os.makedirs(P('code'), exist_ok=True); shutil.copy(os.path.abspath(__file__), P('code/derived_results_265.py'))
print('\nbackup of the replaced files:', B)
