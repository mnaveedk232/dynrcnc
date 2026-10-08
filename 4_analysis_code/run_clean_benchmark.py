#!/usr/bin/env python3
"""
Run DynRCNC on the CLEAN benchmark and produce the four result files
====================================================================
This runs the model exactly as the manuscript pipeline does, but on the
PBC-fixed 7-protein clean benchmark. It writes:

    clean_all_predictions.csv      per-pair predictions (265 rows)
    clean_summary.csv              per-protein + COMBINED confusion + metrics
    clean_lopocv.csv               leave-one-protein-out cross-validation
    clean_threshold_sensitivity.csv  MCC at |dddG| cutoffs 0.5-1.5

It calls dynrcnc.py's own main() so nothing about the method changes; it only
points the paths at the clean folder and renames the outputs with a clean_ prefix.

FOLDER LAYOUT expected (this script sits inside 00_clean_benchmark/):
    run_clean_benchmark.py     <- this file
    dynrcnc.py                 <- the model
    benchmark_265pairs.ddg
    1stn_md/ 1stn.N 1stn.E md_first_frame.gro md_reduced.xtc
    1bni_md/ 2lzm_md/ 1pga_md/ 1csp_md/ 2rn2_md/ 2ci2_md/   (same pattern)

NOTE on filenames: dynrcnc.py expects <pdb>_md_first_frame.gro and
<pdb>_md_reduced.xtc inside each <pdb>_md/ folder. The clean folder stores them
as md_first_frame.gro and md_reduced.xtc. This script makes the expected names
available automatically (via symlinks) without touching your originals, so you
do not have to rename anything.

Run:
    python3 run_clean_benchmark.py
Needs: numpy pandas networkx MDAnalysis mdtraj scipy
"""

import os, sys, shutil, importlib.util

HERE = os.path.dirname(os.path.abspath(__file__))
PROTEINS = ['1stn', '1bni', '2lzm', '1pga', '1csp', '2rn2', '2ci2']
OUT = os.path.join(HERE, 'clean_out')

# ------------------------------------------------------------------
# 1. Make the filenames dynrcnc.py expects, without renaming originals.
#    md_first_frame.gro  ->  <pdb>_md_first_frame.gro   (symlink)
#    md_reduced.xtc      ->  <pdb>_md_reduced.xtc       (symlink)
# ------------------------------------------------------------------
def ensure_expected_names():
    for p in PROTEINS:
        d = os.path.join(HERE, f'{p}_md')
        pairs = [('md_first_frame.gro', f'{p}_md_first_frame.gro'),
                 ('md_reduced.xtc',     f'{p}_md_reduced.xtc')]
        for src, dst in pairs:
            src_path = os.path.join(d, src)
            dst_path = os.path.join(d, dst)
            if os.path.exists(src_path) and not os.path.exists(dst_path):
                try:
                    os.symlink(src, dst_path)          # relative symlink inside the folder
                except OSError:
                    shutil.copy(src_path, dst_path)    # fallback: copy if symlink not allowed

# ------------------------------------------------------------------
# 2. Import dynrcnc.py and run its main() with our paths.
# ------------------------------------------------------------------
def main():
    ensure_expected_names()
    os.makedirs(OUT, exist_ok=True)

    # dynrcnc.py reads --base/--ddg/--out from argv, so set them here.
    sys.argv = ['dynrcnc.py',
                '--base', HERE,
                '--ddg', os.path.join(HERE, 'benchmark_265pairs.ddg'),
                '--out', OUT]
    spec = importlib.util.spec_from_file_location('dynrcnc', os.path.join(HERE, 'dynrcnc.py'))
    dm = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(dm)
    dm.main()

    # 3. Copy outputs up with a clean_ prefix, next to this script.
    for name in ['all_predictions', 'summary', 'lopocv', 'threshold_sensitivity']:
        src = os.path.join(OUT, f'{name}.csv')
        if os.path.exists(src):
            shutil.copy(src, os.path.join(HERE, f'clean_{name}.csv'))
    print("\nWrote: clean_all_predictions.csv, clean_summary.csv, "
          "clean_lopocv.csv, clean_threshold_sensitivity.csv")


if __name__ == '__main__':
    main()
