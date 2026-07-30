# DynRCNC

Predicts epistasis in protein double mutations, meaning whether the two substitutions act independently (additive) or interfere with each other (non-additive). DynRCNC reads a residue contact network and a molecular dynamics trajectory of the wild type, applies six ordered decision rules, and reports the rule that decided each pair.

On a benchmark of 271 experimentally measured pairs from eight proteins it reaches an MCC of 0.557.

## Install

Python 3.9 or newer.

```bash
git clone https://github.com/mnaveedk232/dynrcnc.git
cd dynrcnc
pip install -r requirements.txt
```

Dependencies are NumPy, pandas, SciPy, NetworkX, MDAnalysis and MDTraj, all pinned in `requirements.txt`.

## Quick start

```bash
python3 dynrcnc.py --base /path/to/data
```

```
--base   directory holding the *_md/ folders. Default ./data, or set DYNRCNC_BASE.
--ddg    benchmark file. Default ./benchmark_271pairs.ddg.
--out    output directory. Default ./output.
```

With the protein folders in `./data` the flags can be left off entirely.

```bash
python3 dynrcnc.py
```

Eight proteins take a few minutes.

## Input

One folder per protein, four files in each:

```
data/
  1stn_md/
    1stn.N                      RING nodes
    1stn.E                      RING edges
    1stn_md_first_frame.gro     topology
    1stn_md_reduced.xtc         trajectory
  1bni_md/  2lzm_md/  1pga_md/  1csp_md/  2rn2_md/  2ci2_md/  1qjp_md/
```

Labels are in `benchmark_271pairs.ddg`, tab separated, one row per pair:

```
Mutation_Double  PDB   pH   Method   DDG_Double  DDG_Single1  DDG_Single2  dddG
I6T,T53F         1PGA  5.2  Thermal  1.5         1.76         4.1          -4.36
```

A pair is non-additive when `|dddG|` is 1.0 kcal/mol or more.

## Contact network

Structures come from the RCSB PDB. Remove waters, ligands and alternate locations, and keep a single chain. Networks are built with RING, either on the web server or from a local install. URLs are listed at the end of this file.

```bash
ring -i 1stn.pdb --out_dir ring_out
mv ring_out/1stn.pdb_ringNodes 1stn_md/1stn.N
mv ring_out/1stn.pdb_ringEdges 1stn_md/1stn.E
```

The benchmark networks use RING 4.0 defaults.

## Trajectory

Only the wild type is simulated. All eight benchmark trajectories are 200 ns. Six were written every 100 ps for 2001 frames, and 1PGA and 1CSP every 40 ps for 5001 frames. Simulation conditions are in the Methods section of the paper.

```bash
gmx trjconv -s md.tpr -f md.xtc -o 1stn_md_first_frame.gro -dump 0
gmx trjconv -s md.tpr -f md.xtc -o 1stn_md_reduced.xtc -dt 100 -pbc mol -center
```

## Output

Per protein, the confusion matrix, then every pair with its experimental class, predicted class, agreement and the rule that fired.

```
RESULTS: 1PGA
  TP=7 TN=1 FP=0 FN=0 N=8
  Sensitivity : 100.0%  (7/7)
  Specificity : 100.0%  (1/1)
  MCC         : 1.000

  ── Non-additive pairs (7) ──
  mutation  pH  Method  dddG          exp         pred  correct               mechanism
  I6T,T53F 5.2 Thermal -4.36 Non-additive Non-additive     True R10_C5_hub_perturbation
```

I6T/T53F is non-additive because residue 6 is a hub in the contact network. That is rule R10, under criterion C5.

Four files are written to the output directory.

```
all_predictions.csv          one row per pair with the mechanism that fired
summary.csv                  metrics per protein and the combined row
lopocv.csv                   leave-one-protein-out folds
threshold_sensitivity.csv    MCC against each tunable threshold
```

Bootstrap confidence intervals, a permutation test, a McNemar test and the generalisation gap are printed at the end of the run. The McNemar test needs the per-pair predictions of the virtual-edge baseline, which are not distributed with the code, so by default it compares against the static clique baseline instead. The run states which baseline it used.

## Criteria

Six criteria are tested in order. The first to fire decides the pair. A pair that fires nothing is additive.

```
C1   the two residues share a 3-clique community in the static network
C2   direct contact plus dynamic coupling
C3   dynamic packing
C4   sequential backbone neighbours with correlated motion
C5   network topology, for example a shared hub
C6   isolated charged residues with electrostatic coupling
```

C1 uses the structure only. C2 to C5 use the trajectory. C6 is physicochemical. Rules map to criteria in `CRITERION_OF_RULE` at the top of `dynrcnc.py`.

## Benchmark

271 pairs from 1STN, 1BNI, 2LZM, 1PGA, 1CSP, 2RN2, 2CI2 and 1QJP.

```
TP = 42   TN = 185   FP = 29   FN = 15
Sensitivity  73.7%
Specificity  86.4%
MCC          0.557    95% CI 0.441 to 0.663
```

Leave-one-protein-out cross-validation, thresholds re-optimised on the training proteins of each fold, gives a weighted MCC of 0.417 and a generalisation gap of 0.141.

Other methods on the same 271 pairs:

```
FoldX                                              0.162
best virtual edge, of 135 threshold combinations   0.145
static contact network, no dynamics                0.079
random forest, cross-validated                     0.029
gradient boosting, cross-validated                 0.032
```

The random forest reaches 0.794 in sample.

## Notes

`fps` in the `PROTEINS` dictionary is the write interval in ps. It works with `SKIP_NS` to drop equilibration frames. `SKIP_NS` is 0 in this version, so nothing is dropped.

Changing RING's hydrogen bond or van der Waals thresholds changes the clique communities, and therefore C1.

Trajectories and contact networks are too large for GitHub. They are available from the authors on request.

## Contents

```
dynrcnc.py                                     the method
benchmark_271pairs.ddg                         benchmark labels
benchmark_271pairs_dynrcnc_predictions.xlsx    all 271 pairs with predictions
requirements.txt                               dependencies
LICENSE                                        MIT
```

The workbook has two sheets. Predictions holds one row per pair with the measured ΔΔΔG, |Z-DCCM|, observed class, predicted class, deciding criterion and outcome. Legend explains the columns.

## Links

- RCSB PDB: https://www.rcsb.org
- RING server: https://ring.biocomputingup.it
- RING download: https://biocomputingup.it/services/download
- GROMACS: https://www.gromacs.org

## Citation

Naveed, M., Ming, D. *DynRCNC: A Hierarchical Framework Integrating Molecular Dynamics with Residue Contact Network Topology to Predict Mutation Epistasis*

## License

MIT. See the LICENSE file.
