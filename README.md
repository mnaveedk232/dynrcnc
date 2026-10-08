# DynRCNC

DynRCNC predicts whether the effects of two mutations in a protein add up (additive) or interfere with each other (non-additive). It needs two inputs for the wild-type protein: a residue contact network built with RING, and a molecular dynamics (MD) trajectory. It then applies six decision criteria, C1 to C6, one after the other. The first criterion that fires decides the pair, and the output states which one it was.

On a benchmark of 265 measured double-mutant pairs from seven soluble proteins, DynRCNC reaches an MCC of 0.568. This value was obtained on the same data on which most thresholds were set, so it is an optimistic figure. The cross-validated value is lower (see "Results at a glance").

This repository holds everything behind the numbers in the paper:

Naveed M, Zhang J, Sajid AQ, Fatima A, Ming D. *Context-Dependent Dynamic Coupling in Residue Contact Networks Underlies Thermodynamic Non-Additivity in Protein Double Mutants.* Biopolymers (accepted).

## What is in the repository

The folders are numbered in the order of the Data availability statement of the paper.

```
1_benchmark_dataset/            The benchmark
    benchmark_265pairs.ddg                       265 pairs with the measured ddG values
    benchmark_265pairs_dynrcnc_predictions.xlsx  the same 265 pairs with the DynRCNC prediction
    external_final.ddg                           47 pairs of four further proteins (external check)

2_RING_contact_files/           RING 4.0 contact networks, one folder per protein
    1STN  1BNI  2LZM  1PGA  1CSP  2RN2  2CI2     each holds <pdb>.N (nodes) and <pdb>.E (edges)

3_residue_level_DCCM_matrices/  Dynamic cross-correlation matrices, one pair of files per protein
    <PDB>_DCCM.csv                               DCCM of the C-alpha atoms
    <PDB>_Z-DCCM.csv                             Z-DCCM, the values the criteria use

4_analysis_code/                All analysis code
5_results/                      Result tables of the paper
6_supporting_data/              Supporting Data File S1 (FoldX, DDGun3D and Rosetta predictions per pair)
7_figure_scripts/               Scripts for the figures
requirements.txt                Python packages
LICENSE                         MIT license
```

**About the matrix files.** Every matrix is square, with one row and one column per residue, labelled with residue name and number (for example `ILE6`). DCCM is computed on the C-alpha atoms over the whole trajectory after superposition on the first frame. Z-DCCM is the DCCM value of a residue pair, standardized against all residue pairs with the same sequence separation. The Z-DCCM value of every benchmark pair is also listed in `5_results/clean_all_predictions.csv`, column `zdccm` (as an absolute value), and the two agree.

**Code in `4_analysis_code/`.**

| Script | What it does |
|---|---|
| `dynrcnc.py` | The method itself (all rules and thresholds) |
| `run_clean_benchmark.py`, `verify_clean.py` | Run the full benchmark and check it |
| `ablation_table7.py`, `ablation_dynamic_gates.py` | Ablation tests |
| `fallthrough_test.py`, `order_test.py`, `skip50_test.py` | Tests of criterion order and fall-through |
| `static_rcnc_265.py` | Static RCNC baseline (no dynamics) |
| `ve_clean.py`, `ve_merge_mcnemar_265.py` | Virtual-edge baseline (135 threshold combinations) and McNemar test |
| `ml_baseline_265.py`, `ml_tuned_265.py` | Machine-learning baselines |
| `stats_265.py`, `mcc_diff_bootstrap_265.py`, `within_protein_perm.py` | Confidence intervals and permutation tests |
| `derived_results_265.py` | Replicate analysis and per-criterion precision |
| `dynrcnc_external.py` | External check on four further proteins |
| `build_ddgun_comparison.py`, `build_rosetta_comparison.py` | Tables with the DDGun3D and Rosetta comparison |
| `dynrcnc_predictions_265.py` | Writes the predictions workbook of folder 1 |
| `export_ring_dccm.py` | Checks the RING files and rebuilds the matrices of folders 2 and 3 |

## Results at a glance

Positive class: non-additive, meaning |ddG(double) minus ddG(single 1) minus ddG(single 2)| is 1.0 kcal/mol or more.

```
265 pairs: 1STN 120, 2LZM 57, 2CI2 29, 1BNI 25, 2RN2 16, 1CSP 10, 1PGA 8
54 non-additive, 211 additive

TP = 40   TN = 184   FP = 27   FN = 14
Sensitivity  74.1 %
Specificity  87.2 %
MCC          0.568    (95% CI 0.443 to 0.687; 0.450 to 0.680 when whole proteins are resampled)
```

Leave-one-protein-out cross-validation, with three thresholds re-fitted on the other six proteins, gives a pair-weighted MCC of 0.422. The remaining thresholds were set on the full benchmark, so read 0.422 as an optimistic bound for a protein the method has not seen. Three independent 200 ns trajectories per protein give an MCC of 0.46 ± 0.10.

Other methods on the same 265 pairs:

```
FoldX                                              0.155
DDGun3D                                            0.121
Best virtual edge (of 135 threshold combinations)  0.123
Static contact network, no dynamics                0.083
Random forest, cross-validated                     0.085
Gradient boosting, cross-validated                 0.030
Logistic regression, cross-validated               0.020
```

The random forest reaches 0.786 when it is tested on its own training data, which shows how strongly it overfits. Rosetta cartesian_ddg was run on four proteins (72 pairs): MCC 0.28 at its best cutoff, against 0.513 for DynRCNC on the same pairs.

As a first test on proteins outside the benchmark, DynRCNC was applied to 1JQZ, 1C9O, 1MJC and 1OIA (47 pairs, `1_benchmark_dataset/external_final.ddg`). The MCC is 0.23 (95% CI -0.09 to 0.47). With only seven non-additive pairs this cannot be told apart from chance.

## Install

Python 3.9 or newer.

```
git clone https://github.com/mnaveedk232/dynrcnc.git
cd dynrcnc
pip install -r requirements.txt
```

`dynrcnc.py` needs NumPy, pandas, SciPy, NetworkX, MDAnalysis and MDTraj. Some of the other scripts also use scikit-learn, matplotlib, openpyxl and Pillow, so install those if a script asks for them. The commands below were written for Linux or macOS. On Windows, run the Python commands as they are (use `python` instead of `python3` if needed) and copy files by hand instead of using the shell loop.

## Check the numbers (no trajectory needed)

The per-pair predictions are in `5_results/clean_all_predictions.csv`. This command rebuilds the confusion matrix and the MCC from that file:

```
python3 -c "
import pandas as pd, math
d = pd.read_csv('5_results/clean_all_predictions.csv')
y = d.exp == 'Non-additive'; p = d.pred == 'Non-additive'
tp, tn, fp, fn = (y & p).sum(), (~y & ~p).sum(), (~y & p).sum(), (y & ~p).sum()
print(len(d), tp, tn, fp, fn, round((tp*tn - fp*fn) / math.sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn)), 3))
"
```

The expected output is `265 40 184 27 14 0.568`. The file also contains, for every pair, the criterion that decided it (column `mechanism`).

## Run DynRCNC on the benchmark

The method needs one folder per protein with four files: RING nodes, RING edges, a topology file (`.gro`) and the trajectory (`.xtc`).

```
data/
  1stn_md/
    1stn.N                      RING nodes
    1stn.E                      RING edges
    1stn_md_first_frame.gro     topology (first frame)
    1stn_md_reduced.xtc         trajectory
  1bni_md/  2lzm_md/  1pga_md/  1csp_md/  2rn2_md/  2ci2_md/   (same four files each)
```

1. The RING files are in `2_RING_contact_files/`. Copy them into this layout (Linux or macOS, run inside the repository folder):

```
for p in 1STN 1BNI 2LZM 1PGA 1CSP 2RN2 2CI2; do l=$(echo $p | tr A-Z a-z); mkdir -p data/${l}_md; cp 2_RING_contact_files/$p/* data/${l}_md/; done
```

2. Add the `.gro` and `.xtc` files of each protein to its folder. The trajectories are too large for GitHub (see "What is not included").

3. Run:

```
python3 4_analysis_code/dynrcnc.py --base data --ddg 1_benchmark_dataset/benchmark_265pairs.ddg --out output
```

`--base` is the folder with the `*_md` folders (default `./data`), `--ddg` is the benchmark file, `--out` is the output folder (default `./output`). The run takes a few minutes and writes four files:

```
all_predictions.csv          one row per pair, with the criterion that decided it
summary.csv                  confusion matrix and MCC per protein, and combined
lopocv.csv                   leave-one-protein-out folds
threshold_sensitivity.csv    MCC for different class boundaries (0.5 to 1.5 kcal/mol)
```

With the original trajectories these four files are identical to `clean_all_predictions.csv`, `clean_summary.csv`, `clean_lopocv.csv` and `clean_threshold_sensitivity.csv` in `5_results/`.

## Check and rebuild the RING and matrix folders

`4_analysis_code/export_ring_dccm.py` recomputes DCCM and Z-DCCM from the trajectories with the same function that `dynrcnc.py` uses. For every protein it checks that the RING files belong to it (same residue numbers and names as the topology) and that |Z-DCCM| of every benchmark pair equals the value in `5_results/clean_all_predictions.csv`. Only proteins that pass are copied to the output folder.

```
python3 4_analysis_code/export_ring_dccm.py --base data --dynrcnc 4_analysis_code/dynrcnc.py --pred 5_results/clean_all_predictions.csv --out checked_output
```

All seven proteins pass with the files in this repository.

## Figures

The scripts in `7_figure_scripts/` read their input from `5_results/` and the other folders of this repository. Run them from inside `7_figure_scripts/`: first the `*_data.py` script (it writes a small csv), then the `plot_*.py` script of the same figure.

```
cd 7_figure_scripts
python3 fig2_data.py     && python3 plot_fig2.py      # Figure 2
python3 fig4_data.py     && python3 plot_fig4.py      # Figure 4 and, with plot_figS3.py, Figure S3
python3 fig5_data.py     && python3 plot_fig5.py      # Figure 5
python3 make_figure1.py                               # Figure 1
```

Figure S1 (RMSD) and Figure S2 (threshold sensitivity) need the trajectories. Put them into `data/` as described above (or point `MD_BASE` to your folder), then run `figS1_data.py` and `plot_figS1.py`, and `figS2_data.py` and `plot_figS2.py`.

`select_figure3_panels.py` lists, for each criterion, all correctly predicted non-additive pairs. The six pairs drawn in Figure 3 are all in that list:

```
Panel  Criterion  Protein  Residues  Double mutant
A      C1         1STN     79, 118   G79S,N118D
B      C2         1BNI     13, 17    Y13A,Y17A
C      C3         2LZM     121, 129  L121A,A129M
D      C4         1BNI     16, 17    T16S,Y17A
E      C5         1PGA     6, 53     I6T,T53F
F      C6         2LZM     16, 135   K16E,K135E
```

The script itself suggests the pair with the largest |ddG| in each criterion, which is not always the pair shown in the paper.

## The six criteria

The criteria are tested in the order C1 to C6. A pair that fires none of them is called additive.

```
C1   both residues belong to the same 3-clique community of the contact network
C2   the residues touch directly in the contact network, and their motion is coupled
C3   dynamic packing: coupled motion in flexible coil or isolated helical regions
C4   sequential backbone neighbors with correlated motion
C5   network topology, for example a hub residue or a shared neighbor
C6   two isolated charged residues
```

C1 uses the structure only, C2 to C5 also use the trajectory, and C6 depends on the residue types. The rules behind each criterion and their thresholds are at the top of `4_analysis_code/dynrcnc.py` (`CRITERION_OF_RULE` and the constants beside it).

## Input formats

The benchmark file is tab separated, one row per measurement:

```
Mutation_Double  PDB   pH   Method   DDG_Double  DDG_Single1  DDG_Single2  dddG
I6T,T53F         1PGA  5.2  Thermal  1.5         1.76         4.1          -4.36
```

A pair is non-additive when |dddG| is 1.0 kcal/mol or more. If the same pair was measured under several conditions (pH, method), each measurement is one row.

The benchmark was compiled from the ProThermDB database and the primary literature cited in the paper.

**RING files.** Structures were taken from the RCSB PDB, with waters, ligands and alternative locations removed and a single chain kept. The networks were built with RING 4.0 using strict distance thresholds and the multiple-interaction policy (see the Methods section of the paper). To make a network for another protein:

```
ring -i 1stn.pdb --out_dir ring_out
mv ring_out/1stn.pdb_ringNodes 1stn_md/1stn.N
mv ring_out/1stn.pdb_ringEdges 1stn_md/1stn.E
```

**Trajectories.** Only the wild-type protein is simulated (200 ns, GROMACS, AMBER99SB-ILDN, TIP3P water; details are in the Methods section). For the analysis, five of the proteins were subsampled to 100 ps per frame (2001 frames) and 1PGA and 1CSP to 40 ps per frame (5001 frames):

```
gmx trjconv -s md.tpr -f md.xtc -o 1stn_md_first_frame.gro -dump 0
gmx trjconv -s md.tpr -f md.xtc -o 1stn_md_reduced.xtc -dt 100 -pbc mol -center
```

For 1PGA and 1CSP use `-dt 40`.

## What is not included

- **The MD trajectories** (`.xtc`) and topology files are several megabytes each and are not stored here. They are available from the authors on request. The Z-DCCM matrices in folder 3 are computed from them, so the criteria can be checked without the trajectories.
- **Replicate trajectories.** The paper reports three independent 200 ns runs per protein. This repository uses the first run for all results. The per-pair Z-DCCM values of all three runs are in `5_results/replicate_zdccm_per_pair.csv`.
- Most scripts in `4_analysis_code/` contain the folder layout of the authors' computer as default paths. If a script cannot find a file, change the path constants at the top of the script.

## Technical notes

`fps` in the `PROTEINS` dictionary of `dynrcnc.py` is the time between saved frames in ps. It works together with `SKIP_NS`, which drops the first nanoseconds of the trajectory as equilibration. In this version `SKIP_NS` is 0, so all frames are used.

Changing the hydrogen-bond or van der Waals thresholds of RING changes the clique communities and therefore criterion C1.

## Links

- RCSB PDB: https://www.rcsb.org
- RING server: https://ring.biocomputingup.it
- RING download: https://biocomputingup.it/services/download
- GROMACS: https://www.gromacs.org

## Contact

Corresponding author: Dengming Ming, College of Biotechnology and Pharmaceutical Engineering, Nanjing Tech University.

## License

MIT. See the file `LICENSE`.
