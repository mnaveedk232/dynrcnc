Supporting Data File S1 - per-pair predictions of DynRCNC, FoldX, Rosetta cartesian_ddg and DDGun3D
Manuscript: Context-Dependent Dynamic Coupling in Residue Contact Networks Underlies Thermodynamic Non-Additivity in Protein Double Mutants

Contents
  FoldX_vs_Benchmark_vs_DynRCNC_265pairs.xlsx          all 265 benchmark pairs (Table 6)
  DDGun3D_vs_Benchmark_vs_DynRCNC_265pairs.xlsx        all 265 benchmark pairs (Section 3.7, Table 6)
  Rosetta_vs_Benchmark_vs_DynRCNC_4proteins.xlsx       72 pairs of 1PGA, 1CSP, 2CI2 and 1BNI (Section 3.7)

Each workbook has a README sheet (definitions), a Per-pair sheet, and a Summary sheet whose numbers are those reported in the manuscript:
  FoldX    sensitivity 24.1%, specificity 89.1%, MCC 0.155
  DDGun3D  sensitivity 25.9%, specificity 85.3%, MCC 0.121
  Rosetta  best cutoff 0.9 REU: sensitivity 96.4%, specificity 25.0%, MCC 0.28; DynRCNC on the same 72 pairs MCC 0.513
  DynRCNC  265 pairs: sensitivity 74.1%, specificity 87.2%, MCC 0.568
Positive class: non-additive (|ddd G| >= 1.0 kcal/mol).
