#!/usr/bin/env python3
"""
DynRCNC
=======
Predicts double-mutation thermodynamic non-additivity in proteins.

Method: Static RCNC network plus MD-derived dynamic coupling (DCCM, RMSF),
        organised as a three-tier hierarchical framework of six decision
        criteria (C1-C6). Each criterion is implemented by one or more
        internal rules; see CRITERION_OF_RULE below for the mapping.

Benchmark: 271 experimentally characterised double-mutation pairs across
           eight proteins.


Usage:
    python3 dynrcnc.py                     # run with default paths below
    python3 dynrcnc.py --base /path        # override the data directory
"""

import os, warnings, time
import numpy as np
import pandas as pd
import networkx as nx
import MDAnalysis as mda
import mdtraj as md
from networkx.algorithms.community import k_clique_communities
from collections import defaultdict, Counter
from scipy.stats import chi2 as chi2_dist
warnings.filterwarnings('ignore')

# ─────────────────────────────────────────────
#  CONFIGURATION
#  Data directory can be overridden without editing the file:
#    python3 dynrcnc.py --base /path/to/data
#    DYNRCNC_BASE=/path/to/data python3 dynrcnc.py
#  The directory must contain one <pdb>_md/ folder per protein.
# ─────────────────────────────────────────────
import argparse as _argparse

def _resolve_base():
    ap = _argparse.ArgumentParser(add_help=False)
    ap.add_argument('--base', default=None)
    ap.add_argument('--ddg', default=None)
    ap.add_argument('--out', default=None)
    args, _ = ap.parse_known_args()
    base = args.base or os.environ.get('DYNRCNC_BASE') or './data'
    return os.path.expanduser(base), args.ddg, args.out

BASE, _DDG_OVERRIDE, _OUT_OVERRIDE = _resolve_base()

PROTEINS = {
    '1STN': {'node_file': f'{BASE}/1stn_md/1stn.N', 'edge_file': f'{BASE}/1stn_md/1stn.E',
             'gro': f'{BASE}/1stn_md/1stn_md_first_frame.gro', 'xtc': f'{BASE}/1stn_md/1stn_md_reduced.xtc', 'fps': 100},
    '1BNI': {'node_file': f'{BASE}/1bni_md/1bni.N', 'edge_file': f'{BASE}/1bni_md/1bni.E',
             'gro': f'{BASE}/1bni_md/1bni_md_first_frame.gro', 'xtc': f'{BASE}/1bni_md/1bni_md_reduced.xtc', 'fps': 100},
    '2LZM': {'node_file': f'{BASE}/2lzm_md/2lzm.N', 'edge_file': f'{BASE}/2lzm_md/2lzm.E',
             'gro': f'{BASE}/2lzm_md/2lzm_md_first_frame.gro', 'xtc': f'{BASE}/2lzm_md/2lzm_md_reduced.xtc', 'fps': 100},
    '1PGA': {'node_file': f'{BASE}/1pga_md/1pga.N', 'edge_file': f'{BASE}/1pga_md/1pga.E',
             'gro': f'{BASE}/1pga_md/1pga_md_first_frame.gro', 'xtc': f'{BASE}/1pga_md/1pga_md_reduced.xtc', 'fps': 100},
    # 1CSP trajectory is sampled at 40 ps/frame (200 ns total, 5001 frames),
    # unlike the 100 ps/frame of the other proteins; fps below reflects this.
    '1CSP': {'node_file': f'{BASE}/1csp_md/1csp.N', 'edge_file': f'{BASE}/1csp_md/1csp.E',
             'gro': f'{BASE}/1csp_md/1csp_md_first_frame.gro', 'xtc': f'{BASE}/1csp_md/1csp_md_reduced.xtc', 'fps': 40},
    '2RN2': {'node_file': f'{BASE}/2rn2_md/2rn2.N', 'edge_file': f'{BASE}/2rn2_md/2rn2.E',
             'gro': f'{BASE}/2rn2_md/2rn2_md_first_frame.gro', 'xtc': f'{BASE}/2rn2_md/2rn2_md_reduced.xtc', 'fps': 100},
    '2CI2': {'node_file': f'{BASE}/2ci2_md/2ci2.N', 'edge_file': f'{BASE}/2ci2_md/2ci2.E',
             'gro': f'{BASE}/2ci2_md/2ci2_md_first_frame.gro', 'xtc': f'{BASE}/2ci2_md/2ci2_md_reduced.xtc', 'fps': 100},
    # 1QJP (OmpA, outer-membrane beta-barrel), added in place of 1OH0.
    '1QJP': {'node_file': f'{BASE}/1qjp_md/1qjp.N', 'edge_file': f'{BASE}/1qjp_md/1qjp.E',
             'gro': f'{BASE}/1qjp_md/1qjp_md_first_frame.gro', 'xtc': f'{BASE}/1qjp_md/1qjp_md_reduced.xtc', 'fps': 100},
}

# The 271-pair benchmark ships with the code (repo root); the per-protein MD
# data lives under BASE, one <pdb>_md/ folder per protein.
DDG_FILE   = _DDG_OVERRIDE or './benchmark_271pairs.ddg'
OUTPUT_DIR = _OUT_OVERRIDE or './output'
# Directory holding the per-pair virtual-edge predictions (ve2_<PDB>.csv) used
# for the McNemar comparison against the best virtual-edge baseline. If these
# files are absent, the code falls back to the static-RCNC clique baseline.
VE_PRED_DIR = os.path.join(BASE, 'virtual_edge')
SKIP_NS    = 0
DDDG_CUT   = 1.0   # kcal/mol threshold for non-additivity

# ── Thresholds (calibrated to biophysical literature) ──────
# ── C1  Clique community (Tier 1, static) ──
C1_COUPLING                  = 0.30   # |Z-DCCM| cutoff for a shared k=3 community
C1_HELIX_FRAC                = 0.75   # helix fraction to classify a helix-dominated community
# ── C2  Direct contact (Tier 2, hybrid) ──
C2_COUPLING                  = 0.45   # |Z-DCCM| cutoff for a direct RING contact
# ── C3  Dynamic packing (Tier 2, MD) ──
C3_SHORT_HELIX_COUPLING      = 0.40   # |Z-DCCM| cutoff, short-range helical adjacency
C3_HELIX_BOTH_ISO_COUPLING   = 1.00   # |Z-DCCM| cutoff, both helical residues isolated
C3_COIL_COUPLING             = 1.20   # |Z-DCCM| cutoff, both residues in coil
C3_HELIX_ONE_ISO_COUPLING    = 1.40   # |Z-DCCM| cutoff, one helical residue isolated
C3_SHORT_HELIX_SEP           = 4      # max sequence separation, short-range helix rule
C3_HELIX_ISO_SEP_MAX         = 30     # max sequence separation, isolated-helix rule
C3_HELIX_BOTH_ISO_SEP_MIN    = 4      # min sequence separation, both-isolated helix rule
C3_HELIX_BOTH_ISO_SEP_MAX    = 30     # max sequence separation, both-isolated helix rule
# ── C4  Sequential backbone (Tier 2, MD) ──
C4_BACKBONE_CORR             = 0.55   # raw DCCM cutoff for sequential neighbours
# ── C5  Network topology (Tier 2, MD) ──
C5_COMMON_NEIGHBOR_COUPLING  = 0.65   # |Z-DCCM| cutoff, shared network neighbour
C5_LINKED_COMM_COUPLING      = 1.00   # |Z-DCCM| cutoff, linked distant communities
C5_HUB_BETWEENNESS           = 0.03   # betweenness-centrality cutoff defining a hub
C5_HUB_SEP_MIN               = 15     # min sequence separation, hub-perturbation rule
C5_LINKED_COMM_SEP_MIN       = 30     # min sequence separation, linked-community rule
CHARGED_AA   = {'R', 'K', 'D', 'E'}
MAJOR_SMALL_AA = {'G', 'A'}
MAJOR_LARGE_AA = {'L', 'I', 'F', 'W', 'Y', 'M'}

# ─────────────────────────────────────────────
#  RULE -> CRITERION MAP
#  predict() returns a fine-grained internal rule name; this dictionary maps
#  each rule to the high-level criterion (C1-C6) used in the manuscript and in
#  Table 2, so that code output and manuscript terminology stay consistent.
#  The 14 active rules are listed in sequential order (R1-R14), grouped under
#  their criterion. The R# in each comment is the manuscript rule index.
# ─────────────────────────────────────────────
CRITERION_OF_RULE = {
    # ── C1  Clique community membership (Tier 1, static) ──
    'R1_C1_same_community':          'C1',
    'R2_C1_helix_community':         'C1',
    'R3_C1_helix_strand_boundary':  'C1',
    # ── C2  Direct contact with dynamic support (Tier 2, hybrid) ──
    'R4_C2_direct_interaction':      'C2',
    # ── C3  Dynamic packing coupling (Tier 2, MD) ──
    'R6_C3_helix_one_isolated':     'C3',
    'R7_C3_helix_both_isolated': 'C3',
    'R5_C3_coil_coupling':              'C3',
    'R8_C3_short_helix_adjacent':     'C3',
    # ── C4  Sequential backbone coupling (Tier 2, MD) ──
    'R9_C4_sequential_backbone':     'C4',
    # ── C5  Network topology coupling (Tier 2, MD) ──
    'R10_C5_hub_perturbation':        'C5',
    'R11_C5_common_neighbor':         'C5',
    'R12_C5_hub_no_dccm':             'C5',
    'R13_C5_linked_community':       'C5',
    # ── C6  Electrostatic coupling (Tier 3, hybrid) ──
    'R14_C6_charged_isolated':       'C6',
    # default
    'none': 'none',
}


# ─────────────────────────────────────────────
#  PREPROCESSING
# ─────────────────────────────────────────────
def preprocess_ddg(ddg_df):
    """No deduplication. Every measured pair is retained as a separate data
    point (following Zhang et al. 2024), giving the full 271-pair benchmark
    reported in the manuscript. Pairs with identical mutation, pH, method and
    dddG are genuine repeated measurements in the benchmark and are kept."""
    return ddg_df.reset_index(drop=True)


# ─────────────────────────────────────────────
#  STEP 1: BUILD NETWORK
# ─────────────────────────────────────────────
def build_network(node_file, edge_file):
    n_df = pd.read_table(node_file, encoding='latin-1')
    e_df = pd.read_table(edge_file, encoding='latin-1')
    n_df.columns = n_df.columns.str.strip()
    e_df.columns = e_df.columns.str.strip()

    sites = {}; energy_index = defaultdict(list); node_sites = []
    G = nx.Graph()
    for _, row in e_df.iterrows():
        s1 = int(row['NodeId1'].split(':')[1])
        s2 = int(row['NodeId2'].split(':')[1])
        itype = row['Interaction'].split(':')[0]
        try:    en = float(row['Energy']) if 'Energy' in e_df.columns else 6.0
        except: en = 6.0
        if itype == 'HBOND':
            if abs(s1-s2) > 4:
                idx = str(sorted([s1,s2])); energy_index[idx].append(en)
                sites[idx]=(s1,s2); node_sites.append(idx)
        else:
            idx = str(sorted([s1,s2])); energy_index[idx].append(en)
            sites[idx]=(s1,s2); node_sites.append(idx)
    for idx in node_sites:
        G.add_edge(sites[idx][0], sites[idx][1], Weight=sum(energy_index[idx]))

    last = int(n_df['Position'].iloc[-1])
    for i in range(1, last):
        G.add_edge(i, i+1)

    communities  = list(k_clique_communities(G, 3))
    r2c          = {r: i for i,c in enumerate(communities) for r in c}
    G_nb         = G.copy()
    for i in range(1, last+1):
        if G_nb.has_edge(i, i+1): G_nb.remove_edge(i, i+1)
    bc = nx.betweenness_centrality(G_nb, normalized=True)

    direct = set()
    for _, row in e_df.iterrows():
        s1 = int(row['NodeId1'].split(':')[1])
        s2 = int(row['NodeId2'].split(':')[1])
        itype = row['Interaction'].split(':')[0]
        if abs(s1-s2) > 2:
            if itype in ('IONIC','PIPISTACK','PICATION'):
                direct.add((min(s1,s2), max(s1,s2)))
            elif itype=='HBOND' and abs(s1-s2)>4:
                direct.add((min(s1,s2), max(s1,s2)))

    linked_comm = defaultdict(set)
    for i in range(len(communities)):
        for j in range(i+1, len(communities)):
            if set(communities[i]) & set(communities[j]):
                linked_comm[i].add(j); linked_comm[j].add(i)

    return communities, r2c, direct, G_nb, bc, linked_comm, last


# ─────────────────────────────────────────────
#  STEP 2: MD FEATURES
# ─────────────────────────────────────────────
def compute_md_features(gro, xtc, fps):
    u      = mda.Universe(gro, xtc)
    resids = np.array([r.resid for r in u.residues])
    rmap   = {r.resid: k for k,r in enumerate(u.residues)}
    start  = int(SKIP_NS * 1000 / fps)

    traj_full = md.load(xtc, top=gro)
    ca_idx    = traj_full.topology.select('name CA')
    traj_ca   = traj_full.atom_slice(ca_idx)[start:]
    traj_ca.superpose(traj_ca, 0)
    coords    = traj_ca.xyz

    delta = coords - coords.mean(axis=0)
    num   = np.einsum('fic,fjc->ij', delta, delta)
    sq    = np.sqrt(np.einsum('fic,fic->i', delta, delta))
    dccm  = np.where(np.outer(sq,sq)>0, num/np.outer(sq,sq), 0.0)

    z     = np.zeros_like(dccm)
    n_res = len(resids)
    for d in range(int(resids.max()-resids.min())+1):
        pairs = [(i,j) for i in range(n_res) for j in range(n_res)
                 if i!=j and abs(resids[i]-resids[j])==d]
        if len(pairs) < 5: continue
        vals = [dccm[i,j] for i,j in pairs]
        mu, sigma = np.mean(vals), max(np.std(vals), 0.01)
        for i,j in pairs: z[i,j] = (dccm[i,j]-mu)/sigma

    rmsf     = coords.std(axis=0).mean(axis=1)
    traj_skip = traj_full[start:]
    dssp_raw  = md.compute_dssp(traj_skip[::10], simplified=True)
    ss = [Counter(dssp_raw[:,k]).most_common(1)[0][0] for k in range(dssp_raw.shape[1])]

    print(f"  MD: {len(resids)} residues, {len(traj_ca)} frames, RMSF={rmsf.mean():.4f}nm")
    return resids, rmap, dccm, z, rmsf, ss


# ─────────────────────────────────────────────
#  STEP 3: PREDICTION (six criteria C1-C6, three tiers)
# ─────────────────────────────────────────────
def get_ss(site, rmap, ss):
    return ss[rmap[site]] if site in rmap else 'C'

def is_charged_orig(mp):
    l = ''.join(c for c in mp if c.isalpha())
    return l[0] in CHARGED_AA if l else False

def is_charged_mut(mp):
    l = ''.join(c for c in mp if c.isalpha())
    return l[-1] in CHARGED_AA if l else False

def has_major_vol_change(mp):
    l = ''.join(c for c in mp if c.isalpha())
    if len(l) < 2: return False
    return ((l[0] in MAJOR_SMALL_AA and l[-1] in MAJOR_LARGE_AA) or
            (l[0] in MAJOR_LARGE_AA and l[-1] in MAJOR_SMALL_AA))


def predict(s1, s2, mut_str, communities, r2c, G_nb, bc,
            direct, linked_comm, protein_last, rmap, ss, dccm, z, rmsf,
            thr=None):
    # The three LOPOCV-tunable thresholds are taken from the optional
    # `thr` dict when supplied (during cross-validation) and otherwise from the
    # module-level defaults.
    thr  = thr or {}
    _zr1 = thr.get('C1_COUPLING', C1_COUPLING)
    _zr4 = thr.get('C3_COIL_COUPLING', C3_COIL_COUPLING)
    _bc  = thr.get('C5_HUB_BETWEENNESS', C5_HUB_BETWEENNESS)
    seq  = abs(s1-s2)
    c1   = r2c.get(s1,-1); c2 = r2c.get(s2,-1)
    same = (c1!=-1 and c2!=-1 and c1==c2)
    ss1  = get_ss(s1,rmap,ss); ss2 = get_ss(s2,rmap,ss)
    pair = (min(s1,s2), max(s1,s2))
    i    = rmap.get(s1,-1); j = rmap.get(s2,-1)
    zv   = abs(z[i,j])    if i>=0 and j>=0 else 0.0
    dc   = abs(dccm[i,j]) if i>=0 and j>=0 else 0.0
    bc1  = bc.get(s1,0); bc2 = bc.get(s2,0)
    n1   = set(G_nb.neighbors(s1)) if s1 in G_nb else set()
    n2   = set(G_nb.neighbors(s2)) if s2 in G_nb else set()
    common      = n1 & n2
    adj_comm    = (c1!=-1 and c2!=-1 and c1!=c2)
    both_coil   = (ss1=='C' and ss2=='C')
    both_helix  = (ss1=='H' and ss2=='H')
    parts       = str(mut_str).split(',')
    ch1_keeps   = is_charged_orig(parts[0]) and is_charged_mut(parts[0]) if parts else False
    ch2_keeps   = is_charged_orig(parts[1]) and is_charged_mut(parts[1]) if len(parts)>1 else False
    vol1        = has_major_vol_change(parts[0]) if parts else False
    vol2        = has_major_vol_change(parts[1]) if len(parts)>1 else False
    vol_sig     = vol1 or vol2
    both_iso    = (c1==-1 and c2==-1)
    one_iso     = (c1==-1 or c2==-1)
    comm_linked = (c1!=-1 and c2!=-1 and c1!=c2 and c2 in linked_comm.get(c1,set()))

    # ── TIER 1: Static Network ──────────────────────────────
    # C1: Clique community membership
    if same:
        helix_frac = sum(1 for r in communities[c1]
                        if get_ss(r,rmap,ss)=='H') / max(len(communities[c1]),1)
        if helix_frac >= C1_HELIX_FRAC:
            return 'Additive', 'R2_C1_helix_community'
        if (ss1=='E' and ss2=='H') or (ss1=='H' and ss2=='E'):
            return 'Additive', 'R3_C1_helix_strand_boundary'
        if zv > _zr1:
            return 'Non-additive', 'R1_C1_same_community'
        return 'Additive', 'R1_C1_same_community_low'

    # C2: Direct contact + dynamic support
    if pair in direct and seq > 2 and zv > C2_COUPLING:
        if seq <= 3 and both_helix:
            pass
        else:
            return 'Non-additive', 'R4_C2_direct_interaction'

    # ── TIER 2: MD Dynamic Filters ──────────────────────────

    # C3: Coil dynamic coupling
    if both_coil and seq > 2 and zv > _zr4:
        return 'Non-additive', 'R5_C3_coil_coupling'

    # C3: Helix one isolated
    if both_helix and one_iso and zv > C3_HELIX_ONE_ISO_COUPLING and seq <= C3_HELIX_ISO_SEP_MAX:
        return 'Non-additive', 'R6_C3_helix_one_isolated'

    # C3: Helix both isolated + volume change
    if (both_helix and both_iso and zv > C3_HELIX_BOTH_ISO_COUPLING and
            C3_HELIX_BOTH_ISO_SEP_MIN < seq <= C3_HELIX_BOTH_ISO_SEP_MAX and vol_sig):
        return 'Non-additive', 'R7_C3_helix_both_isolated'

    # Short helix adjacent community
    if seq <= C3_SHORT_HELIX_SEP and both_helix and adj_comm and zv > C3_SHORT_HELIX_COUPLING:
        return 'Non-additive', 'R8_C3_short_helix_adjacent'

    # C4: Sequential + high DCCM
    if seq <= 2 and dc > C4_BACKBONE_CORR and both_iso:
        return 'Non-additive', 'R9_C4_sequential_backbone'

    # C5: Common neighbor
    if len(common) > 0 and seq > 5 and not adj_comm and zv > C5_COMMON_NEIGHBOR_COUPLING:
        return 'Non-additive', 'R11_C5_common_neighbor'

    # C5: Hub betweenness
    if (bc1 > _bc or bc2 > _bc) and seq > C5_HUB_SEP_MIN:
        if zv < 0.15:
            return 'Additive', 'R12_C5_hub_no_dccm'
        return 'Non-additive', 'R10_C5_hub_perturbation'

    # C5: Linked community
    if comm_linked and zv > C5_LINKED_COMM_COUPLING and seq > C5_LINKED_COMM_SEP_MIN:
        return 'Non-additive', 'R13_C5_linked_community'

    # ── TIER 3: Physicochemical ─────────────────────────────
    # C6: Electrostatic coupling
    if ch1_keeps and ch2_keeps and both_iso and seq >= 8:
        return 'Non-additive', 'R14_C6_charged_isolated'

    return 'Additive', 'none'


# ─────────────────────────────────────────────
#  METRICS
# ─────────────────────────────────────────────
def compute_metrics(df):
    yt = (df['exp'] =='Non-additive').astype(int).values
    yp = (df['pred']=='Non-additive').astype(int).values
    tp=int(((yt==1)&(yp==1)).sum()); tn=int(((yt==0)&(yp==0)).sum())
    fp=int(((yt==0)&(yp==1)).sum()); fn=int(((yt==1)&(yp==0)).sum())
    sens = tp/(tp+fn) if tp+fn>0 else 0
    spec = tn/(tn+fp) if tn+fp>0 else 0
    prec = tp/(tp+fp) if tp+fp>0 else 0
    f1   = 2*tp/(2*tp+fp+fn) if 2*tp+fp+fn>0 else 0
    bacc = (sens+spec)/2
    md_  = np.sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn))
    mcc  = (tp*tn-fp*fn)/md_ if md_>0 else 0
    return dict(TP=tp, TN=tn, FP=fp, FN=fn, N=len(df),
                NA=int(yt.sum()), ADD=int((yt==0).sum()),
                Sensitivity=sens, Specificity=spec, Precision=prec,
                F1=f1, BalancedAcc=bacc, MCC=mcc)


def bootstrap_ci(df, key, n=1000, seed=42):
    rng  = np.random.default_rng(seed)
    vals = [compute_metrics(df.iloc[rng.integers(0,len(df),len(df))])[key] for _ in range(n)]
    return np.percentile(vals, 2.5), np.percentile(vals, 97.5)


def permutation_test(yt, yp, n=1000, seed=42):
    rng = np.random.default_rng(seed)
    obs = (yp==yt).mean()
    return sum(1 for _ in range(n) if (yp==rng.permutation(yt)).mean() >= obs) / n


def load_virtual_edge_predictions(ve_dir, ve_z=2.75, ve_d=25):
    """
    Load per-pair predictions of the best virtual-edge configuration
    (|Z-DCCM| > ve_z, Cα distance <= ve_d) from ve2_<PDB>.csv files.
    Returns {(PDB, min_resid, max_resid): 1 if predicted non-additive else 0}.
    These predictions depend only on residue position, so they are later
    broadcast to every measurement row of the same residue pair.
    """
    import glob
    lookup = {}
    for f in glob.glob(os.path.join(ve_dir, 've2_*.csv')):
        pdb = os.path.basename(f).replace('ve2_', '').replace('.csv', '')
        if pdb == '1OH0':
            continue
        try:
            d = pd.read_csv(f, keep_default_na=False)
        except Exception:
            continue
        sub = d[(d['zthr'] == ve_z) & (d['dist'] == ve_d)]
        for _, r in sub.iterrows():
            key = (pdb,) + tuple(sorted((int(r['s1']), int(r['s2']))))
            lookup[key] = 1 if r['pred'] == 'NA' else 0
    return lookup


def mcnemar_test(yt, yp, base_pred):
    """
    McNemar test of DynRCNC (yp) against a baseline model (base_pred),
    relative to the experimental labels (yt). All arrays are aligned 0/1
    vectors of equal length. b = pairs where DynRCNC is correct and the
    baseline is wrong; c = the reverse.
    """
    yt = np.asarray(yt); yp = np.asarray(yp); yb = np.asarray(base_pred)
    b = int(((yp == yt) & (yb != yt)).sum())
    c = int(((yp != yt) & (yb == yt)).sum())
    if b + c == 0:
        return 1.0, b, c
    return float(chi2_dist.sf((abs(b - c) - 1) ** 2 / (b + c), df=1)), b, c


def threshold_sensitivity(protein_dfs):
    all_df = pd.concat(protein_dfs.values()).copy()
    results = []
    for t in np.linspace(0.5, 1.5, 10):
        df_t = all_df.copy()
        df_t['exp'] = df_t['dddG'].abs().apply(
            lambda x: 'Non-additive' if x >= t else 'Additive')
        m = compute_metrics(df_t)
        results.append({'threshold': round(t,2), 'MCC': m['MCC']})
    return pd.DataFrame(results)


def lopocv(protein_cfgs, ddg_df):
    """Nested LOPOCV; the three tunable thresholds are re-optimized per fold."""
    import itertools, math as _math

    GRID = {'C1_COUPLING': [0.20,0.30,0.40], 'C3_COIL_COUPLING': [1.00,1.20,1.40], 'C5_HUB_BETWEENNESS': [0.02,0.03,0.05]}

    def mcc_score(tp,tn,fp,fn):
        d = _math.sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn))
        return (tp*tn-fp*fn)/d if d>0 else 0.0

    print('\n  Running GENUINE NESTED LOPOCV ...')
    all_pdbs = list(protein_cfgs.keys())

    print('  Pre-loading all protein features...')
    pfeatures = {}
    for pdb, cfg in protein_cfgs.items():
        try:
            comm,r2c_,direct_,G_nb_,bc_,lc_,pl_ = build_network(cfg['node_file'],cfg['edge_file'])
            res_,rmap_,dccm_,z_,rmsf_,ss_ = compute_md_features(cfg['gro'],cfg['xtc'],cfg['fps'])
            pfeatures[pdb] = (comm,r2c_,direct_,G_nb_,bc_,lc_,pl_,res_,rmap_,dccm_,z_,rmsf_,ss_)
            print(f'    {pdb}: {len(res_)} residues, {len(z_)} zdccm dim OK')
        except Exception as e:
            print(f'    {pdb}: FAILED: {e}')

    rows = []
    for test_pdb in all_pdbs:
        if test_pdb not in pfeatures: continue
        train_pdbs = [p for p in all_pdbs if p!=test_pdb and p in pfeatures]
        print(f'\n  Fold: held-out={test_pdb} | train={train_pdbs}')

        best_params = None; best_mcc = -99.0
        combos = list(itertools.product(*GRID.values()))
        keys   = list(GRID.keys())
        print(f'    Grid search ({len(combos)} combos)...', end='', flush=True)

        for combo in combos:
            params = dict(zip(keys,combo))
            tp=tn=fp=fn=0
            for tr_pdb in train_pdbs:
                (comm,r2c_,direct_,G_nb_,bc_,lc_,pl_,res_,rmap_,dccm_,z_,rmsf_,ss_) = pfeatures[tr_pdb]
                for _, row in ddg_df[ddg_df['PDB']==tr_pdb].iterrows():
                    parts = str(row['Mutation_Double']).split(',')
                    nums  = [int(''.join(c for c in p if c.isdigit())) for p in parts if any(c.isdigit() for c in p)]
                    if len(nums)!=2: continue
                    s1,s2 = nums
                    exp  = 'Non-additive' if abs(float(row['dddG']))>=DDDG_CUT else 'Additive'
                    pred,_ = predict(s1,s2,row['Mutation_Double'],comm,r2c_,G_nb_,bc_,direct_,lc_,pl_,rmap_,ss_,dccm_,z_,rmsf_,
                                     thr=params)
                    if   exp=='Non-additive' and pred=='Non-additive': tp+=1
                    elif exp=='Additive'     and pred=='Additive':     tn+=1
                    elif exp=='Additive'     and pred=='Non-additive': fp+=1
                    else:                                               fn+=1
            tr_mcc = mcc_score(tp,tn,fp,fn)
            if tr_mcc > best_mcc: best_mcc=tr_mcc; best_params=params

        print(f' best={best_params} train_MCC={best_mcc:.3f}')

        (comm,r2c_,direct_,G_nb_,bc_,lc_,pl_,res_,rmap_,dccm_,z_,rmsf_,ss_) = pfeatures[test_pdb]
        tp=tn=fp=fn=0
        for _, row in ddg_df[ddg_df['PDB']==test_pdb].iterrows():
            parts = str(row['Mutation_Double']).split(',')
            nums  = [int(''.join(c for c in p if c.isdigit())) for p in parts if any(c.isdigit() for c in p)]
            if len(nums)!=2: continue
            s1,s2=nums
            exp  = 'Non-additive' if abs(float(row['dddG']))>=DDDG_CUT else 'Additive'
            pred,_ = predict(s1,s2,row['Mutation_Double'],comm,r2c_,G_nb_,bc_,direct_,lc_,pl_,rmap_,ss_,dccm_,z_,rmsf_,
                             thr=best_params)
            if   exp=='Non-additive' and pred=='Non-additive': tp+=1
            elif exp=='Additive'     and pred=='Additive':     tn+=1
            elif exp=='Additive'     and pred=='Non-additive': fp+=1
            else:                                               fn+=1

        # Test-fold metrics come directly from the manually counted
        # confusion matrix below.
        test_mcc = mcc_score(tp,tn,fp,fn)
        test_sens = tp/(tp+fn) if tp+fn>0 else 0
        test_spec = tn/(tn+fp) if tn+fp>0 else 0
        test_bacc = (test_sens+test_spec)/2
        N_test    = len(ddg_df[ddg_df['PDB']==test_pdb])
        print(f'    → Test MCC={test_mcc:.3f} Sens={test_sens*100:.1f}% Spec={test_spec*100:.1f}%')
        rows.append({'PDB':test_pdb,'N':N_test,'MCC':test_mcc,
                     'Sensitivity':test_sens,'Specificity':test_spec,'BalancedAcc':test_bacc,
                     'best_ZDCCM_R1':best_params['C1_COUPLING'],
                     'best_ZDCCM_R4':best_params['C3_COIL_COUPLING'],
                     'best_BC_THRESH':best_params['C5_HUB_BETWEENNESS']})

    return pd.DataFrame(rows)


# ─────────────────────────────────────────────
#  PER-PROTEIN RUN
# ─────────────────────────────────────────────
def run_protein(pdb, cfg, ddg_df):
    print(f"\n{'─'*55}\n  {pdb}\n{'─'*55}")
    communities, r2c, direct, G_nb, bc, linked_comm, last = build_network(cfg['node_file'], cfg['edge_file'])
    print(f"  Communities: {len(communities)}")
    resids, rmap, dccm, z, rmsf, ss = compute_md_features(cfg['gro'], cfg['xtc'], cfg['fps'])

    sub = ddg_df[ddg_df['PDB']==pdb].copy()
    if len(sub)==0: print("  No DDG data."); return None

    rows = []
    for _, row in sub.iterrows():
        parts = str(row['Mutation_Double']).split(',')
        nums  = [int(''.join(c for c in p if c.isdigit())) for p in parts if any(c.isdigit() for c in p)]
        if len(nums)!=2:
            print(f"  [WARN] {pdb}: skipping malformed mutation string "
                  f"'{row['Mutation_Double']}' (could not parse two residue numbers)")
            continue
        s1,s2 = nums
        # Warn (do not silently zero) if a mutation site is absent from the
        # MD residue map, e.g. a missing/non-standard residue in the structure.
        if s1 not in rmap or s2 not in rmap:
            missing = [r for r in (s1,s2) if r not in rmap]
            print(f"  [WARN] {pdb}: residue(s) {missing} for '{row['Mutation_Double']}' "
                  f"not in MD structure; Z-DCCM/DCCM will be treated as 0 for this pair")
        dddg  = float(row['dddG'])
        exp   = 'Non-additive' if abs(dddg)>=DDDG_CUT else 'Additive'
        pred, mech = predict(s1,s2,row['Mutation_Double'],communities,r2c,G_nb,bc,
                             direct,linked_comm,last,rmap,ss,dccm,z,rmsf)
        i_,j_ = rmap.get(s1,-1), rmap.get(s2,-1)
        rows.append({'PDB':pdb, 'mutation':row['Mutation_Double'],
                     'pH':row.get('pH',''), 'Method':row.get('Method',''),
                     'dddG':dddg, 'exp':exp, 'pred':pred, 'correct':exp==pred,
                     'mechanism':mech, 'seq_dist':abs(s1-s2),
                     'res_i':s1, 'res_j':s2,
                     'comm_i':r2c.get(s1,-1), 'comm_j':r2c.get(s2,-1),
                     'zdccm':round(abs(z[i_,j_]),4) if i_>=0 and j_>=0 else 0,
                     'ss_i':get_ss(s1,rmap,ss), 'ss_j':get_ss(s2,rmap,ss)})

    df = pd.DataFrame(rows)
    m  = compute_metrics(df)

    print(f"\n  RESULTS: {pdb}")
    print(f"  TP={m['TP']} TN={m['TN']} FP={m['FP']} FN={m['FN']} N={m['N']}")
    print(f"  Sensitivity : {m['Sensitivity']*100:.1f}%  ({m['TP']}/{m['NA']})")
    print(f"  Specificity : {m['Specificity']*100:.1f}%  ({m['TN']}/{m['ADD']})")
    print(f"  Precision   : {m['Precision']*100:.1f}%")
    print(f"  F1          : {m['F1']:.3f}")
    print(f"  MCC         : {m['MCC']:.3f}")

    cols = ['mutation','pH','Method','dddG','exp','pred','correct',
            'mechanism','seq_dist','comm_i','comm_j','zdccm','ss_i','ss_j']

    print(f"\n  ── Non-additive pairs ({m['NA']}) ──")
    print(df[df['exp']=='Non-additive'][cols].to_string(index=False))

    print(f"\n  ── Additive pairs ({m['ADD']}) ──")
    print(df[df['exp']=='Additive'][cols].to_string(index=False))

    print(f"\n  Mechanism breakdown:")
    print(df[df['pred']=='Non-additive']['mechanism'].value_counts().to_string())

    return df, m


# ─────────────────────────────────────────────
#  MAIN
# ─────────────────────────────────────────────
def main():
    t0 = time.time()
    print(f"\n{'='*62}")
    print(f"  DynRCNC")
    print(f"  14 active rules | 3 tiers | 8 proteins")
    print(f"{'='*62}")

    os.makedirs(OUTPUT_DIR, exist_ok=True)
    ddg_df = pd.read_csv(DDG_FILE, sep='\t', on_bad_lines='skip')
    ddg_df.columns = ddg_df.columns.str.strip()
    ddg_df = preprocess_ddg(ddg_df)
    print(f"\nDDG loaded: {len(ddg_df)} pairs")

    protein_dfs = {}; summary = []; skipped = []

    for pdb, cfg in PROTEINS.items():
        missing = [k for k in ['node_file','edge_file','gro','xtc'] if not os.path.exists(cfg[k])]
        if missing:
            print(f"\n  {pdb}: SKIPPED, missing: {missing}")
            skipped.append(pdb); continue
        try:
            result = run_protein(pdb, cfg, ddg_df)
            if result is None: skipped.append(pdb); continue
            df, m = result
            protein_dfs[pdb] = df
            summary.append({'PDB':pdb, **m})
        except Exception as ex:
            import traceback
            print(f"\n  {pdb}: ERROR: {ex}"); traceback.print_exc()
            skipped.append(pdb)

    if not protein_dfs:
        print("\nNo proteins processed."); return

    all_df = pd.concat(protein_dfs.values(), ignore_index=True)
    comb   = compute_metrics(all_df)
    summary.append({'PDB':'COMBINED', **comb})

    yt = (all_df['exp'] =='Non-additive').astype(int).values
    yp = (all_df['pred']=='Non-additive').astype(int).values

    print(f"\n  Computing bootstrap CIs...")
    ci     = {k: bootstrap_ci(all_df,k) for k in ['MCC','F1','Sensitivity','Specificity']}
    p_perm = permutation_test(yt, yp)

    # McNemar test against the best virtual-edge configuration (manuscript comparison).
    # Virtual-edge predictions are loaded per residue pair and broadcast to every
    # measurement row, so the test spans all benchmark rows. If the ve2_*.csv files
    # are not found, fall back to the static-RCNC clique baseline and say so.
    ve_lookup = load_virtual_edge_predictions(VE_PRED_DIR)
    if ve_lookup:
        yv = []
        for _, r in all_df.iterrows():
            key  = (r['PDB'],) + tuple(sorted((int(r['res_i']), int(r['res_j']))))
            yv.append(ve_lookup.get(key, 0))
        p_mcn, b_mcn, c_mcn = mcnemar_test(yt, yp, np.array(yv))
        mcn_baseline = f"virtual edge (Z>2.75, d<=25 A); b={b_mcn}, c={c_mcn}"
    else:
        y_static = np.array([1 if (ci_!=-1 and cj_!=-1 and ci_==cj_) else 0
                             for ci_,cj_ in zip(all_df['comm_i'].values, all_df['comm_j'].values)])
        p_mcn, b_mcn, c_mcn = mcnemar_test(yt, yp, y_static)
        mcn_baseline = (f"static RCNC clique baseline (ve2_*.csv not found in {VE_PRED_DIR}); "
                        f"b={b_mcn}, c={c_mcn}")

    protein_cfgs_filtered = {pdb: PROTEINS[pdb] for pdb in protein_dfs}
    lodf   = lopocv(protein_cfgs_filtered, ddg_df)
    ts_df  = threshold_sensitivity(protein_dfs)
    gen_gap = abs(comb['MCC'] - (lodf['MCC']*lodf['N']).sum()/lodf['N'].sum())

    # ── COMBINED RESULTS ────────────────────────────────────
    print(f"\n{'='*62}")
    print(f"  COMBINED: {len(protein_dfs)} proteins, {len(all_df)} pairs")
    print(f"{'='*62}")
    print(f"  TP={comb['TP']} TN={comb['TN']} FP={comb['FP']} FN={comb['FN']}")
    for met in ['Sensitivity','Specificity','Precision','F1','MCC']:
        lo,hi = ci.get(met,(None,None))
        ci_str = f"  95%CI [{lo:.3f}–{hi:.3f}]" if lo is not None else ""
        print(f"  {met:14}: {comb[met]*100:.2f}%{ci_str}")

    print(f"\n  Statistical Tests:")
    print(f"  Permutation p    : {p_perm:.4f} {'✓' if p_perm<0.05 else '✗'}")
    print(f"  McNemar p        : {p_mcn:.4f}  [vs {mcn_baseline}]")
    print(f"  Generalization Gap (MCC): {gen_gap:.4f}")

    # ── PER-PROTEIN TABLE ───────────────────────────────────
    print(f"\n{'─'*78}")
    print(f"  {'PDB':8} {'N':>4} {'NA':>4} {'ADD':>4} {'Sens%':>7} {'Spec%':>7} {'F1':>6} {'MCC':>6}")
    print(f"{'─'*78}")
    for r in summary:
        print(f"  {r['PDB']:8} {r['N']:>4} {r['NA']:>4} {r['ADD']:>4} "
              f"{r['Sensitivity']*100:>7.1f} {r['Specificity']*100:>7.1f} "
              f"{r['F1']:>6.3f} {r['MCC']:>6.3f}")

    # ── LOPOCV ──────────────────────────────────────────────
    print(f"\n  GENUINE NESTED LOPOCV:")
    print(f"  {'Test':8} {'N':>4} {'Sens%':>7} {'Spec%':>7} {'MCC':>6} {'BestR1':>8} {'BestR4':>8} {'BestBC':>8}")
    for _,r in lodf.iterrows():
        print(f"  {r['PDB']:8} {r['N']:>4} {r['Sensitivity']*100:>7.1f} "
              f"{r['Specificity']*100:>7.1f} {r['MCC']:>6.3f} "
              f"{r.get('best_ZDCCM_R1','?'):>8} {r.get('best_ZDCCM_R4','?'):>8} {r.get('best_BC_THRESH','?'):>8}")
    w = lodf['N']/lodf['N'].sum()
    w_mcc = (lodf['MCC']*w).sum()
    print(f"  Weighted LOPOCV MCC: {w_mcc:.3f}")

    # ── THRESHOLD SENSITIVITY ───────────────────────────────
    print(f"\n  Threshold Sensitivity (MCC at different |dddG| cutoffs):")
    for _,r in ts_df.iterrows():
        bar = '█'*int(r['MCC']*20)
        print(f"    {r['threshold']:.1f} kcal/mol: MCC={r['MCC']:.3f} {bar}")

    # ── SAVE CSV ────────────────────────────────────────────
    all_df.to_csv(f'{OUTPUT_DIR}/all_predictions.csv', index=False)
    pd.DataFrame(summary).to_csv(f'{OUTPUT_DIR}/summary.csv', index=False)
    lodf.to_csv(f'{OUTPUT_DIR}/lopocv.csv', index=False)
    ts_df.to_csv(f'{OUTPUT_DIR}/threshold_sensitivity.csv', index=False)

    if skipped: print(f"\n  Skipped proteins: {skipped}")
    print(f"\n  Results saved → {os.path.abspath(OUTPUT_DIR)}/")
    print(f"  Runtime: {time.time()-t0:.1f}s")
    print(f"{'='*62}\n")


if __name__ == '__main__':
    main()
