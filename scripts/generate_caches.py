"""
generate_caches.py — Pre-compute expensive notebook computations.

Run from the workshop root:
    conda run -n cheminf python scripts/generate_caches.py

Outputs:
    data/cache_umap_coords.npy       — UMAP 2D embedding of all 1 303 compounds
    data/cache_sim_matrix_50.npy     — 50×50 Tanimoto similarity matrix (Part E subset)
"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

import numpy as np
import pandas as pd
from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit.Chem.rdFingerprintGenerator import GetMorganGenerator
import umap
from tqdm import tqdm
from chem_utils import numpy_to_rdkit_fp

# ── Constants (must match Cell A1) ───────────────────────────────────────────
RANDOM_SEED = 42
FP_RADIUS   = 2
FP_NBITS    = 2048
N_SUBSET    = 50

DATA_DIR = os.path.join(os.path.dirname(__file__), '..', 'data')

UMAP_CACHE = os.path.join(DATA_DIR, 'cache_umap_coords.npy')
SIM_CACHE  = os.path.join(DATA_DIR, 'cache_sim_matrix_50.npy')

# ── Load & standardise ────────────────────────────────────────────────────────
print("Loading dataset...")
df = pd.read_csv(os.path.join(DATA_DIR, 'malaria_box_afrodb_combined.csv'))

RDLogger.DisableLog('rdApp.*')
_lfc = rdMolStandardize.LargestFragmentChooser(preferOrganic=True)

def smiles_to_mol(smi):
    if not isinstance(smi, str) or not smi.strip():
        return None
    return Chem.MolFromSmiles(smi.strip())

df['mol'] = df['SMILES'].apply(smiles_to_mol)
df['mol'] = df['mol'].apply(lambda m: _lfc.choose(m) if m else None)
df_clean  = df[df['mol'].notna()].copy().reset_index(drop=True)
print(f"  Molecules after standardisation: {len(df_clean)}")

# ── Fingerprints ──────────────────────────────────────────────────────────────
print("Computing fingerprints...")
gen = GetMorganGenerator(radius=FP_RADIUS, fpSize=FP_NBITS)
fp_matrix = np.vstack([
    gen.GetFingerprintAsNumPy(m).astype(np.uint8)
    for m in tqdm(df_clean['mol'], desc='  ECFP4')
])
print(f"  fp_matrix shape: {fp_matrix.shape}")

# ── Similarity matrix (50-compound subset) ────────────────────────────────────
if os.path.exists(SIM_CACHE):
    print(f"  [SKIP] {SIM_CACHE} already exists.")
else:
    print("Computing 50×50 Tanimoto similarity matrix...")
    fp_list    = [numpy_to_rdkit_fp(fp_matrix[i]) for i in range(N_SUBSET)]
    sim_matrix = np.array([
        DataStructs.BulkTanimotoSimilarity(fp_list[i], fp_list)
        for i in range(N_SUBSET)
    ])
    np.save(SIM_CACHE, sim_matrix)
    print(f"  ✅ Saved {SIM_CACHE}  shape={sim_matrix.shape}")

# ── UMAP ──────────────────────────────────────────────────────────────────────
if os.path.exists(UMAP_CACHE):
    print(f"  [SKIP] {UMAP_CACHE} already exists.")
else:
    print("Computing UMAP (this takes ~30 s)...")
    reducer = umap.UMAP(
        n_neighbors  = 15,
        min_dist     = 0.1,
        metric       = 'jaccard',
        n_components = 2,
        random_state = RANDOM_SEED
    )
    umap_coords = reducer.fit_transform(fp_matrix)
    np.save(UMAP_CACHE, umap_coords)
    print(f"  ✅ Saved {UMAP_CACHE}  shape={umap_coords.shape}")

print("\nDone. Both cache files are ready.")
