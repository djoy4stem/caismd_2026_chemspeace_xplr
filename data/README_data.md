# Dataset Guide

The raw dataset files are **included in this repository** — no downloading is needed.
Before opening the notebook, generate the processed CSVs by running:

```bash
# From the workshop root, with the cheminf conda env active:
python scripts/prepare_datasets.py
```

---

## Files in this folder

```
data/
├── MalariaBox400compoundsDec2014.xls   ← Raw: MMV Malaria Box (400 compounds)
├── AfroDB_3D.sdf                       ← Raw: AfroDb (954 natural products, 3D conformers)
├── pubchem_aid2302_2k.csv              ← Raw: PubChem AID 2302 screen (2 000 compounds, binary labels)
├── cache_umap_coords.npy               ← Pre-computed: UMAP 2D embedding (1303 × 2, float32)
├── cache_sim_matrix_50.npy             ← Pre-computed: 50×50 Tanimoto similarity matrix (float64)
└── README_data.md                      ← This file
```

> ⚙️ `malaria_box.csv`, `afrodb_subset.csv`, and `malaria_box_afrodb_combined.csv` are not listed
> above because they are **generated files** — produced by running `scripts/prepare_datasets.py`
> from the raw sources above. Run that script once before opening the notebook.

### Pre-computed cache files

`cache_umap_coords.npy` and `cache_sim_matrix_50.npy` are committed to the repository so that
workshop participants on slow hardware (or Google Colab) can skip the expensive computations
and load results in < 1 second.

**What they contain:**

| File | Cell | Shape | Size | Time saved |
|---|---|---|---|---|
| `cache_umap_coords.npy` | G1 | (1303, 2) float32 | ~10 KB | ~30 s laptop, ~2 min slow hardware |
| `cache_sim_matrix_50.npy` | E1 | (50, 50) float64 | ~20 KB | ~0.1 s (negligible, but consistent) |

**How to regenerate** (e.g. after changing `RANDOM_SEED` or `FP_RADIUS`):
```bash
# Delete the stale file(s), then either:
rm data/cache_umap_coords.npy data/cache_sim_matrix_50.npy

# Option A — run the standalone script:
conda run -n cheminf python scripts/generate_caches.py

# Option B — re-run cells E1 and G1 in the executed notebook;
# the cache-or-compute logic will detect the missing file and recompute.
```

---

## Dataset 1 — MMV Malaria Box (`MalariaBox400compoundsDec2014.xls`)

**What it is:**  
400 structurally diverse compounds assembled by the Medicines for Malaria Venture (MMV) and GSK,
all with confirmed antiplasmodial activity against *Plasmodium falciparum* (IC50 ≤ 1 µM). One of
the most widely cited open datasets in antiparasitic drug discovery.

**Source:** MMV supporting information page  
https://www.mmv.org/mmv-open/malaria-box/malaria-box-supporting-information

**File format:** Excel (`.xls`), single sheet `vortex_sheet`, 400 rows × 18 columns

**Key columns:**

| Column | Description |
|---|---|
| `HEOS_COMPOUND_ID` | MMV compound identifier (e.g. `MMV019066`) |
| `Smiles` | SMILES string |
| `EC50_nM` | Anti-*P. falciparum* 3D7 EC50 in nM |
| `percent_inh @ 2 uM` | % inhibition at 2 µM |
| `percent_inh @ 5 uM` | % inhibition at 5 µM |
| `Molecular_Weight` | Precomputed MW |
| `ALogP` | Precomputed ALogP |
| `Num_H_Donors` | H-bond donors |
| `Ro5_ViolationCount` | Lipinski Rule of 5 violations |
| `Set` | `Drug-like` or `Probe-like` |
| `source` | Originating institution (GSK, StJude, etc.) |
| `ChEMBL_NTD_ID` | ChEMBL-NTD identifier |

**License:** CC BY 3.0 — cite as: Spangenberg *et al.*, *Malar. J.* (2013).

---

## Dataset 2 — AfroDb (`AfroDB_3D.sdf`)

**What it is:**  
Dataset S1 from the original AfroDb publication (Ntie-Kang *et al.*, *PLoS ONE*, 2013).
954 African natural products with 3D conformations optimised using QikProp (Schrödinger).

**Source:**  
Ntie-Kang F, *et al.* "AfroDb: A Select Highly Potent and Diverse Natural Product Library
from African Medicinal Plants." *PLoS ONE* 8(10): e78085 (2013).  
https://doi.org/10.1371/journal.pone.0078085

> **Note:** The AfroDb website (`african-compounds.org`) is offline as of 2025. The data
> is preserved in the paper's Supporting Information (Dataset S1).

**File format:** SDF with 3D coordinates and QikProp property annotations, 954 valid molecules

**Key SDF properties (QikProp descriptors):**

| Property | Description |
|---|---|
| `s_m_entry_name` | AfroDb entry name (e.g. `AfroDb.1`) |
| `r_qp_mol_MW` | Molecular weight |
| `r_qp_QPlogPo/w` | Predicted octanol/water logP |
| `r_qp_PSA` | Polar surface area |
| `r_qp_donorHB` | H-bond donors |
| `r_qp_accptHB` | H-bond acceptors |
| `r_qp_QPlogS` | Predicted aqueous solubility |
| `r_qp_QPlogBB` | Blood-brain barrier permeability |
| `r_qp_PercentHumanOralAbsorption` | Predicted human oral absorption (%) |
| `i_qp_#rotor` | Rotatable bonds |
| `i_qp_RuleOfFive` | Lipinski Rule of 5: pass (1) / fail (0) |
| `i_qp_CNS` | CNS activity prediction |

**No bioactivity values** are included — structures and predicted ADMET properties only.

---

## Dataset 3 — PubChem AID 2302 (`pubchem_aid2302_2k.csv`)

**What it is:**  
A 2 000-compound random subset of PubChem BioAssay AID 2302: a whole-cell *Plasmodium falciparum*
Dd2-strain growth inhibition assay run at 2 µM. Compounds are labelled **Active** or **Inactive**,
making this dataset ideal for activity-cliff detection (Extension 3) — the binary labels give a
clear cliff signal that the narrow pIC50 range of the Malaria Box cannot provide.

**Source:** PubChem BioAssay — https://pubchem.ncbi.nlm.nih.gov/bioassay/2302

**File format:** CSV, 2 000 rows

**Columns:**

| Column | Description |
|---|---|
| `CID` | PubChem Compound ID |
| `SMILES` | Isomeric SMILES from PubChem |
| `activity` | `"Active"` or `"Inactive"` (2 µM cut-off) |

**Used in:** Extension 3 only — loaded directly in the Extension 3 code cell.

---

## Processed CSV files (notebook-ready)

Generated by `prepare_datasets.py`. Columns are harmonised across all sources:

| Column | Source |
|---|---|
| `COMPOUND_ID` | MMV ID or AfroDb entry name |
| `SMILES` | Canonical SMILES (RDKit-standardised) |
| `pIC50` | −log₁₀(IC50 in M); `NaN` for AfroDb (no activity data) |
| `source` | `"MMV"` or `"AfroDb"` |
| `MW` | Molecular weight |
| `LogP` | Calculated logP |
| `HBD` | H-bond donors |
| `HBA` | H-bond acceptors |
| `TPSA` | Topological polar surface area |
| `RotBonds` | Rotatable bond count |

---

## Licences and citation

| Dataset | Licence | Citation |
|---|---|---|
| MMV Malaria Box | CC BY 3.0 | Spangenberg *et al.*, *Malar. J.* 12, 175 (2013) |
| AfroDb | Academic use | Ntie-Kang *et al.*, *PLoS ONE* 8, e78085 (2013) |
| PubChem AID 2302 | Public domain | PubChem BioAssay AID 2302, NCBI |
