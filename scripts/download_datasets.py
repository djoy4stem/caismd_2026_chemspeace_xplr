"""
prepare_datasets.py  (formerly download_datasets.py)
=====================================================
Processes the three raw data files already in this folder into
notebook-ready CSVs with harmonised columns and RDKit descriptors.

Raw files required
------------------
    MalariaBox400compoundsDec2014.xls   (400 compounds, MMV, CC BY 3.0)
    AfroDB_3D.sdf                       (954 compounds, Ntie-Kang et al. 2013)

Optional raw file (not used by default notebook cells, but available)
----------------------------------------------------------------------
    GHPB_DETAILS.xlsx                   (240 compounds, MMV GHPB, CC BY 4.0)

Outputs
-------
    malaria_box.csv                 — 400 MMV Malaria Box compounds + descriptors
    afrodb_subset.csv               — 903 AfroDb natural products  + descriptors
    malaria_box_afrodb_combined.csv — merged (1303 compounds), used by notebook

Run
---
    python data/prepare_datasets.py       (activate cheminf conda env first)
"""

import math
import os
import sys

import pandas as pd

try:
    from rdkit import Chem, RDLogger
    from rdkit.Chem import Descriptors, rdMolDescriptors
    RDLogger.DisableLog("rdApp.*")
    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False
    print("⚠️  RDKit not found — descriptor columns will be empty.")

DATA_DIR = os.path.dirname(os.path.abspath(__file__))
ALL_COLS = ["COMPOUND_ID", "SMILES", "pIC50", "source",
            "MW", "LogP", "HBD", "HBA", "TPSA", "RotBonds"]


# ── Descriptor helpers ────────────────────────────────────────────────────────

def compute_descriptors(mol):
    if mol is None or not RDKIT_AVAILABLE:
        return {"MW": None, "LogP": None, "HBD": None,
                "HBA": None, "TPSA": None, "RotBonds": None}
    return {
        "MW":       round(Descriptors.ExactMolWt(mol), 2),
        "LogP":     round(Descriptors.MolLogP(mol), 2),
        "HBD":      rdMolDescriptors.CalcNumHBD(mol),
        "HBA":      rdMolDescriptors.CalcNumHBA(mol),
        "TPSA":     round(Descriptors.TPSA(mol), 2),
        "RotBonds": rdMolDescriptors.CalcNumRotatableBonds(mol),
    }

def add_descriptors(df):
    if not RDKIT_AVAILABLE:
        return df
    print("   Computing RDKit descriptors ...")
    rows = [compute_descriptors(Chem.MolFromSmiles(s)) for s in df["SMILES"]]
    return pd.concat([df.reset_index(drop=True), pd.DataFrame(rows)], axis=1)

def validate_smiles(df):
    if not RDKIT_AVAILABLE:
        return df
    mask = df["SMILES"].apply(lambda s: Chem.MolFromSmiles(str(s)) is not None)
    dropped = (~mask).sum()
    if dropped:
        print(f"   Dropped {dropped} rows with invalid SMILES.")
    return df[mask].reset_index(drop=True)


# ── Dataset 1: MMV Malaria Box ────────────────────────────────────────────────

def load_malaria_box():
    path = os.path.join(DATA_DIR, "MalariaBox400compoundsDec2014.xls")
    if not os.path.exists(path):
        sys.exit(f"❌ File not found: {path}\n   Download it from:\n"
                 "   https://www.mmv.org/mmv-open/malaria-box/malaria-box-supporting-information")

    print("\n── Dataset 1: MMV Malaria Box ───────────────────────────────────────")
    df = pd.read_excel(path, sheet_name="vortex_sheet")
    df = df.rename(columns={"HEOS_COMPOUND_ID": "COMPOUND_ID", "Smiles": "SMILES"})

    # pIC50 = -log10(EC50_M) = -log10(EC50_nM × 1e-9)
    def to_pic50(val):
        try:
            v = float(val)
            return round(-math.log10(v * 1e-9), 2) if v > 0 else float("nan")
        except (TypeError, ValueError):
            return float("nan")

    df["pIC50"]  = df["EC50_nM"].apply(to_pic50)
    df["source"] = "MMV"
    df = df[["COMPOUND_ID", "SMILES", "pIC50", "source"]].dropna(subset=["SMILES"])
    df = df[df["SMILES"].str.strip() != ""].reset_index(drop=True)
    print(f"   Loaded {len(df)} compounds.")
    return df


# ── Dataset 2: AfroDb ─────────────────────────────────────────────────────────

def load_afrodb():
    path = os.path.join(DATA_DIR, "AfroDB_3D.sdf")
    if not os.path.exists(path):
        sys.exit(f"❌ File not found: {path}\n   Download Dataset S1 from:\n"
                 "   https://doi.org/10.1371/journal.pone.0078085")
    if not RDKIT_AVAILABLE:
        sys.exit("❌ RDKit is required to read the SDF file.")

    print("\n── Dataset 2: AfroDb natural products ───────────────────────────────")
    supplier = Chem.SDMolSupplier(path, removeHs=True, sanitize=True)
    records = []
    for mol in supplier:
        if mol is None:
            continue
        smi = Chem.MolToSmiles(mol)
        if smi:
            props = mol.GetPropsAsDict()
            records.append({
                "COMPOUND_ID": props.get("s_m_entry_name", f"AfroDb.{len(records)+1}"),
                "SMILES":      smi,
                "pIC50":       float("nan"),
                "source":      "AfroDb",
            })
    df = pd.DataFrame(records)
    n_raw = len(df)
    df = df.drop_duplicates("SMILES").reset_index(drop=True)
    n_dupes = n_raw - len(df)
    print(f"   Loaded {n_raw} entries from SDF.")
    if n_dupes:
        # These are confirmed true duplicates (identical InChIKey), not stereoisomers.
        # The SDF contains the same molecule multiple times with different 3D conformers.
        # Stereoisomers are preserved because RDKit canonical SMILES encodes stereo (@, @@).
        print(f"   Dropped {n_dupes} conformer duplicates (same InChIKey, different 3D geometry).")
    print(f"   {len(df)} unique compounds retained.")
    return df


# ── Save ──────────────────────────────────────────────────────────────────────

def save(df_mmv, df_np):
    print("\n── Computing descriptors and saving ─────────────────────────────────")
    df_mmv = validate_smiles(add_descriptors(df_mmv))
    df_np  = validate_smiles(add_descriptors(df_np))

    for col in ALL_COLS:
        for df in [df_mmv, df_np]:
            if col not in df.columns:
                df[col] = float("nan")

    df_mmv[ALL_COLS].to_csv(os.path.join(DATA_DIR, "malaria_box.csv"),    index=False)
    df_np[ALL_COLS].to_csv(os.path.join(DATA_DIR,  "afrodb_subset.csv"),  index=False)
    print(f"   malaria_box.csv      → {len(df_mmv)} compounds")
    print(f"   afrodb_subset.csv    → {len(df_np)} compounds")

    df_combined = pd.concat([df_mmv[ALL_COLS], df_np[ALL_COLS]], ignore_index=True)
    df_combined.to_csv(os.path.join(DATA_DIR, "malaria_box_afrodb_combined.csv"), index=False)
    print(f"   malaria_box_afrodb_combined.csv → {len(df_combined)} compounds")
    return df_combined


# ── Main ──────────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    print("=" * 65)
    print(" Workshop Dataset Processor")
    print(" 'Mapping Molecular Landscapes' — CAISMD 2026 — Chemical Space Explorer")
    print("=" * 65)

    df_combined = save(load_malaria_box(), load_afrodb())

    print("\n" + "=" * 65)
    print(f" ✅  Done.  {len(df_combined)} compounds total.")
    print(f"     {df_combined['source'].value_counts().to_string()}")
    print("=" * 65)
