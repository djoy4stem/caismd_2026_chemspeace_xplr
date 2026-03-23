"""
prepare_pubchem_assay.py
========================
Downloads PubChem BioAssay AID 2302 (GSK Inhibition of P. falciparum Dd2,
whole-cell LDH assay) and creates a balanced 2 000-compound training sample
for the functional_group_profiling.ipynb Extension J workshop.

Assay summary (as of 2026):
  - Total tested compounds : 13 456
  - Active  (≥ 50% inhibition at 2 µM, IC50 < 2 µM): 7 921
  - Inactive                                          : 5 461
  - Active/total ratio ~ 58.9 %

Sample selected here (RANDOM_SEED = 42):
  - 1 180 Active  + 820 Inactive  = 2 000 total  (same 59/41 ratio)

Output
------
  data/pubchem_aid2302_2k.csv
  Columns:  CID, SMILES, activity ('Active' / 'Inactive'), outcome_int (1/2)

Usage
-----
  python src/prepare_pubchem_assay.py
  # or, from the notebook:
  # %run ../src/prepare_pubchem_assay.py
"""

import os
import time
import json
import random
import requests
import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger

RDLogger.DisableLog('rdApp.*')

RANDOM_SEED  = 42
TOTAL_SAMPLE = 2_000
AID          = 2302
BASE_URL     = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"

# Ratio mirrors the actual assay composition  (7921 active / 13382 tested ≈ 59.2 %)
ACTIVE_FRAC  = 7_921 / (7_921 + 5_461)        # ≈ 0.592
N_ACTIVE     = round(TOTAL_SAMPLE * ACTIVE_FRAC)   # 1 184
N_INACTIVE   = TOTAL_SAMPLE - N_ACTIVE             #   816

OUTPUT_DIR   = os.path.join(os.path.dirname(__file__), "..", "data")
OUTPUT_PATH  = os.path.join(OUTPUT_DIR, "pubchem_aid2302_2k.csv")

CHUNK_SIZE   = 200    # PUG-REST SMILES fetch chunk (keep < 200 to stay within URL limits)
SLEEP_S      = 0.35   # polite delay between API calls

# SSL verification — disable for conda envs that lack system CA bundle
# (PubChem is a US-NIH / NLM server; disabling verify is acceptable in a workshop context)
VERIFY_SSL   = False
import urllib3
urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)


# ─────────────────────────────────────────────────────────────────────────────
# Helpers
# ─────────────────────────────────────────────────────────────────────────────

def get_cids_for_activity(aid: int, activity_type: str) -> list[int]:
    """
    Return all CIDs for the given activity type ('active' or 'inactive')
    in a PubChem bioassay using the PUG-REST API.
    """
    url = f"{BASE_URL}/assay/aid/{aid}/cids/JSON?cids_type={activity_type}"
    r = requests.get(url, timeout=60, verify=VERIFY_SSL)
    r.raise_for_status()
    data = r.json()
    return data["InformationList"]["Information"][0]["CID"]


def fetch_smiles_batch(cids: list[int]) -> dict[int, str]:
    """
    Fetch isomeric SMILES for a batch of CIDs via PUG-REST.
    Returns {cid: smiles}.  CIDs with no SMILES are omitted.
    """
    cid_str = ",".join(str(c) for c in cids)
    url = f"{BASE_URL}/compound/cid/{cid_str}/property/IsomericSMILES,CanonicalSMILES/JSON"
    try:
        r = requests.get(url, timeout=60, verify=VERIFY_SSL)
        r.raise_for_status()
        props = r.json().get("PropertyTable", {}).get("Properties", [])
        # PubChem returns "SMILES" for isomeric, "ConnectivitySMILES" for canonical
        return {int(p["CID"]): p["SMILES"] for p in props if "SMILES" in p}
    except Exception as exc:
        print(f"  ⚠️  Batch fetch failed ({exc}) — skipping {len(cids)} CIDs")
        return {}


def fetch_smiles_for_cids(cids: list[int], desc: str = "") -> dict[int, str]:
    """
    Fetch SMILES for all CIDs in chunks, with a polite delay.
    """
    result = {}
    total  = len(cids)
    for i in range(0, total, CHUNK_SIZE):
        chunk = cids[i: i + CHUNK_SIZE]
        batch = fetch_smiles_batch(chunk)
        result.update(batch)
        pct = min(100, (i + len(chunk)) / total * 100)
        print(f"  {desc} {i + len(chunk)}/{total}  ({pct:.0f}%)  — {len(result)} SMILES so far",
              end="\r")
        time.sleep(SLEEP_S)
    print()  # newline after \r
    return result


def smiles_is_valid(smi: str) -> bool:
    """Return True if RDKit can parse the SMILES."""
    try:
        mol = Chem.MolFromSmiles(smi)
        return mol is not None
    except Exception:
        return False


# ─────────────────────────────────────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────────────────────────────────────

def main():
    rng = random.Random(RANDOM_SEED)

    # ── 1. Get all CIDs per activity class ──────────────────────────────────
    print(f"Fetching active CIDs for AID {AID}…")
    active_cids   = get_cids_for_activity(AID, "active")
    print(f"  Active CIDs   : {len(active_cids):,}")
    time.sleep(SLEEP_S)

    print(f"Fetching inactive CIDs for AID {AID}…")
    inactive_cids = get_cids_for_activity(AID, "inactive")
    print(f"  Inactive CIDs : {len(inactive_cids):,}")
    time.sleep(SLEEP_S)

    # ── 2. Random sample (oversample to account for invalid SMILES) ──────────
    OVERSAMPLE = 1.5   # fetch 50 % extra so we can filter invalids and still hit target
    sampled_active   = rng.sample(active_cids,   min(len(active_cids),   round(N_ACTIVE   * OVERSAMPLE)))
    sampled_inactive = rng.sample(inactive_cids, min(len(inactive_cids), round(N_INACTIVE * OVERSAMPLE)))

    print(f"\nSampled (with 50% overshoot): {len(sampled_active)} active, {len(sampled_inactive)} inactive")

    # ── 3. Fetch SMILES ──────────────────────────────────────────────────────
    print("\nFetching SMILES for active sample…")
    smi_active   = fetch_smiles_for_cids(sampled_active,   desc="Active  ")

    print("Fetching SMILES for inactive sample…")
    smi_inactive = fetch_smiles_for_cids(sampled_inactive, desc="Inactive")

    # ── 4. Build & filter DataFrame ─────────────────────────────────────────
    rows = []
    for cid, smi in smi_active.items():
        if smiles_is_valid(smi):
            rows.append({"CID": cid, "SMILES": smi, "activity": "Active",   "outcome_int": 1})
    for cid, smi in smi_inactive.items():
        if smiles_is_valid(smi):
            rows.append({"CID": cid, "SMILES": smi, "activity": "Inactive", "outcome_int": 2})

    df_all = pd.DataFrame(rows)
    print(f"\nValid molecules: {(df_all['activity']=='Active').sum()} active, "
          f"{(df_all['activity']=='Inactive').sum()} inactive")

    # Trim to exactly TOTAL_SAMPLE maintaining ratio
    df_active   = df_all[df_all['activity'] == 'Active'  ].sample(n=min(N_ACTIVE,   (df_all['activity']=='Active').sum()),   random_state=RANDOM_SEED)
    df_inactive = df_all[df_all['activity'] == 'Inactive'].sample(n=min(N_INACTIVE, (df_all['activity']=='Inactive').sum()), random_state=RANDOM_SEED)

    df_final = pd.concat([df_active, df_inactive]).sample(frac=1, random_state=RANDOM_SEED).reset_index(drop=True)

    # ── 5. Save ──────────────────────────────────────────────────────────────
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    df_final.to_csv(OUTPUT_PATH, index=False)

    act_n   = (df_final['activity'] == 'Active').sum()
    inact_n = (df_final['activity'] == 'Inactive').sum()
    print(f"\n✅  Saved {len(df_final)} compounds to {OUTPUT_PATH}")
    print(f"   Active  : {act_n} ({act_n/len(df_final)*100:.1f}%)")
    print(f"   Inactive: {inact_n} ({inact_n/len(df_final)*100:.1f}%)")
    print(f"   Source  : PubChem AID {AID} — GSK Inhibition of P. falciparum Dd2")


if __name__ == "__main__":
    main()
