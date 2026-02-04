import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski, AllChem
from meeko import MoleculePreparation
from vina import Vina
import os

# --- CONFIG ---
# We will audit the annotated hits file
INPUT_FILE = "halo_top_hits_annotated.csv" 
RECEPTOR = "halo_run/receptor.pdbqt"
CENTER = [3.43, -0.29, 17.69]
BOX_SIZE = [22, 22, 22]

print("👮 Starting FULL Forensic Audit (Stress Test)...")

# 1. Load Data
try:
    df = pd.read_csv(INPUT_FILE)
    print(f"📂 Loading {len(df)} candidates from: {INPUT_FILE}")
except FileNotFoundError:
    print(f"❌ Error: {INPUT_FILE} not found. Make sure you ran the analysis step!")
    exit()

print("-" * 80)
print(f"{'Name':<40} | {'Orig':<7} | {'Stress':<7} | {'Status'}")
print("-" * 80)

audit_results = []

for index, row in df.iterrows():
    name = row['name']
    smiles = row['smiles']
    orig_score = row['affinity']
    
    # 2. Validation
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        continue

    # 3. Prep Ligand (3D)
    mol_h = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol_h, AllChem.ETKDGv3())
    prep = MoleculePreparation()
    prep.prepare(mol_h)
    pdbqt_orig = prep.write_pdbqt_string()

    try:
        # 4. VINA STRESS TEST (Engine 1)
        # Crank exhaustiveness up to 32 (Standard is 8)
        v = Vina(sf_name='vina')
        v.set_receptor(RECEPTOR)
        v.compute_vina_maps(center=CENTER, box_size=BOX_SIZE)
        v.set_ligand_from_string(pdbqt_orig)
        
        v.dock(exhaustiveness=32, n_poses=1)
        stress_score = v.score()[0]
        
        # 5. Status Logic
        delta = stress_score - orig_score
        
        if stress_score < -9.5: 
            status = "💎 ELITE"
        elif stress_score < -9.0: 
            status = "🏆 Excellent"
        elif delta < -0.3: 
            status = "✅ Improved"
        elif delta > 0.5:  
            status = "⚠️ Unstable"
        else:              
            status = "🔹 Stable"

        print(f"{name:<40} | {orig_score:<7.2f} | {stress_score:<7.2f} | {status}")

        audit_results.append({
            "name": name,
            "original_score": orig_score,
            "stress_score": stress_score,
            "status": status
        })

    except Exception as e:
        print(f"{name:<40} | {orig_score:<7} | {'ERR':<7} | ❌ Error")

# Export
final_df = pd.DataFrame(audit_results)
final_df.to_csv("halo_audit_final.csv", index=False)
print("-" * 80)
print(f"✅ Audit Complete. Saved to 'halo_audit_final.csv'")