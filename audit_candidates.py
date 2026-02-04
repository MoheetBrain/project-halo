import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski, AllChem
from meeko import MoleculePreparation
from vina import Vina
import os

# --- CONFIG ---
# We will audit the annotated hits file we created earlier
INPUT_FILE = "halo_top_hits_annotated.csv" 
RECEPTOR = "halo_run/receptor.pdbqt"
CENTER = [3.43, -0.29, 17.69]
BOX_SIZE = [22, 22, 22]

print("👮 Starting FULL Forensic Audit (Stress Test)...")
print(f"📂 Loading candidates from: {INPUT_FILE}")

try:
    df = pd.read_csv(INPUT_FILE)
except FileNotFoundError:
    print(f"❌ Error: {INPUT_FILE} not found.")
    exit()

print(f"⚡ Queue size: {len(df)} candidates")
print("-" * 80)
print(f"{'Name':<40} | {'Orig':<7} | {'Stress':<7} | {'Vinardo':<7} | {'Status'}")
print("-" * 80)

audit_results = []

for index, row in df.iterrows():
    name = row['name']
    smiles = row['smiles']
    orig_score = row['affinity']
    
    # 1. Validation (Skip bad molecules)
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        print(f"{name:<40} | {orig_score:<7} | {'---':<7} | {'---':<7} | ❌ Invalid SMILES")
        continue

    # 2. Prep Ligand
    mol_h = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol_h, AllChem.ETKDGv3())
    prep = MoleculePreparation()
    prep.prepare(mol_h)
    pdbqt_orig = prep.write_pdbqt_string()

    try:
        # 3. VINA STRESS TEST (Engine 1)
        # We crank exhaustiveness up to 32 (Standard is 8)
        v = Vina(sf_name='vina')
        v.set_receptor(RECEPTOR)
        v.compute_vina_maps(center=CENTER, box_size=BOX_SIZE)
        v.set_ligand_from_string(pdbqt_orig)
        
        v.dock(exhaustiveness=32, n_poses=1)
        stress_score = v.score()[0]
        
        # Get the pose string for Vinardo
        docked_pose = v.poses(n_poses=1)

        # 4. VINARDO CHECK (Engine 2)
        # Fix: Manually strip MODEL/ENDMDL tags to prevent crash
        clean_pose = []
        for line in docked_pose.splitlines():
            if not line.startswith("MODEL") and not line.startswith("ENDMDL"):
                clean_pose.append(line)
        clean_pose_str = "\n".join(clean_pose)

        v_ortho = Vina(sf_name='vinardo')
        v_ortho.set_receptor(RECEPTOR)
        v_ortho.compute_vina_maps(center=CENTER, box_size=BOX_SIZE)
        v_ortho.set_ligand_from_string(clean_pose_str)
        vinardo_score = v_ortho.score()[0]
        
        # Status Logic
        delta = stress_score - orig_score
        if stress_score < -9.0: status = "💎 ELITE"
        elif delta < -0.5: status = "✅ Improved"
        elif delta > 0.5:  status = "⚠️ Unstable"
        else:              status = "🔹 Stable"

        print(f"{name:<40} | {orig_score:<7.2f} | {stress_score:<7.2f} | {vinardo_score:<7.2f} | {status}")

        audit_results.append({
            "name": name,
            "smiles": smiles,
            "original_score": orig_score,
            "stress_score": stress_score,
            "vinardo_score": vinardo_score,
            "status": status
        })

    except Exception as e:
        print(f"{name:<40} | {orig_score:<7} | {'ERR':<7} | {'ERR':<7} | ❌ {str(e)[:20]}")

# Export final audited list
final_df = pd.DataFrame(audit_results)
final_df = final_df.sort_values(by="stress_score", ascending=True)
final_df.to_csv("halo_audit_final.csv", index=False)
print("-" * 80)
print(f"✅ Audit Complete. Saved {len(final_df)} verified hits to 'halo_audit_final.csv'")