import pandas as pd
import glob
from rdkit import Chem
from rdkit.Chem import AllChem
from meeko import MoleculePreparation
from vina import Vina

# --- CONFIGURATION ---
RECEPTOR_FILE = "halo_run/receptor.pdbqt"
LIBRARY_FILE = "halo_candidates.csv"
RESULTS_FILE = "halo_screening_results.csv"

# The coordinates we found earlier (VN1R1 Pocket)
CENTER = [3.43, -0.29, 17.69]
BOX_SIZE = [22, 22, 22]  # Slightly tighter box for speed

print(f"🚀 Starting Halo-Dock Screen on {LIBRARY_FILE}...")

# 1. Load the Library
df = pd.read_csv(LIBRARY_FILE)
results = []

# 2. Setup Vina Engine
v = Vina(sf_name='vina')
v.set_receptor(RECEPTOR_FILE)
v.compute_vina_maps(center=CENTER, box_size=BOX_SIZE)

print(f"⚡ Vina Grid Computed. Processing {len(df)} candidates...")
print("-" * 60)
print(f"{'ID':<15} | {'Name':<35} | {'Affinity':<10}")
print("-" * 60)

# 3. The Docking Loop
for index, row in df.iterrows():
    name = row['name']
    smiles = row['smiles']
    cand_id = row['id']
    
    try:
        # A. Build 3D Molecule
        mol = Chem.MolFromSmiles(smiles)
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol, AllChem.ETKDGv3()) # Generate 3D coords
        
        # B. Prepare PDBQT (Meeko)
        preparator = MoleculePreparation()
        preparator.prepare(mol)
        pdbqt_string = preparator.write_pdbqt_string()
        
        # C. Dock It
        v.set_ligand_from_string(pdbqt_string)
        v.dock(exhaustiveness=8, n_poses=1) # Fast screening mode
        
        # D. Get Score
        score = v.score()[0] # Best energy
        
        print(f"{cand_id:<15} | {name:<35} | {score:.2f}")
        
        results.append({
            'id': cand_id,
            'name': name,
            'smiles': smiles,
            'affinity': score
        })
        
    except Exception as e:
        print(f"❌ Error on {name}: {e}")

# 4. Save & Sort
print("-" * 60)
results_df = pd.DataFrame(results)
# Sort by best score (lowest number is better)
results_df = results_df.sort_values(by='affinity', ascending=True)
results_df.to_csv(RESULTS_FILE, index=False)

best_cand = results_df.iloc[0]
print(f"🏆 WINNER: {best_cand['name']} ({best_cand['affinity']:.2f} kcal/mol)")
print(f"💾 Saved full results to {RESULTS_FILE}")