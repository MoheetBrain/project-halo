import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors  # <--- The missing tool

print("🏭 Starting 'Project Halo' Chemical Library Generator...")

# 1. Define the Building Blocks (The "Legos")

# A. The Anchors (Adamantane Variations)
anchors = [
    ("Adamantane-1", "C12CC3CC(C1)CC(C3)C2"), 
    ("Adamantane-2", "C1C2CC3CC1CC(C2)C3"),
    ("Memantine-Core", "CC12CC3CC(C1)(C)CC(C3)C2"), 
    ("Rimantadine-Core", "CC(C12CC3CC(C1)CC(C3)C2)")
]

# B. The Linkers (The Chain)
linkers = [
    ("Direct", ""),
    ("Methyl", "C"),
    ("Ethyl", "CC"),
    ("Propyl", "CCC"),
    ("Ether", "O"),
    ("ThioEther", "S"),
    ("Amide", "NC(=O)"),
    ("ReverseAmide", "C(=O)N")
]

# C. The Heads (The Signaling Trigger)
heads = [
    ("Hedione-Core", "C1CC(=O)C(C1)CC(=O)OC"),
    ("Jasmone-Core", "C1CC(=O)C(C1)CC=CC"),
    ("Lactone-Head", "C1CC2(CCC1)OCC(=O)O2"), 
    ("Cyclohexanone", "C1CCCCC1=O")
]

library = []

print(f"⚗️  Mixing {len(anchors)} Anchors x {len(linkers)} Linkers x {len(heads)} Heads...")

count = 0
for anchor_name, anchor_smi in anchors:
    for linker_name, linker_smi in linkers:
        for head_name, head_smi in heads:
            
            # Construct the candidate
            if linker_name == "Direct":
                candidate_smi = f"{anchor_smi}{head_smi}"
            else:
                candidate_smi = f"{anchor_smi}{linker_smi}{head_smi}"
            
            # Name it
            id_tag = f"Halo_Gen1_{count:04d}"
            full_name = f"{anchor_name}-{linker_name}-{head_name}"
            
            # Validate with RDKit
            mol = Chem.MolFromSmiles(candidate_smi)
            if mol:
                # Calculate basic properties using Descriptors
                mw = Descriptors.MolWt(mol)
                logp = Descriptors.MolLogP(mol) # <--- Fixed line
                
                # Lipinski Rule Check
                if mw < 600: 
                    library.append({
                        "id": id_tag,
                        "name": full_name,
                        "smiles": candidate_smi,
                        "mw": mw,
                        "logp": logp
                    })
                    count += 1

# 3. Export
df = pd.DataFrame(library)
filename = "halo_candidates.csv"
df.to_csv(filename, index=False)

print("-" * 40)
print(f"✅ Library Generation Complete!")
print(f"📦 Total Valid Candidates: {len(df)}")
print(f"📂 Saved to: {filename}")
print(f"🧪 Average MW: {df['mw'].mean():.2f}")
print("-" * 40)
print("👉 Next Step: Run the screen.")