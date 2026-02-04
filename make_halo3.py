import os
from rdkit import Chem
from rdkit.Chem import AllChem
from meeko import MoleculePreparation, PDBQTWriterLegacy

# The Adamantane-Hybrid Candidate
name = "HaloV3"
smiles = "COC(=O)CC1C(CCC1=O)C23CC4CC(C3)CC(C4)C2"

print(f"💎 Generating {name} from SMILES: {smiles}")

# 1. Generate 3D Structure
mol = Chem.MolFromSmiles(smiles)
mol = Chem.AddHs(mol)
AllChem.EmbedMolecule(mol, AllChem.ETKDG())
AllChem.MMFFOptimizeMolecule(mol)

# 2. Write to PDBQT using Meeko
preparator = MoleculePreparation()
setups = preparator.prepare(mol)
writer = PDBQTWriterLegacy()
pdbqt_string = writer.write_string(setups[0])[0]

out_file = f"halo_run/{name}.pdbqt"
with open(out_file, "w") as f:
    f.write(pdbqt_string)

print(f"✅ Success! Saved to {out_file}")
