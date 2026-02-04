
from rdkit import Chem
from rdkit.Chem import AllChem
from meeko import MoleculePreparation
from vina import Vina

# The Winning Recipe
smiles = "CC12CC3CC(C1)(C)CC(C3)C2CCCC1CC2(CCC1)OCC(=O)O2" # Memantine-Propyl-Lactone
name = "Halo_Winner_929"

print(f"💎 Re-Docking the Champion: {name}")

# 1. Build
mol = Chem.MolFromSmiles(smiles)
mol = Chem.AddHs(mol)
AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())
prep = MoleculePreparation()
prep.prepare(mol)
pdbqt = prep.write_pdbqt_string()

# 2. Dock
v = Vina(sf_name="vina")
v.set_receptor("halo_run/receptor.pdbqt")
v.compute_vina_maps(center=[3.43, -0.29, 17.69], box_size=[22, 22, 22])
v.set_ligand_from_string(pdbqt)
v.dock(exhaustiveness=32, n_poses=1) # High precision run
v.write_poses("halo_run/winner_929.pdbqt", n_poses=1, overwrite=True)

print("✅ Saved to halo_run/winner_929.pdbqt")

