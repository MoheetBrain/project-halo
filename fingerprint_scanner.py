import math

def parse_pdbqt(filename, is_ligand=False):
    atoms = []
    with open(filename, 'r') as f:
        for line in f:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                # PDBQT columns:
                # 30-38 = X, 38-46 = Y, 46-54 = Z
                try:
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    
                    if is_ligand:
                        atoms.append((x, y, z))
                    else:
                        # Extract Residue Info for Receptor
                        # 17-20 = Residue Name, 22-26 = Residue ID
                        res_name = line[17:20].strip()
                        res_id = line[22:26].strip()
                        atoms.append({'x': x, 'y': y, 'z': z, 'res': f"{res_name}{res_id}"})
                except:
                    continue
            # Stop ligand parsing after the first model (best pose)
            if is_ligand and line.startswith("ENDMDL"):
                break
    return atoms

def calc_dist(a1, a2):
    return math.sqrt((a1[0]-a2['x'])**2 + (a1[1]-a2['y'])**2 + (a1[2]-a2['z'])**2)

print("🕵️‍♂️ Scanning Interaction Fingerprint...")

# 1. Load the Best Pose
ligand_atoms = parse_pdbqt("halo_run/HaloV3_seed3_out.pdbqt", is_ligand=True)
print(f"✅ Loaded Ligand: {len(ligand_atoms)} atoms")

# 2. Load the Receptor
receptor_atoms = parse_pdbqt("halo_run/receptor.pdbqt", is_ligand=False)
print(f"✅ Loaded Receptor: {len(receptor_atoms)} atoms")

# 3. Calculate Contacts (Cutoff = 4.0 Angstroms)
contacts = set()
cutoff = 4.0

print(f"🔍 Searching for contacts within {cutoff} Å...")

for l_atom in ligand_atoms:
    for r_atom in receptor_atoms:
        dist = calc_dist(l_atom, r_atom)
        if dist <= cutoff:
            contacts.add(r_atom['res'])

# 4. The Verdict
print("-" * 30)
print("🧬 INTERACTING RESIDUES:")
sorted_contacts = sorted(list(contacts))
print(", ".join(sorted_contacts))
print("-" * 30)

# Check for Key Hydrophobic Residues (Based on Homology)
# Note: These are hypothetical "hotspots" for V1Rs. 
# If we see lots of Leucine (LEU), Valine (VAL), Phenylalanine (PHE), it's good.
hydrophobic_hits = [r for r in contacts if "LEU" in r or "VAL" in r or "PHE" in r or "ILE" in r]
print(f"💪 Hydrophobic Anchors Hit: {len(hydrophobic_hits)}")