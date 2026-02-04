#!/usr/bin/env python3
import argparse, os, sys, subprocess
import numpy as np
import pandas as pd
from Bio.PDB import MMCIFParser, PDBIO
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors
from meeko import MoleculePreparation, PDBQTWriterLegacy

# -----------------------------
# Robust Utilities
# -----------------------------

def run_cmd(cmd, check=True):
    """Runs system commands and captures errors visibly."""
    p = subprocess.run(cmd, capture_output=True, text=True)
    if check and p.returncode != 0:
        print(f"❌ SYSTEM ERROR during: {' '.join(cmd)}")
        print(f"STDERR: {p.stderr}")
        return None
    return p

def cif_to_pdb(cif_path, pdb_path):
    print(f"🔄 Converting CIF -> PDB...")
    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure("VN1R1", cif_path)
    io = PDBIO()
    io.set_structure(structure)
    io.save(pdb_path)

def pdb_to_pdbqt_receptor(pdb_path, pdbqt_path):
    print(f"🧪 Preparing Receptor (PDB -> PDBQT)...")
    run_cmd(["obabel", pdb_path, "-O", pdbqt_path, "-xr", "-xh", "--partialcharge", "gasteiger"])

def estimate_pocket_center_gpcr_heuristic(pdb_path):
    coords = []
    with open(pdb_path, "r") as f:
        for line in f:
            if line.startswith("ATOM"):
                coords.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
    coords = np.array(coords)
    z_thr = np.percentile(coords[:, 2], 70)
    top = coords[coords[:, 2] >= z_thr]
    return top.mean(axis=0)

def smiles_to_pdbqt_ligand(smiles, pdbqt_path):
    """Converts SMILES to a 3D PDBQT using RDKit and Meeko."""
    mol = Chem.MolFromSmiles(smiles)
    if not mol: return False, "Invalid SMILES"
    
    mol = Chem.AddHs(mol)
    if AllChem.EmbedMolecule(mol, AllChem.ETKDG()) != 0:
        return False, "3D Embedding Failed"
    
    AllChem.MMFFOptimizeMolecule(mol)
    
    preparator = MoleculePreparation()
    setups = preparator.prepare(mol)
    writer = PDBQTWriterLegacy()
    pdbqt_str, ok, err = writer.write_string(setups[0])
    
    if not ok: return False, f"Meeko Error: {err}"
    
    with open(pdbqt_path, "w") as f:
        f.write(pdbqt_str)
    return True, None

def parse_vina_score(log_txt):
    """Robustly pulls the top affinity score from a Vina log."""
    if not os.path.exists(log_txt): return None
    with open(log_txt, "r") as f:
        for line in f:
            parts = line.split()
            if len(parts) >= 2 and parts[0] == "1":
                try: return float(parts[1])
                except: continue
    return None

# -----------------------------
# Main Execution Loop
# -----------------------------

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--receptor-cif", required=True)
    ap.add_argument("--ligand-csv", required=True)
    ap.add_argument("--workdir", default="halo_run")
    ap.add_argument("--exhaustiveness", type=int, default=16)
    args = ap.parse_args()

    os.makedirs(args.workdir, exist_ok=True)
    
    # 1. Prep Receptor
    receptor_pdb = os.path.join(args.workdir, "receptor.pdb")
    receptor_pdbqt = os.path.join(args.workdir, "receptor.pdbqt")
    if not os.path.exists(receptor_pdbqt):
        cif_to_pdb(args.receptor_cif, receptor_pdb)
        pdb_to_pdbqt_receptor(receptor_pdb, receptor_pdbqt)

    center = estimate_pocket_center_gpcr_heuristic(receptor_pdb)
    print(f"🎯 Pocket Targeted at: {center}")

    # 2. Prep & Dock Ligands
    df = pd.read_csv(args.ligand_csv)
    results = []

    for idx, row in df.iterrows():
        name = row['name']
        smiles = row['smiles']
        print(f"\n🚀 Processing: {name}...")
        
        lig_pdbqt = os.path.join(args.workdir, f"{name}.pdbqt")
        out_pdbqt = os.path.join(args.workdir, f"{name}_out.pdbqt")
        log_txt = os.path.join(args.workdir, f"{name}.log")

        # Prep Ligand
        ok, err = smiles_to_pdbqt_ligand(smiles, lig_pdbqt)
        if not ok:
            print(f"   ⚠️ Skipping {name}: {err}")
            continue

        # Dock
        cmd = [
            "vina", "--receptor", receptor_pdbqt, "--ligand", lig_pdbqt,
            "--center_x", f"{center[0]:.3f}", "--center_y", f"{center[1]:.3f}", "--center_z", f"{center[2]:.3f}",
            "--size_x", "25", "--size_y", "25", "--size_z", "25",
            "--exhaustiveness", str(args.exhaustiveness), "--out", out_pdbqt, "--log", log_txt
        ]
        run_cmd(cmd, check=False)

        # Result
        score = parse_vina_score(log_txt)
        if score is not None:
            print(f"   ✅ Result: {score} kcal/mol")
            results.append({"name": name, "affinity": score})
        else:
            print(f"   ⚠️ Docking failed for {name}")

    # 3. Final Summary
    if results:
        res_df = pd.DataFrame(results).sort_values("affinity")
        res_df.to_csv(os.path.join(args.workdir, "final_results.csv"), index=False)
        print(f"\n🏆 WINNER: {res_df.iloc[0]['name']} at {res_df.iloc[0]['affinity']} kcal/mol")

if __name__ == "__main__":
    main()