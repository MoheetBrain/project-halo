from rdkit import Chem
from rdkit.Chem import GraphDescriptors, QED, RDConfig
import os, sys

# The Contender
name = "HaloV3"
smiles = "COC(=O)CC1C(CCC1=O)C23CC4CC(C3)CC(C4)C2"

print(f"⚗️  Analyzing Feasibility for: {name}")
print(f"🔗 SMILES: {smiles}")
print("-" * 40)

# 1. Sanitize Check (The "Fantasy" Filter)
mol = Chem.MolFromSmiles(smiles)
if mol is None:
    print("❌ CRITICAL FAIL: Molecule is chemically impossible (Invalid SMILES).")
    sys.exit()

try:
    Chem.SanitizeMol(mol)
    print("✅ Chemistry Check: PASSED (Bonds are valid)")
except:
    print("❌ CRITICAL FAIL: Valence errors detected (e.g. Carbon with 5 bonds).")
    sys.exit()

# 2. Complexity Analysis (The "Cost" Filter)
# BertzCT is a measure of topological complexity.
# < 500 = Easy
# 500 - 800 = Medium
# > 800 = Hard/Expensive
bertz = GraphDescriptors.BertzCT(mol)

# QED (Drug-likeness) - 0.0 (Bad) to 1.0 (Ideal)
qed = QED.qed(mol)

print(f"📊 Bertz Complexity: {bertz:.2f}")
print(f"💊 Drug-Likeness (QED): {qed:.2f}")

# 3. The Verdict
print("-" * 40)
if bertz < 500:
    print("🟢 VERDICT: EASY to Synthesize. (Likely cheap precursors)")
elif bertz < 800:
    print("🟡 VERDICT: MEDIUM Difficulty. (Requires standard lab setup)")
else:
    print("🔴 VERDICT: HARD / EXPENSIVE. (Complex ring systems detected)")

# Specific check for Adamantane linkage
if "C23CC4CC(C3)CC(C4)C2" in smiles:
    print("ℹ️  Note: Adamantane cage detected. This is a common, stable building block.")
    