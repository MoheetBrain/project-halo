import pandas as pd

print("📊 Analyzing Halo Screening Results...")

# 1. Load the raw results
try:
    df = pd.read_csv("halo_screening_results.csv")
except FileNotFoundError:
    print("❌ Error: 'halo_screening_results.csv' not found. Did the screen finish?")
    exit()

# 2. Parse the Name to get Components
def parse_name(row):
    # Heuristic parsing based on your naming convention
    core = "Unknown"
    head = "Unknown"
    linker = "Direct"
    
    name = row['name']
    
    if "Memantine" in name: core = "Memantine"
    elif "Rimantadine" in name: core = "Rimantadine"
    elif "Adamantane" in name: core = "Adamantane"
    
    if "Lactone" in name: head = "Lactone"
    elif "Hedione" in name: head = "Hedione"
    elif "Jasmone" in name: head = "Jasmone"
    elif "Cyclohexanone" in name: head = "Cyclohexanone"
    
    if "Propyl" in name: linker = "Propyl"
    elif "Ethyl" in name: linker = "Ethyl"
    elif "Methyl" in name: linker = "Methyl"
    elif "Amide" in name: linker = "Amide"
    elif "ReverseAmide" in name: linker = "RevAmide"
    elif "ThioEther" in name: linker = "ThioEther"
    elif "Ether" in name: linker = "Ether"
    
    return pd.Series([core, linker, head])

df[['Core', 'Linker', 'Head']] = df.apply(parse_name, axis=1)

# 3. Sort by Affinity (Best First)
df_sorted = df.sort_values(by='affinity', ascending=True)

# 4. Print the "Top 5 Synthesis Candidates"
print("-" * 60)
print(f"{'Rank':<5} | {'Name':<40} | {'Score':<8} | {'Core':<12}")
print("-" * 60)

for i in range(5):
    row = df_sorted.iloc[i]
    print(f"{i+1:<5} | {row['name']:<40} | {row['affinity']:<8.2f} | {row['Core']:<12}")

print("-" * 60)

# Export for records
df_sorted.head(20).to_csv("halo_top_hits_annotated.csv", index=False)
print("✅ Saved top 20 candidates to 'halo_top_hits_annotated.csv'")
