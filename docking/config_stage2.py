import pandas as pd
import os

# Path to your CSV file
input_csv = "/home/bio/TrmD_drug_design_project/REINVENT/results" 
output_dir = "/home/bio/TrmD_drug_design_project/docking/ligands_rev"  
os.makedirs(output_dir, exist_ok=True) 

# Read the CSV file with proper delimiter and encoding
try:
    data = pd.read_csv(input_csv, delimiter=",", encoding="utf-8-sig")  # Adjust delimiter if needed
except Exception as e:
    print(f"Error reading CSV: {e}")
    exit()

# Check if 'Smiles' column exists
if 'Smiles' not in data.columns:
    raise ValueError("CSV file must have a 'Smiles' column")

# Initialize counters
total_molecules = len(data)
successful_conversions = 0
failed_conversions = []

# Iterate through the rows
for index, row in data.iterrows():
    smiles = row['Smiles']
    molecule_id = row.get('ChEMBL ID', f"mol_{index+1}")  # Default ID if not provided
    output_pdbqt = os.path.join(output_dir, f"{molecule_id}.pdbqt")
    
    # Run OpenBabel to convert SMILES to PDBQT
    result = os.system(f'obabel -:"{smiles}" -O {output_pdbqt} --gen3d -p')
    
    # Check if the conversion was successful
    if result == 0:  # OpenBabel returns 0 for success
        successful_conversions += 1
    else:
        failed_conversions.append(molecule_id)

# Summary
print(f"Total molecules in CSV: {total_molecules}")
print(f"Successfully converted molecules: {successful_conversions}")
print(f"Failed conversions: {len(failed_conversions)}")

if failed_conversions:
    print("Failed molecule IDs:")
    for mol in failed_conversions:
        print(mol)
