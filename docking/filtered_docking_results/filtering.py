


import pandas as pd
import os

# Paths to input files
docking_results_file = '/home/bio/TrmD_drug_design_project/docking/docking_results_sorted.csv'
smiles_file = '/home/bio/TrmD_drug_design_project/ligand_library/raw/smiles_mol_without_duplicates.csv'  # Replace with your SMILES file path
output_csv = '/home/bio/TrmD_drug_design_project/docking/filtered_docking_results/smiles_below_8.csv'

# Create output directory if needed
output_dir = os.path.dirname(output_csv)
os.makedirs(output_dir, exist_ok=True)

# Read docking results file
try:
    docking_results = pd.read_csv(docking_results_file, delimiter=",", encoding="utf-8-sig")
except Exception as e:
    print(f"Error reading docking results file: {e}")
    exit()

# Read SMILES file
try:
    smiles_data = pd.read_csv(smiles_file, delimiter=";", encoding="utf-8-sig")
except Exception as e:
    print(f"Error reading SMILES file: {e}")
    exit()

# Check for necessary columns
if 'Protein Code' not in docking_results.columns:
    raise ValueError("Docking results file must contain a 'score' column")
if 'Protein Code' not in docking_results.columns:
    raise ValueError("Docking results file must contain a 'ChEMBL ID' column")
if 'Smiles' not in smiles_data.columns:
    raise ValueError("SMILES file must contain a 'Smiles' column")
if 'ChEMBL ID' not in smiles_data.columns:
    raise ValueError("SMILES file must contain a 'ChEMBL ID' column")

# Filter docking results with score <= -7
filtered_docking_results = docking_results[docking_results['Docking Score'] <= -8]

# Match with SMILES data
filtered_smiles = smiles_data[smiles_data['ChEMBL ID'].isin(filtered_docking_results['Protein Code'])]

# Save the filtered SMILES to a new CSV file
try:
    filtered_smiles.to_csv(output_csv, index=False, encoding="utf-8-sig")
    print(f"Filtered SMILES saved to: {output_csv}")
except Exception as e:
    print(f"Error saving filtered SMILES: {e}")
