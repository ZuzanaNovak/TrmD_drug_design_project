import pandas as pd
import os

# Path to your CSV file
input_csv = "/Users/ladislavbucek/Desktop/UK/Mgr./1./Meet-U/TrmD_drug_design_project/ligand library/raw/smiles_mol_without_duplicates.csv"  # Replace with your CSV file path
output_dir = "output_pdbqt"  # Directory to save PDBQT files
os.makedirs(output_dir, exist_ok=True)  # Create the output directory if it doesn't exist

# Read the CSV file with proper delimiter and encoding
try:
    data = pd.read_csv(input_csv, delimiter=";", encoding="utf-8-sig")
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
optimization_steps = 1000  # Adjust based on molecule complexity
force_field = "MMFF94"  # or "UFF" if necessary

# Iterate through the rows
for index, row in data.iterrows():
    smiles = row['Smiles']
    molecule_id = row.get('ChEMBL ID', f"mol_{index+1}")  # Default ID if not provided
    
    # Define file paths
    initial_pdb = os.path.join(output_dir, f"{molecule_id}_initial.pdb")
    optimized_pdb = os.path.join(output_dir, f"{molecule_id}_optimized.pdb")
    pdbqt_file = os.path.join(output_dir, f"{molecule_id}.pdbqt")

    # Step 1: Generate 3D coordinates with hydrogens
    os.system(f'obabel -:"{smiles}" -O {initial_pdb} --gen3d -p')


    # Step 2: Optimize geometry
    os.system(f'obminimize -ff {force_field} -n {optimization_steps} {initial_pdb} -O {optimized_pdb}')

    # Step 3: Convert optimized PDB to PDBQT
    result = os.system(f'obabel {optimized_pdb} -O {pdbqt_file}')
    
    os.remove(initial_pdb)
    #os.remove(optimized_pdb)

    print(f"Processed {molecule_id}: {pdbqt_file}")

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
