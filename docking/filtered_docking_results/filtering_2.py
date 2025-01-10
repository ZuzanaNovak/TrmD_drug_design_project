import pandas as pd


# Path to your input CSV file
input_csv = '/home/bio/TrmD_drug_design_project/docking/filtered_docking_results/smiles_below_8.csv'  # Replace with the actual path to your CSV file
output_csv = '/home/bio/TrmD_drug_design_project/docking/filtered_docking_results/just_smiles_below_8.csv'  # Replace with the desired output path for the extracted SMILES

# Read the CSV file
try:
    data = pd.read_csv(input_csv, delimiter=",", encoding="utf-8-sig")
except Exception as e:
    print(f"Error reading CSV file: {e}")
    exit()

# Check if the Smiles column exists
if 'Smiles' not in data.columns:
    print("No column named 'Smiles' found in the dataset.")
    exit()

# Extract the Smiles column
smiles_data = data['Smiles']

# Save the extracted Smiles to a new CSV file
try:
    smiles_data.to_csv(output_csv, index=False, header=True, encoding="utf-8-sig")
    print(f"SMILES strings successfully extracted and saved to: {output_csv}")
except Exception as e:
    print(f"Error saving SMILES to file: {e}")

