import os
import pandas as pd

# Paths
ligand_dir = "/home/bio/TrmD_drug_design_project/REINVENT/results/reinvent1_output_pdbqt_v2"
receptor = "/home/bio/test3/pocket_v6.pdbqt"              # Receptor file
output_dir = "/home/bio/TrmD_drug_design_project/docking_v2/redocking_reinvent/redocking_output_v6"  # Output directory
csv_file = "/home/bio/TrmD_drug_design_project/docking_v2/docking_reinvent/docking_reinvent_sorted_v6.csv"  # File containing the docking scores

# Grid box parameters
center_x, center_y, center_z = 45.650, 1.391, 11.014
size_x, size_y, size_z = 20, 20, 20

# Exhaustiveness parameter for Vina
exhaustiveness = 32

# Read the docking scores and extract the top 500 molecules
scores_df = pd.read_csv(csv_file)
top_molecules = scores_df.nsmallest(500, "Docking Score")

# Ensure the output directory exists
os.makedirs(output_dir, exist_ok=True)

# Process each top molecule
for _, row in top_molecules.iterrows():
    ligand_file =  row["Protein Code"]  
    ligand_path = os.path.join(ligand_dir, ligand_file)

    if os.path.exists(ligand_path):  # Check if the ligand file exists
        # Generate unique output names
        ligand_name = os.path.splitext(ligand_file)[0]
        config_file = os.path.join(output_dir, f"{ligand_name}_config.txt")
        output_file = os.path.join(output_dir, f"{ligand_name}_output.pdbqt")

        # Create a config file for this ligand
        config_content = f"""
receptor = {receptor}
ligand = {ligand_path}
center_x = {center_x}
center_y = {center_y}
center_z = {center_z}
size_x = {size_x}
size_y = {size_y}
size_z = {size_z}
exhaustiveness = {exhaustiveness}
out = {output_file}
"""
        with open(config_file, "w") as f:
            f.write(config_content)

        # Run AutoDock Vina
        os.system(f"vina --config {config_file}")
    else:
        print(f"Warning: Ligand file {ligand_file} not found.")
