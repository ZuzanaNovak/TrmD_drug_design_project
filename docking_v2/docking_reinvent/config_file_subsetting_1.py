import os

ligand_dir = "/home/bio/TrmD_drug_design_project/REINVENT/subsets/subset_1"  # Specify subset directory here
receptor = "/home/bio/test3/pocket_v6.pdbqt"
output_dir = "/home/bio/TrmD_drug_design_project/docking_v2/docking_reinvent/output_reinvent_docking_v6_subset_1"
center_x, center_y, center_z = 45.650, 1.391, 11.014
size_x, size_y, size_z = 20, 20, 20

# Create output directory if it doesn't exist
os.makedirs(output_dir, exist_ok=True)

for ligand_file in os.listdir(ligand_dir):
    if ligand_file.endswith(".pdbqt"):
        ligand_path = os.path.join(ligand_dir, ligand_file)

        # Generate unique output names
        ligand_name = os.path.splitext(ligand_file)[0]
        config_file = os.path.join(output_dir, f"{ligand_name}.txt")
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
        out = {output_file}
        """
        
        config_file = os.path.join(output_dir, f"{ligand_name}_config.txt")
        with open(config_file, "w") as f:
            f.write(config_content)
        
        # Run AutoDock Vina
        os.system(f"vina --config {config_file}")

print(f"Docking completed for ligands in {ligand_dir}.")
