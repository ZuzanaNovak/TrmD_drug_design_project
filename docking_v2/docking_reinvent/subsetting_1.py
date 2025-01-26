import os
import shutil

ligand_dir = "/home/bio/TrmD_drug_design_project/REINVENT/results/reinvent1_output_pdbqt_v2"
output_base = "/home/bio/TrmD_drug_design_project/REINVENT/subsets"
os.makedirs(output_base, exist_ok=True)

# Create subset directories
subset_dirs = [
    os.path.join(output_base, f"subset_{i}") for i in range(1, 4)
]
for subset_dir in subset_dirs:
    os.makedirs(subset_dir, exist_ok=True)

# Get list of all ligands
ligand_files = [f for f in os.listdir(ligand_dir) if f.endswith(".pdbqt")]

# Split into three roughly equal subsets
for i, ligand_file in enumerate(ligand_files):
    subset_index = i % 3  # Round-robin distribution
    shutil.copy(os.path.join(ligand_dir, ligand_file), subset_dirs[subset_index])

print(f"Split {len(ligand_files)} ligands into 3 subsets.")
