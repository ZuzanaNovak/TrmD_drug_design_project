import os
import re
import csv

# Directory containing the .pdbqt files
input_directories = [
    "/home/bio/TrmD_drug_design_project/docking_v2/docking_reinvent/output_reinvent_docking_v6_subset_1",
    "/home/bio/TrmD_drug_design_project/docking_v2/docking_reinvent/output_reinvent_docking_v6_subset_2",
    "/home/bio/TrmD_drug_design_project/docking_v2/docking_reinvent/output_reinvent_docking_v6_subset_3"
]
output_csv = "/home/bio/TrmD_drug_design_project/docking_v2/docking_reinvent/docking_reinvent_sorted_v6.csv"

# Regular expression to capture the docking score from .pdbqt files
score_pattern = re.compile(r"REMARK VINA RESULT:\s+([-+]?\d*\.\d+|\d+)")

# Store scores in a list
scores = []

# Process each file in the directory
for input_directory in input_directories:
    for filename in os.listdir(input_directory):
        if filename.endswith(".pdbqt"):
            file_path = os.path.join(input_directory, filename)
            with open(file_path, 'r') as file:
                content = file.read()
                match = score_pattern.search(content)
                if match:
                    score = float(match.group(1))
                    # Remove the suffix "_output.pdbqt" from the filename
                    molecule_name = filename.replace("_output", "")
                    scores.append((molecule_name, score))

# Sort the scores by the docking score (ascending order)
sorted_scores = sorted(scores, key=lambda x: x[1])

# Write the sorted scores to a CSV file
with open(output_csv, 'w', newline='') as csvfile:
    csv_writer = csv.writer(csvfile)
    csv_writer.writerow(["Protein Code", "Docking Score"])  # Header
    csv_writer.writerows(sorted_scores)

print(f"Docking scores saved to {output_csv}")
