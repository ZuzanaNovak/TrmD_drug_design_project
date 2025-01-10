import os
import csv

def filter_and_save_vina_results(input_dir, output_csv):
    scores_and_files = []

    # Iterate through the files in the input directory
    for file_name in os.listdir(input_dir):
        if file_name.endswith('output.pdbqt'):
            file_path = os.path.join(input_dir, file_name)

            # Extract docking score from the file
            with open(file_path, 'r') as file:
                for line in file:
                    if "REMARK VINA RESULT:" in line:
                        score = float(line.split(":")[1].split()[0])

                        # Extract protein code from the file name (assuming it's part of the name)
                        protein_code = file_name.split('_')[0]
                        scores_and_files.append((protein_code, score))
                        break

    # Sort by docking score in ascending order
    scores_and_files.sort(key=lambda x: x[1])

    # Write all results to a CSV file
    with open(output_csv, 'w', newline='') as csvfile:
        csvwriter = csv.writer(csvfile)
        csvwriter.writerow(['Protein Code', 'Docking Score'])  # Write header

        for protein_code, score in scores_and_files:
            csvwriter.writerow([protein_code, score])

# Define input and output paths
input_directory = "/home/bio/TrmD_drug_design_project/docking/docking_results_first"
output_csv_file = "/home/bio/TrmD_drug_design_project/docking/docking_results_sorted.csv"

# Call the function
filter_and_save_vina_results(input_directory, output_csv_file)
