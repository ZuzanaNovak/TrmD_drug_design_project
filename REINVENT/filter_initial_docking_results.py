import os
import shutil

scores_and_files = []
def filter_and_copy_lowest_vina_results(input_dir, output_dir, top_n=3):

    for file_name in os.listdir(input_dir):
        if file_name.endswith('output.pdbqt'):
            file_path = os.path.join(input_dir, file_name)

            with open(file_path, 'r') as file:
                for line in file:
                    if "REMARK VINA RESULT:" in line:
                        score = float(line.split(":")[1].split()[0])
                        scores_and_files.append((score, file_path))
                        break

    scores_and_files.sort(key=lambda x: x[0])
    lowest_files = scores_and_files[:top_n]


    for _, file_path in lowest_files:
        shutil.copy(file_path, output_dir)

input_directory = "/home/bio/TrmD_drug_design_project/docking/docking_results_first"  
output_directory = "/home/bio/TrmD_drug_design_project/docking/docking_results_filtered"  

filter_and_copy_lowest_vina_results(input_directory, output_directory, top_n=3)
