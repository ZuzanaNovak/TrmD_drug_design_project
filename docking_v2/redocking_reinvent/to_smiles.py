import os
import csv
import subprocess

# Directory containing PDBQT files
input_directory = "redocking_output_v6"  # Replace with your directory
output_csv = "output_smiles1.csv"

# Open the output CSV file for writing
with open(output_csv, "w", newline="") as csvfile:
    csv_writer = csv.writer(csvfile, delimiter=";")
    csv_writer.writerow(["Name", "SMILES"])  # Write the header

    # Loop through all files in the input directory
    for filename in os.listdir(input_directory):
        if filename.endswith("_output.pdbqt"):
            # Extract the name without "_output.pdbqt"
            name = filename.split("_output.pdbqt")[0]

            # Full path to the PDBQT file
            pdbqt_path = os.path.join(input_directory, filename)

            # Use Open Babel to convert only the first conformer to SMILES
            try:
                result = subprocess.run(
                    ["obabel", pdbqt_path, "-osmi", "--stop", "1"],
                    capture_output=True,
                    text=True,
                    check=True,
                )

                # Get the SMILES string from the output
                smiles = result.stdout.strip()

                # Write to the CSV file
                if smiles:
                    csv_writer.writerow([name, smiles])
                else:
                    print(f"No SMILES generated for file: {filename}")

            except subprocess.CalledProcessError as e:
                print(f"Error processing file {filename}: {e}")

print(f"Conversion complete. SMILES saved in {output_csv}")
