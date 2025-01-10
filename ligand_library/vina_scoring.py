import os
import sys
import tempfile
import subprocess
from rdkit import Chem
from rdkit.Chem import AllChem

def prepare_molecule(smiles):
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, AllChem.ETKDG())
    mol_file = tempfile.NamedTemporaryFile(suffix=".pdbqt", delete=False).name
    Chem.MolToPDBFile(mol, mol_file.replace(".pdbqt", ".pdb"))
    
    # Convert to PDBQT using OpenBabel
    subprocess.run(f"obabel {mol_file.replace('.pdbqt', '.pdb')} -O {mol_file}", shell=True)
    
    return mol_file

def dock_molecule(mol_file, receptor_file, vina_config):
    docking_output = tempfile.NamedTemporaryFile(suffix=".pdbqt", delete=False).name
    docking_log = tempfile.NamedTemporaryFile(suffix=".log", delete=False).name

    command = f"vina --receptor {receptor_file} --ligand {mol_file} --config {vina_config} --out {docking_output} --log {docking_log}"
    subprocess.run(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    
    # Parse docking score from log file
    with open(docking_log, 'r') as f:
        for line in f:
            if "REMARK VINA RESULT" in line:
                score = float(line.split()[3])
                return score
    return None

if __name__ == "__main__":
    smiles = sys.argv[1]
    receptor_file = sys.argv[2]
    vina_config = sys.argv[3]

    mol_file = prepare_molecule(smiles)
    score = dock_molecule(mol_file, receptor_file, vina_config)

    if score is not None:
        print(score)
    else:
        print("0.0")
