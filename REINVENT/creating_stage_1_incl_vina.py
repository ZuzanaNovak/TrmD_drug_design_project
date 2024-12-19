import os
import shutil
import subprocess
import tempfile
from rdkit import Chem
from rdkit.Chem import AllChem
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

import reinvent
from reinvent.notebooks import load_tb_data, plot_scalars, get_image, create_mol_grid


#doplnit název wd
wd = "___"
shutil.rmtree(wd, ignore_errors=True)
os.makedirs(wd, exist_ok=True)
os.chdir(wd)

# Global configuration parameters
global_parameters = """
run_type = "reinforcement_learning"
device = "cpu"
tb_logdir = "tb_stage1"
json_out_config = "_stage1.json"
"""

# Define the prior and agent model paths
prior_filename = os.path.join(reinvent.__path__[0], "..", "priors", "reinvent.prior")
agent_filename = prior_filename

parameters = f"""
[parameters]

prior_file = "{prior_filename}"
agent_file = "{agent_filename}"
summary_csv_prefix = "stage1"

batch_size = 100

use_checkpoint = false
"""

# Reinforcement Learning strategy
learning_strategy = """
[learning_strategy]

type = "dap"
sigma = 128
rate = 0.0001
"""

# Custom Docking Scoring Function
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

# Add Docking Scoring Component to REINVENT Configuration
receptor_file = "path/to/receptor.pdbqt"
vina_config = "path/to/vina_config.txt"

stages = f"""
[[stage]]

max_score = 1.0
max_steps = 300

chkpt_file = 'stage1.chkpt'

[stage.scoring]
type = "geometric_mean"

[[stage.scoring.component]]
[stage.scoring.component.custom]
name = "DockingScore"
weight = 1.0

[stage.scoring.component.custom.endpoint]
type = "external"
script = "python vina_scoring.py"
args = ["{{smiles}}", "{receptor_file}", "{vina_config}"]

[[stage.scoring.component]]
[stage.scoring.component.QED]

[[stage.scoring.component.QED.endpoint]]
name = "QED"
weight = 0.6

[[stage.scoring.component]]
[stage.scoring.component.NumAtomStereoCenters]

[[stage.scoring.component.NumAtomStereoCenters.endpoint]]
name = "Stereo"
weight = 0.4

transform.type = "left_step"
transform.low = 0
"""

# Write the combined configuration to a TOML file
config = global_parameters + parameters + learning_strategy + stages
toml_config_filename = "stage1.toml"

with open(toml_config_filename, "w") as tf:
    tf.write(config)

print("Configuration file 'stage1.toml' generated successfully.")
