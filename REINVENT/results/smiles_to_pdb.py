from rdkit import Chem
from rdkit.Chem import AllChem

def smiles_to_pdb(smiles, output_file):

    try:
        # Generate a molecule object from the SMILES string
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            print(f"Invalid SMILES: {smiles}")
            return
        # Add hydrogens
        mol = Chem.AddHs(mol)
        # Generate 3D coordinates
        AllChem.EmbedMolecule(mol, AllChem.ETKDG())
        AllChem.UFFOptimizeMolecule(mol)  # Optimize geometry
        # Write to a PDB file
        Chem.MolToPDBFile(mol, output_file)
        print(f"PDB file saved to: {output_file}")
    except Exception as e:
        print(f"Error during conversion: {e}")

# Example usage
if __name__ == "__main__":
    smiles = "Fc1cc(CCc2cccc(-c3cccnc3)c2)cc(C(F)(F)F)c1"  
    output_file = "output33.pdb"
    smiles_to_pdb(smiles, output_file)
