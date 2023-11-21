from mordred import Calculator, descriptors
from rdkit import Chem
from rdkit.Chem import AllChem
from openbabel import pybel


def make_mol(br_smiles: str):
    br_mol = Chem.MolFromSmiles(br_smiles)
    patt = Chem.MolFromSmiles("Br")
    repl = Chem.MolFromSmiles("N")
    matches = br_mol.GetSubstructMatches(patt)

    n_mol = AllChem.ReplaceSubstructs(br_mol, patt, repl)[0]
    n_mol.UpdatePropertyCache(strict=False)
    Chem.SanitizeMol(
        n_mol,
        Chem.SanitizeFlags.SANITIZE_FINDRADICALS
        | Chem.SanitizeFlags.SANITIZE_KEKULIZE
        | Chem.SanitizeFlags.SANITIZE_SETAROMATICITY
        | Chem.SanitizeFlags.SANITIZE_SETCONJUGATION
        | Chem.SanitizeFlags.SANITIZE_SETHYBRIDIZATION
        | Chem.SanitizeFlags.SANITIZE_SYMMRINGS,
        catchErrors=True,
    )
    n_mol_h = Chem.AddHs(n_mol)
    # AllChem.EmbedMolecule(n_mol_h)
    # AllChem.MMFFOptimizeMolecule(n_mol_h)
    n_mol_h_2d_molblock = Chem.MolToMolBlock(n_mol_h)
    mol = pybel.readstring("mol", n_mol_h_2d_molblock)
    mol.make3D(steps=200)

    return matches[0][0], mol.write("mol")

'''
a, m = make_mol("c1ccccc1CCBr")
mol = Chem.MolFromMolBlock(m,removeHs=False)
print(Chem.MolToMolBlock(mol))
df = calc.pandas(mols)
'''

