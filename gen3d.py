from mordred import Calculator, descriptors
from rdkit import Chem
from rdkit.Chem import AllChem
from openbabel import pybel
from molio import read_mol_block, write_mol_file
import numpy as np
from scipy.spatial.transform import Rotation


def make_mol(br_smiles: str, standlize: bool = True):
    br_mol = Chem.MolFromSmiles(br_smiles)
    patt = Chem.MolFromSmiles("Br")
    repl = Chem.MolFromSmiles("N")
    matches = br_mol.GetSubstructMatches(patt)
    Nidx = matches[0][0]

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
    mol_block = mol.write("mol")

    mol = Chem.MolFromMolBlock(mol_block, removeHs=False)

    neighbors = mol.GetAtomWithIdx(Nidx).GetNeighbors()
    Cidx = -1
    for atom in neighbors:
        if atom.GetSymbol() != "H":
            Cidx = atom.GetIdx()


    if standlize:
        e, c, p, t = read_mol_block(mol_block)
        c = c - c[Nidx]
        pole = np.cross(c[Cidx], np.array([0, 0, 1]))
        pole = pole / np.linalg.norm(pole)
        costheta = np.dot(c[Cidx], np.array([0, 0, 1])) / np.linalg.norm(c[Cidx])
        theta = np.arccos(costheta)
        rotevec = theta * pole
        r = Rotation.from_rotvec(rotevec)
        c = r.apply(c)
        c = c+np.array([0,0,30])
        mol_block = write_mol_file(e, c, p, t)

    return Nidx, Cidx, mol_block  # 0-index

if __name__ =="__main__":
    a, b, m = make_mol("CC(C)c1ccccc1CCBr", standlize=True)
    print(a, b)
    # print(m)
    with open("te-t.mol", "w") as f:
        f.write(m)
