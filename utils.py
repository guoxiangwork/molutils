from rdkit import Chem
from rdkit.Chem import AllChem
from openbabel import openbabel, pybel



def smi2mol(smiles):
    random_seed=42
    mol=AllChem.AddHs(Chem.MolFromSmiles(smiles))
    ncheck=AllChem.EmbedMolecule(mol,randomSeed=random_seed)
    if ncheck==-1:
        obConversion=openbabel.OBConversion()
        obConversion.SetInAndOutFormats("smi","mol")
        mol2=pybel.readstring('smi',smiles)
        mol2.make3D()
        #mol2.removeh()
        block=mol2.write('mol')
        mol3=Chem.MolFromMolBlock(block)
        mol=AllChem.AddHs(mol3,addCoords=True)
    '''
    obConversion=openbabel.OBConversion()
    obConversion.SetInAndOutFormats("smi","mol")
    mol2=pybel.readstring('smi',smiles)
    mol2.make3D()
    #mol2.removeh()
    block=mol2.write('mol')
    mol3=Chem.MolFromMolBlock(block)
    mol=AllChem.AddHs(mol3,addCoords=True)
    '''

    AllChem.MMFFOptimizeMolecule(mol)

    return Chem.MolToMolBlock(mol)