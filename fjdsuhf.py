from atom import Molecule,Bond
from rdkit import Chem
from rdkit.Chem import AllChem
from openbabel import openbabel,pybel
import numpy as np

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

mb = smi2mol('BrCCCCC')
print(mb,file=open('mol1.mol','w'))
basismol_block='OpenBabel09052308583D\nui\n\n 61 62  0  0  0  0  0  0  0  0999 V2000\n    0.1451   -0.2331   -0.4365 N   0  0  0  0  0  0  0  0  0  0  0  0\n    1.8094   -0.3421    0.0197 P   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.6534    0.6976    0.7692 P   0  0  0  0  0  0  0  0  0  0  0  0\n    2.2193   -2.1046   -0.1139 C   0  0  0  0  0  0  0  0  0  0  0  0\n    1.2274   -3.0709    0.1011 C   0  0  0  0  0  0  0  0  0  0  0  0\n    3.5524   -2.5065   -0.2789 C   0  0  0  0  0  0  0  0  0  0  0  0\n    1.5676   -4.4206    0.1296 C   0  0  0  0  0  0  0  0  0  0  0  0\n    3.8846   -3.8568   -0.2411 C   0  0  0  0  0  0  0  0  0  0  0  0\n    2.8933   -4.8150   -0.0380 C   0  0  0  0  0  0  0  0  0  0  0  0\n    0.7944   -5.1715    0.2825 H   0  0  0  0  0  0  0  0  0  0  0  0\n    4.9203   -4.1619   -0.3757 H   0  0  0  0  0  0  0  0  0  0  0  0\n    3.1548   -5.8707   -0.0114 H   0  0  0  0  0  0  0  0  0  0  0  0\n    2.7679    0.4667   -1.2961 C   0  0  0  0  0  0  0  0  0  0  0  0\n    3.2077    1.7817   -1.0895 C   0  0  0  0  0  0  0  0  0  0  0  0\n    3.0091   -0.1489   -2.5335 C   0  0  0  0  0  0  0  0  0  0  0  0\n    3.8653    2.4737   -2.1017 C   0  0  0  0  0  0  0  0  0  0  0  0\n    3.6681    0.5447   -3.5432 C   0  0  0  0  0  0  0  0  0  0  0  0\n    4.0940    1.8550   -3.3288 C   0  0  0  0  0  0  0  0  0  0  0  0\n    4.2067    3.4926   -1.9302 H   0  0  0  0  0  0  0  0  0  0  0  0\n    3.8500    0.0614   -4.5008 H   0  0  0  0  0  0  0  0  0  0  0  0\n    4.6112    2.3933   -4.1204 H   0  0  0  0  0  0  0  0  0  0  0  0\n   -2.2813   -0.0227    1.1198 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -2.3010   -1.2127    1.8627 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -3.4918    0.5646    0.7327 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -3.5097   -1.8104    2.2015 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -4.6989   -0.0280    1.0888 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -4.7111   -1.2135    1.8200 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -3.5162   -2.7337    2.7782 H   0  0  0  0  0  0  0  0  0  0  0  0\n   -5.6349    0.4319    0.7773 H   0  0  0  0  0  0  0  0  0  0  0  0\n   -5.6580   -1.6721    2.0978 H   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.9238    2.3709    0.1206 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.1662    2.8272   -0.9665 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -1.7078    3.2894    0.8359 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.2067    4.1665   -1.3433 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -1.7481    4.6268    0.4537 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.9984    5.0683   -0.6354 C   0  0  0  0  0  0  0  0  0  0  0  0\n    0.3848    4.5042   -2.1925 H   0  0  0  0  0  0  0  0  0  0  0  0\n   -2.3654    5.3279    1.0117 H   0  0  0  0  0  0  0  0  0  0  0  0\n   -1.0286    6.1154   -0.9287 H   0  0  0  0  0  0  0  0  0  0  0  0\n    3.0260    2.2639   -0.1271 H   0  0  0  0  0  0  0  0  0  0  0  0\n    4.3305   -1.7616   -0.4491 H   0  0  0  0  0  0  0  0  0  0  0  0\n   -1.3594   -1.6711    2.1733 H   0  0  0  0  0  0  0  0  0  0  0  0\n    0.4689    2.1321   -1.5156 H   0  0  0  0  0  0  0  0  0  0  0  0\n   -2.2920    2.9567    1.6950 H   0  0  0  0  0  0  0  0  0  0  0  0\n   -3.4907    1.4696    0.1272 H   0  0  0  0  0  0  0  0  0  0  0  0\n    0.1889   -2.7583    0.2144 H   0  0  0  0  0  0  0  0  0  0  0  0\n    2.6847   -1.1762   -2.7014 H   0  0  0  0  0  0  0  0  0  0  0  0\n    1.3403    0.6420    2.2507 Cr  0  0  0  0  0  0  0  0  0  0  0  0\n    2.5448    1.5604    4.2365 C   0  0  0  0  0  0  0  0  0  0  0  0\n    1.8528    2.5024    3.2617 C   0  0  0  0  0  3  0  0  0  0  0  0\n    3.0727    2.1441    5.0092 H   0  0  0  0  0  0  0  0  0  0  0  0\n    3.3359    0.9780    3.7274 H   0  0  0  0  0  0  0  0  0  0  0  0\n    0.9864    2.9239    3.7271 H   0  0  0  0  0  0  0  0  0  0  0  0\n    2.5260    3.2867    2.9848 H   0  0  0  0  0  0  0  0  0  0  0  0\n    1.0034   -0.5366    4.0982 C   0  0  0  0  0  3  0  0  0  0  0  0\n    1.8165   -1.0748    3.5815 H   0  0  0  0  0  0  0  0  0  0  0  0\n    0.5718   -1.2679    4.7934 H   0  0  0  0  0  0  0  0  0  0  0  0\n    1.6063    0.5862    4.9411 C   0  0  0  0  0  0  0  0  0  0  0  0\n    2.1578    0.1217    5.7701 H   0  0  0  0  0  0  0  0  0  0  0  0\n    0.7889    1.1568    5.4130 H   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.1496   -0.4391   -1.3697 H   0  0  0  0  0  0  0  0  0  0  0  0\n  1  2  1  0  0  0  0\n  1  3  1  0  0  0  0\n  1 61  1  0  0  0  0\n  2  4  1  0  0  0  0\n  2 13  1  0  0  0  0\n  3 22  1  0  0  0  0\n  3 31  1  0  0  0  0\n  4  5  2  0  0  0  0\n  4  6  1  0  0  0  0\n  5  7  1  0  0  0  0\n  5 46  1  0  0  0  0\n  6  8  2  0  0  0  0\n  6 41  1  0  0  0  0\n  7  9  2  0  0  0  0\n  7 10  1  0  0  0  0\n  8  9  1  0  0  0  0\n  8 11  1  0  0  0  0\n  9 12  1  0  0  0  0\n 13 14  2  0  0  0  0\n 13 15  1  0  0  0  0\n 14 16  1  0  0  0  0\n 14 40  1  0  0  0  0\n 15 17  2  0  0  0  0\n 15 47  1  0  0  0  0\n 16 18  2  0  0  0  0\n 16 19  1  0  0  0  0\n 17 18  1  0  0  0  0\n 17 20  1  0  0  0  0\n 18 21  1  0  0  0  0\n 22 23  2  0  0  0  0\n 22 24  1  0  0  0  0\n 23 25  1  0  0  0  0\n 23 42  1  0  0  0  0\n 24 26  2  0  0  0  0\n 24 45  1  0  0  0  0\n 25 27  2  0  0  0  0\n 25 28  1  0  0  0  0\n 26 27  1  0  0  0  0\n 26 29  1  0  0  0  0\n 27 30  1  0  0  0  0\n 31 32  2  0  0  0  0\n 31 33  1  0  0  0  0\n 32 34  1  0  0  0  0\n 32 43  1  0  0  0  0\n 33 35  2  0  0  0  0\n 33 44  1  0  0  0  0\n 34 36  2  0  0  0  0\n 34 37  1  0  0  0  0\n 35 36  1  0  0  0  0\n 35 38  1  0  0  0  0\n 36 39  1  0  0  0  0\n 49 50  1  0  0  0  0\n 49 51  1  0  0  0  0\n 49 52  1  0  0  0  0\n 49 58  1  0  0  0  0\n 50 53  1  0  0  0  0\n 50 54  1  0  0  0  0\n 55 56  1  0  0  0  0\n 55 57  1  0  0  0  0\n 55 58  1  0  0  0  0\n 58 59  1  0  0  0  0\n 58 60  1  0  0  0  0\nM  END\n'
print(basismol_block,file=open('base.mol','w'))
mol1= Molecule.from_molblock(mb)
base = Molecule.from_molblock(basismol_block)

br_idx = [(idx,atom) for idx,atom in enumerate(mol1.atoms) if atom.symbol=='Br']
br_nei = [bond for bond in mol1.bonds if br_idx[0][1] in bond.atoms] 

mol1.translate(np.array([0,0,-20]))

mol_join = mol1+base
newBond = Bond(1,[mol_join.atoms[1],mol_join.atoms[17]],order=1)
mol_join.add_bonds([newBond])
print(mol_join,file=open('join.mol','w'))