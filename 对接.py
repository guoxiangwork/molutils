from atom import Molecule, Bond
from rdkit import Chem
from rdkit.Chem import AllChem
from openbabel import openbabel, pybel
import numpy as np
from utils import smi2mol

# smi = open('remove_duplication-562.smi','r').read().split()[7]
# mb = smi2mol(smi)
# print(smi)
# print(mb, file=open(f"7.mol", "w"))
# grop = Molecule.from_molblock(mb)
# print(grop, file=open(f"g-7.mol", "w"))
# #grop2 = grop.rotate(np.array([0,0,0]),np.array([0,0,1]),angle=np.pi/3)
# #print(grop2, file=open(f"g2-537.mol", "w"))
# basismol_block = open('base.mol','r').read()
# base = Molecule.from_molblock(basismol_block)
# mol_join = base.substitute(60, grop, 1)

# print(mol_join, file=open(f"join-7.mol", "w"))




smiles_list = open('re-wh.smi','r').read().split()
for idx,smi in enumerate(smiles_list):
    mb = smi2mol(smi)
    grop = Molecule.from_molblock(mb)

    basismol_block = open('base-rdg.mol','r').read()
    base = Molecule.from_molblock(basismol_block)
    mol_join = base.substitute(60, grop, 1)

    print(mol_join, file=open(f"re-wh/rdg/join-{idx}.mol", "w"))
