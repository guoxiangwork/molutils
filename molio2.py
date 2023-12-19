import numpy as np
from atom import Atom,Bond


def read_atomline(atomline:str):
    line = atomline.split()
    coordinate  = np.array([float(x) for x in line[:3]])
    symbol = line[4]

    return  symbol,coordinate


def read_bondline(bondline:str):
    line = bondline.split()


def read_molblcok(molblock:str):
    lines = molblock.splitlines()

    num_atoms = int(lines[3].split()[0])
    num_bonds = int(lines[3].split()[1])
    lines_atoms = lines[4 : 4 + num_atoms]
    lines_bonds = lines[4 + num_atoms : 4 + num_atoms + num_bonds]


    atoms = []
    for idx,line in enumerate(lines_atoms):
        l = line.split()
        coordinate  = np.array([float(x) for x in l[:3]])
        symbol = line[4]
        atoms.append(Atom(idx,symbol,coordinate))

    bonds = []
    for idx,line in enumerate(lines_bonds):
        l = line.split()
        a1 = int(l[0])
        a2 = int(l[1])
        a3 = int(l[2])
        bonds.append(Bond(idx,[atoms[a1],atoms[a2]],a3))

    return atoms, bonds