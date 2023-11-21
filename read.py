from typing import List
import numpy as np

def read_molfile(mol_file):
    with open(mol_file) as f:
        lines = f.readlines()

    line1 = lines[0].strip()
    line2 = lines[1].strip()

    num_atoms = int(lines[3].split()[0])
    num_bonds = int(lines[3].split()[1])

    lines_atoms = lines[4 : 4 + num_atoms]
    lines_bonds = lines[5 + num_atoms : 4 + num_atoms + num_bonds]
    e, c = read_atoms_lines(lines_atoms)
    p, t = read_bonds_lines(lines_bonds)
    return e, c, p, t


def read_molblock(mol_block:str):
    lines = mol_block.splitlines()

    line1 = lines[0].strip()
    line2 = lines[1].strip()

    num_atoms = int(lines[3].split()[0])
    num_bonds = int(lines[3].split()[1])

    lines_atoms = lines[4 : 4 + num_atoms]
    lines_bonds = lines[5 + num_atoms : 4 + num_atoms + num_bonds]
    e, c = read_atoms_lines(lines_atoms)
    p, t = read_bonds_lines(lines_bonds)
    return e, c, p, t


def read_atoms_lines(lines_atoms: List[str]):
    elements = []
    coordinates = []
    for line in lines_atoms:
        l = line.split()
        elements.append(l[3])
        coordinates.append([float(x) for x in l[0:3]])

    elements = np.array(elements, dtype="<U2")
    coordinates = np.array(coordinates)
    return elements, coordinates


def read_bonds_lines(lines_bonds: List[str]):
    bonds_pair = []
    bonds_type = []
    for line in lines_bonds:
        l = line.split()
        bonds_type.append(int(l[2]))
        bonds_pair.append([int(x) for x in l[0:2]])

    bonds_type = np.array(bonds_type)
    bonds_pair = np.array(bonds_pair)
    return bonds_pair, bonds_type