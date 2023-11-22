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


def read_mol_block(mol_block:str):
    lines = mol_block.splitlines()

    num_atoms = int(lines[3].split()[0])
    num_bonds = int(lines[3].split()[1])

    lines_atoms = lines[4 : 4 + num_atoms]
    lines_bonds = lines[4 + num_atoms : 4 + num_atoms + num_bonds]
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





def write_mol_file(e,c,p,t):

    lines = []
    lines.append("made by guoxiang python script")
    lines.append("3D molecuel")
    lines.append("")
    lines.append(
        f" {len(e)} {len(t)}  0  0  0  0  0  0  0  0999 V2000"
    )
    atoms_lines = [
        f"{item[1][0]:>9.4f} {item[1][1]:>9.4f} {item[1][2]:>9.4f} {item[0]:<2s}  0  0  0  0  0  0  0  0  0  0  0  0"
        for item in zip(e, c)
    ]
    bonds_lines = [
        f"{item[0][0]:>3d}{item[0][1]:>3d}{item[1]:>3d}  0  0  0  0"
        for item in zip(p, t)
    ]
    lines += atoms_lines
    lines += bonds_lines
    lines.append("M  END")
    return "\n".join(lines)



if __name__=="__main__":
    with open('te-f.mol') as f:
        mol_block = f.read()
    e,c,p,t = read_mol_block(mol_block)
    print(p)
    m=write_mol_file(e,c,p,t)
    print(m)