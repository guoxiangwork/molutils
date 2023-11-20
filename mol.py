from pathlib import Path
from typing import List
import numpy as np
from scipy.spatial.transform import Rotation


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


class MolFile:
    def __init__(self, mol_file_path) -> None:
        if type(mol_file_path) == str:
            self.mol_file_path = Path(mol_file_path)
        elif type(mol_file_path) == Path:
            self.mol_file_path = mol_file_path

        else:
            raise

        self.elements, self.coordinates, self.bond_pair, self.bond_type = read_molfile(
            self.mol_file_path
        )


def write_mol_file(m: MolFile):
    lines = []
    lines.append("made by guoxiang python script")
    lines.append("3D molecuel")
    lines.append("")
    lines.append(
        f" {len(m.elements)} {len(m.bond_type)}  0  0  0  0  0  0  0  0999 V2000"
    )
    atoms_lines = [
        f"{item[1][0]:>9.4f} {item[1][1]:>9.4f} {item[1][2]:>9.4f} {item[0]:<2s}  0  0  0  0  0  0  0  0  0  0  0  0"
        for item in zip(m.elements, m.coordinates)
    ]
    bonds_lines = [
        f"{item[0][0]:>4d}  {item[0][1]:>4d}  {item[1]:>1d}  0  0  0  0"
        for item in zip(m.bond_pair, m.bond_type)
    ]
    lines += atoms_lines
    lines += bonds_lines
    lines.append("M  END")
    return "\n".join(lines)


def transition(m: MolFile, vec: np.ndarray):
    m.coordinates = m.coordinates + vec

    return m


def rotation(m: MolFile, pole, point: np.ndarray, angle):
    old_coordiantes = m.coordinates
    coordinates = m.coordinates - point
    pole_norm = pole / np.sqrt(np.linalg.norm(pole, 2))

    theta = angle / 180 * np.pi

    r = Rotation.from_rotvec(theta * pole_norm)
    c = r.apply(coordinates)

    return c + point




def combine_molfile(m1, m2, newbonds: List[list]):
    m1 = MolFile(m1)
    m2 = MolFile(m2)
    elements = np.concatenate([m1.elements, m2.elements])
    coordinates = np.concatenate([m1.coordinates, m2.coordinates], axis=0)
    bond_pair = np.concatenate([m1.bond_pair, m2.bond_pair + len(m1.elements)], axis=0)
    bond_type = np.concatenate([m1.bond_type, m2.bond_type])

    l1 = []
    l2 = []
    for item in newbonds:
        l1.append([item[0], item[1] + len(m1.elements)])
        l2.append(item[2])

    l1 = np.array(l1)
    l2 = np.array(l2)

    bond_pair = np.concatenate([bond_pair, l1], axis=0)
    bond_type = np.concatenate([bond_type, l2], axis=0)

    m1.bond_type = bond_type
    m1.elements = elements
    m1.coordinates = coordinates
    m1.bond_pair = bond_pair

    return write_mol_file(m1)


if __name__ == "__main__":

    m3 = combine_molfile('1.mol', '1.mol', [[1, 1, 1]])
    print(m3)

