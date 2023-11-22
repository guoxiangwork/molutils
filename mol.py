from pathlib import Path
from typing import List
import numpy as np
from scipy.spatial.transform import Rotation
from .molio import *



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



def transition(m: MolFile, vec: np.ndarray):
    m.coordinates = m.coordinates + vec

    return m

def transition_molblock(mol_block:str):
    e,c,p,t=read_mol_block(mol_block)


def rotation(m: MolFile, pole, point: np.ndarray, angle):
    coordinates = m.coordinates - point
    pole_norm = pole / np.linalg.norm(pole)

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

