import numpy as np
from typing import List
from scipy.spatial.transform import Rotation
from copy import deepcopy


class Atom:
    def __init__(self, symbol: str, index: int, cartesian_coordinate: np.ndarray) -> None:
        self.symbol = symbol
        self.index = index
        self.coordinate = cartesian_coordinate

    def translate(self, vector: np.ndarray):
        self.coordinate += vector

        return self

    def rotate(self, center: np.ndarray, axis: np.ndarray, angle: float):
        self.translate(center)
        r = Rotation.from_rotvec(angle * axis / np.linalg.norm(axis))
        self.coordinate = r.apply(self.coordinate)
        self.translate(-center)

        return self


class Bond:
    def __init__(self, index: int, atoms: List[Atom], order: int) -> None:
        self.index = index
        self.atoms = atoms
        self.order = order


class Molecule:
    def __init__(self, atoms: List[Atom], bonds: List[Bond]) -> None:
        self.atoms = atoms
        self.bonds = bonds

    def remove_atoms(self, atomsindices: List[int]):
        self.atoms = [self.atoms[i] for i in range(len(self.atoms)) if i not in atomsindices]
        self.bonds = [
            self.bonds[i]
            for i in range(len(self.bonds))
            if (self.bonds[i].atoms[0].index not in atomsindices) and (self.bonds[i].atoms[1].index not in atomsindices)
        ]

    def remove_bonds(self, bondsindices: List[int]):
        self.bonds = [self.bonds[i] for i in range(len(self.bonds)) if i not in bondsindices]

    def add_atoms(self, atoms: List[Atom]):
        self.add_atoms += atoms

    def add_bonds(self, bonds: List[Bond]):
        self.bonds += bonds

    def translate(self, vector: np.ndarray):
        self.atoms = [atom.translate(vector) for atom in self.atoms]
        return self

    def rotate(self, center: np.ndarray, axis: np.ndarray, angle: float):
        self.atoms = [atom.rotate(center, axis, angle) for atom in self.atoms]

    @classmethod
    def from_molfile(cls, molfile):
        with open(molfile, "r") as f:
            molblock = f.read()
        return cls.from_molblock(molblock)

    @classmethod
    def from_molblock(cls, molblock: str):
        lines = molblock.splitlines()

        num_atoms = int(lines[3].split()[0])
        num_bonds = int(lines[3].split()[1])
        lines_atoms = lines[4 : 4 + num_atoms]
        lines_bonds = lines[4 + num_atoms : 4 + num_atoms + num_bonds]

        atoms = []
        for idx, line in enumerate(lines_atoms):
            l = line.split()
            coordinate = np.array([float(x) for x in l[:3]])
            symbol = l[3]
            atoms.append(Atom(symbol, idx, coordinate))
        bonds = []
        for idx, line in enumerate(lines_bonds):
            l = line.split()
            a1 = int(l[0])
            a2 = int(l[1])
            a3 = int(l[2])
            bonds.append(Bond(idx, [atoms[a1 - 1], atoms[a2 - 1]], a3))

        return cls(atoms, bonds)

    def __add__(self, other):
        # new_atoms = deepcopy(self.atoms) + deepcopy(other.atoms)
        # new_bonds = deepcopy(self.bonds) + deepcopy(other.bonds)

        new_atoms = self.atoms + other.atoms
        new_bonds = self.bonds + other.bonds
        return Molecule(new_atoms, new_bonds)

    def __str__(self) -> str:
        lines = [
            "Generated by molutils",
            " OpenBabel11202318493D",
            "",
            f"{len(self.atoms):>3d}{len(self.bonds):>3d}  0  0  0  0  0  0  0  0999 V2000",
        ]
        for atom in self.atoms:
            lines.append(
                f" {atom.coordinate[0]:>9.4f} {atom.coordinate[1]:>9.4f} {atom.coordinate[2]:>9.4f} {atom.symbol:<2s}  0  0  0  0  0  0  0  0  0  0  0  0"
            )

        for bond in self.bonds:
            lines.append(
                f"{self.atoms.index(bond.atoms[0])+1:>3d}{self.atoms.index(bond.atoms[1])+1:>3d}{bond.order:>3d}  0  0  0  0"
            )

        return "\n".join(lines)


if __name__ == "__main__":
    a = Molecule.from_molfile("1.mol")
    b = deepcopy(a)
    b = b.translate(np.array([0, 0, 10]))
    c = b + a
    bond = Bond(23, [c.atoms[13], c.atoms[13 + 64]], 878)
    c.add_bonds([bond])
    print(c, file=open("2.mol", "w"))