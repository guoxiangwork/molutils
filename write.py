from .mol import MolFile



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