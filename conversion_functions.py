from openff.units.elements import SYMBOLS
from rdkit.Chem import rdmolfiles
import numpy as np


def conf_to_xyz_string(conformer, molecule):
    atoms = [
                    {
                        "element": SYMBOLS[atom.atomic_number],
                        "x": conformer[index, 0],
                        "y": conformer[index, 1],
                        "z": conformer[index, 2],
                    }
                    for index, atom in enumerate(molecule.atoms)
                ]
    xyz = f'{molecule.n_atoms}\n{molecule.to_smiles()}\n'
    for row in atoms:
        xyz += f"{row['element']}\t{np.around(row['x'].magnitude,decimals=6)}\t{np.around(row['y'].magnitude,decimals=6)}\t{np.around(row['z'].magnitude, decimals=6)}\n"
    return xyz

#def xyz_string_to_conf(xyz_string, conformer):
    
