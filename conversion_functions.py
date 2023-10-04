from openff.units.elements import SYMBOLS
from openff.toolkit import Molecule
from rdkit.Chem import rdmolfiles
import numpy as np
from openff.units import unit


def conf_to_xyz_string(conformer: unit.units.Quantity, 
                       molecule: "Molecule") -> str:
    """
    Converts confermer to xyz string.

    Parameters
    ----------
    conformer
        The conformer to generate the xyz string from. 
    molecule
        The molecule object containing the conformer. 
        
    Returns
    -------
        xyz string of the conformer.
    """
    #extract the co-ordinates from the conformer array and the atoms from the openff molecule object. 
    atoms = [
                    {
                        "element": SYMBOLS[atom.atomic_number],
                        "x": conformer[index, 0],
                        "y": conformer[index, 1],
                        "z": conformer[index, 2],
                    }
                    for index, atom in enumerate(molecule.atoms)
                ]
    #format the top of the xyz file with the number of atoms and the smiles string
    xyz = f'{molecule.n_atoms}\n{molecule.to_smiles()}\n'
    for row in atoms:
        xyz += f"{row['element']}\t{np.around(row['x'].magnitude,decimals=8)}\t{np.around(row['y'].magnitude,decimals=8)}\t{np.around(row['z'].magnitude, decimals=8)}\n"
    return xyz

def xyz_string_to_conf(xyz_string, 
                       molecule) -> unit.units.Quantity:
    """
    Takes the xyz string in and generates a conformer object compatabile with the openff Molecule object with units Angstrom.

    Parameters
    ----------
    xyz_string
        The xyz string of the conformer.
    molecule
        The molecule object containing the conformer. 
    
    Returns
    -------
        The conformer object.
    """
    #skip the heading of the xyz_string 
    xyz_list = np.array([row.split('\t') for row in xyz_string.split('\n')[2:]][:-1])

    conformer_from_xyz = np.zeros((molecule.n_atoms, 3))

    for atom_index, coordinates in enumerate(xyz_list):
        #append to conformer array, skip the element in the xyz file
        conformer_from_xyz[atom_index, :] = coordinates[1:] 

    conformer_from_xyz = conformer_from_xyz * unit.angstrom

    return conformer_from_xyz

