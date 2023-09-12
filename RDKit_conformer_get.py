from openff.toolkit import Molecule
from rdkit.Chem import AllChem
from rdkit import Chem
import numpy
from openff.units import unit


def generate(molecule, no_conformers = 10):

    rdmol = molecule.to_rdkit()
    AllChem.EmbedMultipleConfs(rdmol, numConfs=no_conformers, randomSeed=42)

    conformers = []

    for i in range(no_conformers):
        conformer = numpy.zeros((rdmol.GetConformer(i).GetNumAtoms(), 3))

        for atom_index, coordinates in enumerate(rdmol.GetConformer(i).GetPositions()):
            conformer[atom_index, :] = coordinates

        conformers.append(conformer * unit.angstrom)

    return conformers

if __name__ == "__main__":
    generate()