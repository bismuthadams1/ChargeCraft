from openff.toolkit import Molecule
from rdkit.Chem import AllChem
import numpy
from openff.units import unit
from openff.recharge.conformers import ConformerGenerator, ConformerSettings


class Conformers:
    """
    class for generating conformers

    """

    @staticmethod
    def generate(
        molecule: "Molecule",
        max_conformers: int = 10,
        generation_type: str = 'openeye'
    ) -> list[numpy.array]:
        """
        Takes the molecule in and calls conformer generation functions. 

        Parameters
        ----------
        molecule
            The molecule to generate the conformers for.
        max_conformers
            The number of conformers to generate.
        generation_type
            The calculator for the conformers which can either openeye (license required) or rdkit
     
        Returns
        -------
            The list of conformers.

        """
        
        if generation_type == 'openeye':
            return Conformers._openeye_gen(molecule, max_conformers)
        elif generation_type == 'rdkit':
            return Conformers._rdkit_gen(molecule, max_conformers)
        else:
            return print('invalid conformer generator option')
    

    @staticmethod
    def _rdkit_gen(
            molecule: "Molecule", 
            max_conformers: int
    ) -> list[numpy.array]:
        """
        Takes the molecule in and generates conformers using RDKit. 

        Parameters
        ----------
        molecule
            The molecule to generate the conformers for.
        max_conformers
            The number of conformers to generate.
     
        Returns
        -------
            The list of conformers.
        """

        rdmol = molecule.to_rdkit()
        AllChem.EmbedMultipleConfs(rdmol, numConfs=max_conformers, randomSeed=42)

        print(f'rdkit has generated {rdmol.GetNumConformers()} conformers from a requested {max_conformers} conformers')

        conformers = []

        for confs in range(max_conformers):
            conformer = numpy.zeros((rdmol.GetConformer(confs).GetNumAtoms(), 3))
            for atom_index, coordinates in enumerate(rdmol.GetConformer(confs).GetPositions()):
                conformer[atom_index, :] = coordinates
            conformers.append(conformer * unit.angstrom)

        #Embed conformer in molecule object
        for conf in conformers:
            molecule.add_conformer(conf)

        return conformers

    @staticmethod
    def _openeye_gen(
        molecule: "Molecule",
        max_conformers: int
    ) -> list[numpy.array]:
        """
        Takes the molecule in and generates conformers using openeye. 

        Parameters
        ----------
        molecule
            The molecule to generate the conformers for.
        max_conformers
            The number of conformers to generate.
     
        Returns
        -------
            The list of conformers.
        """
        conformers = ConformerGenerator.generate(
        molecule, ConformerSettings(max_conformers=max_conformers))
        
        print(f'openeye has generated {len(conformers)} conformers from a requested {max_conformers} conformers')

        #Embed conformer in molecule object
        for conf in conformers:
            molecule.add_conformer(conf)

        return conformers
