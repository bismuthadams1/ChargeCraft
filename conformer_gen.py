from openff.toolkit import Molecule
from rdkit.Chem import AllChem
from openff.toolkit import Molecule
from rdkit import Chem
import numpy
from openff.units import unit
from openff.recharge.conformers import ConformerGenerator, ConformerSettings


class Conformers:
    """
    class for generating conformers

    """

    @classmethod
    def _generate(
        cls,
        molecule: "Molecule",
        max_conformers: int = 10,
        generation_type: str = 'openeye'
    ) -> list[numpy.array]:
        """
        Takes the molecule in and generates conformers. 

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
            return cls._rdkit_gen(molecule, max_conformers)
        elif generation_type == 'rdkit':
            return cls._openeye_gen(molecule, max_conformers)
        else:
            return 'invalid conformer gen'
    

    @classmethod
    def _rdkit_gen(
            cls,
            molecule: "Molecule", 
            max_conformers: int
    ) -> list[numpy.array]:

        rdmol = molecule.to_rdkit()
        AllChem.EmbedMultipleConfs(rdmol, numConfs=max_conformers, randomSeed=42)

        conformers = []

        for confs in range(max_conformers):
            conformer = numpy.zeros((rdmol.GetConformer(confs).GetNumAtoms(), 3))
            for atom_index, coordinates in enumerate(rdmol.GetConformer(confs).GetPositions()):
                conformer[atom_index, :] = coordinates
            conformers.append(conformer * unit.angstrom)

        return conformers

    @classmethod
    def _openeye_gen(
        cls,
        molecule: "Molecule",
        max_conformers: int
    ) -> list[numpy.array]:
        
        conformers = ConformerGenerator.generate(
        molecule, ConformerSettings(max_conformers=max_conformers))


        return conformers
    
