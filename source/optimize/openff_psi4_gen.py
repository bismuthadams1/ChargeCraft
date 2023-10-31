"""Compute ESP and electric field data using Psi4"""
import os
import subprocess
from typing import TYPE_CHECKING, Optional, Tuple

import jinja2
import numpy
from openff.units import unit
from openff.units.elements import SYMBOLS
from openff.utilities import get_data_file_path, temporary_cd
from psi4.core import Wavefunction

from openff.recharge.esp import ESPGenerator, ESPSettings
from openff.recharge.esp.exceptions import Psi4Error
from source.utilities.conversion_functions import conf_to_xyz_string

if TYPE_CHECKING:
    from openff.toolkit import Molecule

from qcelemental.models import Molecule as QCMolecule
from psi4.core import Molecule as Psi4Molecule

import psi4

CWD = os.getcwd()

class CustomPsi4ESPGenerator:
    """An class which will compute the electrostatic potential of
    a molecule using Psi4.
    """

    @classmethod
    def _generate_input(
        cls,
        molecule: "Molecule",
        conformer: unit.Quantity,
        settings: ESPSettings,
        minimize: bool,
        compute_esp: bool,
        compute_field: bool,
        dynamic_level: int = 1
    ) -> str:
        """Generate the input files for Psi4.

        Parameters
        ----------
        molecule
            The molecule to generate the ESP for.
        conformer
            The conformer of the molecule to generate the ESP for.
        settings
            The settings to use when generating the ESP.
        minimize
            Whether to energy minimize the conformer prior to computing the ESP using
            the same level of theory that the ESP will be computed at.
        compute_esp
            Whether to compute the ESP at each grid point.
        compute_field
            Whether to compute the field at each grid point.

        Returns
        -------
            The contents of the input file.
        """
        # Compute the total formal charge on the molecule.
        # Trust that it's in units of elementary charge.
        formal_charge = sum(atom.formal_charge for atom in molecule.atoms).m

        # Compute the spin multiplicity
        total_atomic_number = sum(atom.atomic_number for atom in molecule.atoms)
        spin_multiplicity = 1 if (formal_charge + total_atomic_number) % 2 == 0 else 2

        # Store the atoms and coordinates in a jinja friendly dict.
        conformer = conformer.to(unit.angstrom).m

        atoms = [
            {
                "element": SYMBOLS[atom.atomic_number],
                "x": conformer[index, 0],
                "y": conformer[index, 1],
                "z": conformer[index, 2],
            }
            for index, atom in enumerate(molecule.atoms)
        ]
        

        #

        # Format the jinja template
        template_path = get_data_file_path(
            os.path.join("psi4", "input.dat"), "openff.recharge"
        )

        with open(template_path) as file:
            template = jinja2.Template(file.read())

        enable_pcm = settings.pcm_settings is not None

        properties = []

        if compute_esp:
            properties.append("GRID_ESP")
        if compute_field:
            properties.append("GRID_FIELD")
        
        properties.extend(["MULLIKEN_CHARGES", "LOWDIN_CHARGES", "DIPOLE", "QUADRUPOLE", "MBIS_CHARGES"])

        template_inputs = {
            "charge": formal_charge,
            "spin": spin_multiplicity,
            "atoms": atoms,
            "basis": settings.basis,
            "method": settings.method,
            "enable_pcm": enable_pcm,
            "dft_settings": settings.psi4_dft_grid_settings.value,
            "minimize": minimize,
            "properties": str(properties),
            "dynamic_level": dynamic_level,
            "return_wfn":"True"
        }

        if enable_pcm:
            template_inputs.update(
                {
                    "pcm_solver": settings.pcm_settings.solver,
                    "pcm_solvent": settings.pcm_settings.solvent,
                    "pcm_radii_set": settings.pcm_settings.radii_model,
                    "pcm_scaling": settings.pcm_settings.radii_scaling,
                    "pcm_area": settings.pcm_settings.cavity_area,
                }
            )

        rendered_template = template.render(template_inputs)
        # Remove the white space after the for loop
        rendered_template = rendered_template.replace("  \n}", "}")

        return rendered_template

    @classmethod
    def generate(
        cls,
        molecule: "Molecule", 
        conformer: unit.Quantity,
        grid: unit.Quantity,
        settings: ESPSettings,
        dynamic_level: int,
        directory: str = CWD,
        minimize: bool = True,
        compute_esp: bool = True,
        compute_field: bool = True,
    ) -> Tuple[unit.Quantity, Optional[unit.Quantity], Optional[unit.Quantity]]:
        # Perform the calculation in a temporary directory


        with temporary_cd(directory):
            # Store the input file.
            input_contents = cls._generate_input(
                molecule, conformer, settings, minimize, compute_esp, compute_field
            )

            # Store the grid to file.
            grid = grid.to(unit.angstrom).m
            numpy.savetxt("grid.dat", grid, delimiter=" ", fmt="%16.10f")

            # Attempt to run the calculation
            psi4_process = subprocess.Popen(
                ["psi4", "input.dat", "output.dat"],
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
            )

            std_output, std_error = psi4_process.communicate()
            exit_code = psi4_process.returncode

            if exit_code != 0:
                raise Psi4Error(std_output.decode(), std_error.decode())

            esp, electric_field = None, None

            if compute_esp:
                esp = (
                    numpy.loadtxt("grid_esp.dat").reshape(-1, 1) * unit.hartree / unit.e
                )
            if compute_field:
                electric_field = (
                    numpy.loadtxt("grid_field.dat")
                    * unit.hartree
                    / (unit.e * unit.bohr)
                )

            with open("final-geometry.xyz") as file:
                output_lines = file.read().splitlines(False)

            final_coordinates = (
                numpy.array(
                    [
                        [
                            float(coordinate)
                            for coordinate in coordinate_line.split()[1:]
                        ]
                        for coordinate_line in output_lines[2:]
                        if len(coordinate_line) > 0
                    ]
                )
                * unit.angstrom
            )


        return final_coordinates, grid, esp, electric_field
  

class Psi4Generate:
    """
    A class that will compute the one electron properties of the wavefunction including ESPs, dipoles, quadroples, grids, milliken, lowdin, and mbis charges.
    """


    @classmethod
    def get_properties(cls,
                molecule: "Molecule",
                conformer: "QCMolecule",
                grid: unit.Quantity,
                settings: ESPSettings,
                #compute_esp: bool,
                #compute_field: bool,
                dynamic_level: int = 1,
                directory: str = CWD,
                ) -> dict:
        
        with temporary_cd(directory):
            grid = grid.to(unit.angstrom).m
            numpy.savetxt("grid.dat", grid, delimiter=" ", fmt="%16.10f")

            formal_charge = sum(atom.formal_charge for atom in molecule.atoms).m

            # Compute the spin multiplicity
            total_atomic_number = sum(atom.atomic_number for atom in molecule.atoms)
            spin_multiplicity = 1 if (formal_charge + total_atomic_number) % 2 == 0 else 2

            #molecule_from_string = conformer.to_string("psi4")

            molecule_psi4 = psi4.geometry(conformer.to_string("psi4"))

            #Ultrafine grid
            psi4.set_options({"DFT_SPHERICAL_POINTS":"590",
                              "DFT_RADIAL_POINTS":"99"})

            #Load in geometry from string, nocom and noreorient to stop automatic  reoerientation
            #molecule_psi4 = psi4.geometry(f"""
            #                        noreorient
            #                       {molecule_from_string}""")
            
            molecule_psi4.set_molecular_charge(formal_charge)
            molecule_psi4.set_multiplicity(spin_multiplicity)

            E, wfn =  psi4.prop('hf/6-31G*', properties=["GRID_ESP",
                                                         "GRID_FIELD",
                                                         "MULLIKEN_CHARGES", 
                                                         "LOWDIN_CHARGES", 
                                                         "DIPOLE", 
                                                         "QUADRUPOLE", 
                                                         "MBIS_CHARGES"], 
                                                         molecule = molecule_psi4,
                                                         return_wfn = True)


            
            esp = (
                numpy.loadtxt("grid_esp.dat").reshape(-1, 1) * unit.hartree / unit.e
            )
        
            electric_field = (
                numpy.loadtxt("grid_field.dat")
                * unit.hartree
                / (unit.e * unit.bohr)
            )

            variable_names = ["MULLIKEN_CHARGES", "LOWDIN_CHARGES", "HF DIPOLE", "HF QUADRUPOLE", "MBIS CHARGES"]
            variables_dictionary = {name: wfn.variable(name) for name in variable_names}

            #xyz_coords = Psi4Molecule.xyz(molecule_psi4)
            xyz_coords = molecule_psi4.print_out_in_angstrom()
            #print(molecule_psi4.print_out_in_angstrom())

            return xyz_coords, grid, esp, electric_field, variables_dictionary
