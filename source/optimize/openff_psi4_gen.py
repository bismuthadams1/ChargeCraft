"""Compute ESP and electric field data using Psi4"""
import os
import subprocess
from typing import TYPE_CHECKING, Optional, Tuple

import jinja2
import numpy
from openff.units import unit
from openff.units.elements import SYMBOLS
from openff.utilities import get_data_file_path, temporary_cd

from openff.recharge.esp import ESPGenerator, ESPSettings
from openff.recharge.esp.exceptions import Psi4Error
from source.utilities.conversion_functions import conf_to_xyz_string

from openff.toolkit import Molecule

from qcelemental.models import Molecule as QCMolecule

import psi4
from psi4.core import GeometryUnits

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
                dynamic_level: int = 1,
                directory: str = CWD,
                ) -> Tuple[unit.Quantity, unit.Quantity, unit.Quantity, unit.Quantity, dict()]:
        
        with temporary_cd(directory):
            grid = grid.to(unit.angstrom).m
            numpy.savetxt("grid.dat", grid, delimiter=" ", fmt="%16.10f")

            formal_charge = sum(atom.formal_charge for atom in molecule.atoms).m

            # Compute the spin multiplicity
            total_atomic_number = sum(atom.atomic_number for atom in molecule.atoms)
            spin_multiplicity = 1 if (formal_charge + total_atomic_number) % 2 == 0 else 2

            conformer_Ang = Molecule.from_qcschema(conformer)

            conformer_Ang = conformer_Ang.conformers[0].to(unit.angstrom).m

    
            conformer_Ang_string = ""
            for index, atom in enumerate(molecule.atoms):
                  conformer_Ang_string += f"{SYMBOLS[atom.atomic_number]}\t{conformer_Ang[index, 0]}\t{conformer_Ang[index, 1]}\t{conformer_Ang[index, 2]}\n"
       
            molecule_psi4 = psi4.geometry(conformer_Ang_string.strip())
            
            molecule_psi4.set_units(GeometryUnits.Angstrom)

            #Ultrafine grid
            psi4.set_options({"DFT_SPHERICAL_POINTS":"590",
                              "DFT_RADIAL_POINTS":"99"})
            
            molecule_psi4.set_molecular_charge(formal_charge)
            molecule_psi4.set_multiplicity(spin_multiplicity)
            #Currently QM settings hard coded in, can get from ESPSettings object.
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

            #variable_names = ["MULLIKEN_CHARGES", "LOWDIN_CHARGES", "HF DIPOLE", "HF QUADRUPOLE", "MBIS CHARGES"]
            #variables_dictionary = {name: wfn.variable(name) for name in variable_names}
            
            variables_dictionary = dict()
            #psi4 computes charges in a.u., elementary charge
            variables_dictionary["MULLIKEN_CHARGES"] = wfn.variable("MULLIKEN_CHARGES") * unit.e
            variables_dictionary["LOWDIN_CHARGES"] = wfn.variable("LOWDIN_CHARGES") * unit.e
            variables_dictionary["MBIS CHARGES"] = wfn.variable("MBIS CHARGES") * unit.e
            #psi4 grab the MBIS multipoles
            variables_dictionary["MBIS DIPOLE"] = wfn.variable("MBIS DIPOLES") * unit.e * unit.bohr_radius
            variables_dictionary["MBIS QUADRUPOLE"] = wfn.variable("MBIS QUADRUPOLES") * unit.e * unit.bohr_radius**2
            variables_dictionary["MBIS OCTOPOLE"] = wfn.variable("MBIS OCTUPOLES") * unit.e * unit.bohr_radius**3
            #psi4 computes n multipoles in a.u, in elementary charge * bohr radius**n
            variables_dictionary["HF DIPOLE"] = wfn.variable("HF DIPOLE") * unit.e * unit.bohr_radius
            variables_dictionary["HF QUADRUPOLE"] = wfn.variable("HF QUADRUPOLE") * unit.e * unit.bohr_radius**2

            #qcelemental.geometry is outputted in bohr, convert to  angstrom
            final_coordinates = (conformer.geometry * unit.bohr).to(unit.angstrom)

            return final_coordinates, grid, esp, electric_field, variables_dictionary
