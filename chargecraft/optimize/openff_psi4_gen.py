"""Compute ESP and electric field data using Psi4"""
import os
import subprocess
from typing import TYPE_CHECKING, Optional, Tuple

import jinja2
import numpy
from openff.units import unit
from openff.units.elements import SYMBOLS
from openff.utilities import get_data_file_path, temporary_cd

from chargecraft.storage.data_classes import ESPSettings, DDXSettings


from openff.recharge.esp.exceptions import Psi4Error
from chargecraft.globals import GlobalConfig

from openff.toolkit import Molecule

from qcelemental.models import Molecule as QCMolecule

import psi4
from psi4.core import GeometryUnits
from psi4.core import Options

CWD = os.getcwd()


class Psi4Generate:
    """
    A class that will compute the one electron properties of the wavefunction including ESPs, dipoles, quadroples, grids, milliken, lowdin, and mbis charges.
    """

   
    @classmethod
    def get_properties(cls,
                molecule: "Molecule",
                conformer: "QCMolecule",
                grid: unit.Quantity,
                settings: ESPSettings ,
                dynamic_level: int = 1,
                directory: str = CWD,
                extra_options: dict[any] | None = None
                ) -> Tuple[unit.Quantity, unit.Quantity, unit.Quantity, unit.Quantity, dict, int]:
        
        with temporary_cd(directory):
            grid = grid.to(unit.angstrom).m
            numpy.savetxt("grid.dat", grid, delimiter=" ", fmt="%16.10f")

            formal_charge = sum(atom.formal_charge for atom in molecule.atoms).m

            # Compute the spin multiplicity
            total_atomic_number = sum(atom.atomic_number for atom in molecule.atoms)
            spin_multiplicity = 1 if (formal_charge + total_atomic_number) % 2 == 0 else 2

            conformer_Ang = Molecule.from_qcschema(conformer)

            conformer_Ang = conformer_Ang.conformers[0].to(unit.angstrom).m

            #TODO check reorientation here
            conformer_Ang_string = "nocom\nnoreorient\n"
            for index, atom in enumerate(molecule.atoms):
                  conformer_Ang_string += f"{SYMBOLS[atom.atomic_number]}\t{conformer_Ang[index, 0]}\t{conformer_Ang[index, 1]}\t{conformer_Ang[index, 2]}\n"
       
            molecule_psi4 = psi4.geometry(conformer_Ang_string.strip())
            print(molecule_psi4)
            molecule_psi4.set_units(GeometryUnits.Angstrom)
            
            #TODO set number of threads somewhere
            #Reset options between runs
            psi4.core.clean_options()
            #Ultrafine grid
            psi4.set_options({"DFT_SPHERICAL_POINTS":"590",
                              "DFT_RADIAL_POINTS":"99"})
            #Set additional options
            if extra_options:
                print(f'setting extra options: {extra_options}')
                psi4.set_options(extra_options)
            #Number of threads should be the number of cores * num of threads per core
            psi4.set_num_threads(GlobalConfig().total_threads())
            psi4.set_memory(GlobalConfig().memory())
            print(f'number of threads is {GlobalConfig().total_threads()}')
            enable_solvent = settings.pcm_settings is not None or settings.ddx_settings is not None
            print(f'settings pcm {settings.pcm_settings}')
            print(f'settings ddx {settings.ddx_settings}')
            print(f'enable solvent: {enable_solvent}')
            if enable_solvent:
                print('setting solvent...')
                if settings.pcm_settings is not None and settings.pcm_settings.solver is not None:
                            psi4.set_options({"pcm":True})
                            psi4.set_options({ "pcm__input":  f"""
                                            Units = Angstrom
                                            Medium {{
                                                SolverType = {settings.pcm_settings.solver}
                                                Solvent = {settings.pcm_settings.solvent}
                                            }}

                                            Cavity {{
                                                RadiiSet = {settings.pcm_settings.radii_model} # Bondi | UFF | Allinger
                                                Type = GePol
                                                Scaling = {settings.pcm_settings.radii_scaling} # radii for spheres scaled by 1.2
                                                Area = {settings.pcm_settings.cavity_area}
                                                Mode = "Implicit"
                                            }}
                                            """} )
                else:
                    #check if dialetric constant is specified or not
                    if settings.ddx_settings.epsilon is not None:
                        print('ddx numeric option')
                        psi4.set_options({"ddx": "true", #supply a solvent here to see if epsilon then gets picked up
                        "DDX_SOLVENT_EPSILON": settings.ddx_settings.epsilon,
                        "DDX_RADII_SET": settings.ddx_settings.radii_set,
                        "DDX_MODEL": settings.ddx_settings.ddx_model})
                    else:
                        print('ddx solvent option')
                        psi4.set_options({"ddx": "true",
                        "DDX_SOLVENT": settings.ddx_settings.solvent,
                        "DDX_RADII_SET": settings.ddx_settings.radii_set,
                        "DDX_MODEL": settings.ddx_settings.ddx_model})
                          

            molecule_psi4.set_molecular_charge(formal_charge)
            molecule_psi4.set_multiplicity(spin_multiplicity)
            #Currently QM settings hard coded in, can get from ESPSettings object.
            try:
                E, wfn =  psi4.energy(f'{settings.method}/{settings.basis}', molecule = molecule_psi4, return_wfn = True)
                psi4.oeprop(wfn,"GRID_ESP","GRID_FIELD","MULLIKEN_CHARGES", "LOWDIN_CHARGES", "DIPOLE","QUADRUPOLE", "MBIS_CHARGES")

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
                print(wfn.variables())

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
                #different indexes for dipole if dft vs hf method
                variables_dictionary["DIPOLE"] = wfn.variable("SCF DIPOLE") * unit.e * unit.bohr_radius
                variables_dictionary["QUADRUPOLE"] = wfn.variable("QUADRUPOLE") * unit.e * unit.bohr_radius**2
                variables_dictionary["ALPHA_DENSITY"] = wfn.Da().to_array()
                variables_dictionary["BETA_DENSITY"] = wfn.Db().to_array()
            
                #qcelemental.geometry is outputted in bohr, convert to  angstrom
                final_coordinates = (conformer.geometry * unit.bohr).to(unit.angstrom)
                #Cleanup scratch files
                psi4.core.clean()


                return final_coordinates, grid, esp, electric_field, variables_dictionary, E
            except Exception as e:
                    print(e)
                    psi4.core.clean()
                    return Psi4Error
  