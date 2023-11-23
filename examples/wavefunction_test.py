import os
import sys
import random
import numpy as np

from tqdm import tqdm
from rdkit.Chem import AllChem
from rdkit.Chem.rdForceFieldHelpers import MMFFOptimizeMolecule

from openff.units.elements import SYMBOLS
from openff.toolkit import Molecule
from openff.recharge.conformers import ConformerGenerator, ConformerSettings
from source.storage.ddx_storage import ESPSettings, DDXSettings
from openff.recharge.esp.storage import MoleculeESPRecord, MoleculeESPStore
from openff.recharge.grids import LatticeGridSettings, GridSettingsType, GridGenerator
from openff.recharge.utilities.molecule import smiles_to_molecule
from qcelemental.models.procedures import OptimizationInput, QCInputSpecification
from openff.units import unit
from qcelemental.models.common_models import Model

sys.path.append('/Users/localadmin/Documents/projects/QM_ESP_Psi4')

from source.optimize.openff_psi4_gen import Psi4Generate
from source.conformers.conformer_gen import Conformers
from source.utilities.conversion_functions import conf_to_xyz_string


import psi4

sys.settrace 

test_mol =  smiles_to_molecule('[H]O[H]')

#generate water test molecule as openff.toolkit.Molecule
test_mol =  smiles_to_molecule('[H]O[H]')
conformer_list = Conformers.generate(test_mol, generation_type='rdkit')
conformer_list[0]
qc_mol =  test_mol.to_qcschema(conformer=0)

#Generate grid.dat file for grid_esp and grid_field
grid_settings = LatticeGridSettings(
        type="fcc", spacing=0.5, inner_vdw_scale=1.4, outer_vdw_scale=2.0
    )


esp_settings = ESPSettings(basis="6-31G*", method="hf", grid_settings=grid_settings)

xyz = conformer_string = conf_to_xyz_string(conformer_list[0], test_mol)

psi4.set_output_file('output.dat')

molecule = psi4.geometry(
"""
3
[H]O[H]
H	-0.81024874	-0.18571238	-0.0
O	-0.00248133	0.36959931	-0.0
H	0.81273007	-0.18388693	0.0
"""
)

molecule.set_molecular_charge(0)
molecule.set_multiplicity(1)

E, wfn = psi4.energy('hf/6-31g*', return_wfn=True, molecule = molecule)
psi4.prop('hf/6-31G*', properties=["MULLIKEN_CHARGES", "LOWDIN_CHARGES", "DIPOLE", "QUADRUPOLE", "MBIS_CHARGES"])  #"GRID_ESP", "GRID_FIELD"
#psi4.core.print_variables()
psi4.core.variables(['MBIS Charges'])
#print(wfn.atomic_point_charges())