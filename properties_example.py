import sys
import os
CWD = os.getcwd()
sys.path.append(os.path.dirname(CWD))

from source.inputSetup.SmilesInputs import ReadInput
from source.optimize.esp_generator_wrapper import ESPGenerator, PropGenerator
from source.conformers.conformer_gen import Conformers
from openff.recharge.utilities.molecule import smiles_to_molecule
from openff.recharge.grids import LatticeGridSettings
from openff.recharge.esp import ESPSettings
from source.storage.storage import MoleculePropRecord, MoleculePropStore



def main():
    
    #Read the .smi input and add to list
    smiles = ReadInput.read_smiles('test_files.smi')

    # Define the grid that the electrostatic properties will be trained on and the
    # level of theory to compute the properties at.
    grid_settings = LatticeGridSettings(
        type="fcc", spacing=0.5, inner_vdw_scale=1.4, outer_vdw_scale=2.0
    )

    esp_settings = ESPSettings(basis="6-31G*", method="hf", grid_settings=grid_settings)

    #Loop through molecules
    for mol in smiles:
        molecule = smiles_to_molecule(mol)
        #Generate the conformers
        conformer_list = Conformers.generate(molecule, generation_type='rdkit', max_conformers=10)
        ESP_gen = PropGenerator(molecule = molecule, conformers = conformer_list, esp_settings = esp_settings, grid_settings = grid_settings, prop_data_store = MoleculePropStore(database_path='prop_test_2.db'))
        ESP_gen.memory = 2e+9 #2gb
        print(f'number of cores is {ESP_gen.ncores}')
        print(f'memory is {ESP_gen.memory}')
        ESP_gen.run_props()

if __name__ == "__main__":
    main()
