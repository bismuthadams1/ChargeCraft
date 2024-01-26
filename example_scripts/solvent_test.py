import chargecraft
from chargecraft.inputSetup.SmilesInputs import ReadInput
from chargecraft.optimize.esp_generator_wrapper import ESPGenerator, PropGenerator
from chargecraft.conformers.conformer_gen import Conformers
from openff.recharge.utilities.molecule import smiles_to_molecule
from openff.recharge.grids import LatticeGridSettings
from chargecraft.storage.data_classes import ESPSettings, PCMSettings
from chargecraft.storage.storage import MoleculePropRecord, MoleculePropStore
import subprocess

def main():
    
    #Read the .smi input and add to list
    smiles = '[H]C([H])([H])SC([H])([H])[H]'
    molecule = smiles_to_molecule(smiles)
    # Define the grid that the electrostatic properties will be trained on and the
    # level of theory to compute the properties at.
    grid_settings = LatticeGridSettings(
        type="fcc", spacing=0.5, inner_vdw_scale=1.4, outer_vdw_scale=2.0
    )
    conformer_list = Conformers.generate(molecule, generation_type='rdkit', max_conformers=1)

    esp_settings_HF = ESPSettings(basis="6-31g*", method="HF", grid_settings=grid_settings) 
    esp_settings_B3LYP =  ESPSettings(basis="6-31g*", method="B3LYP", grid_settings=grid_settings) 
    esp_settings_B3LYP_solvent =  ESPSettings(basis="6-31g*", method="B3LYP", grid_settings=grid_settings, pcm_settings=PCMSettings(solvent='water')) 
    #Loop through molecules
    #Generate the conformers
    ESP_gen_hf = PropGenerator(molecule = molecule, 
                            conformers = conformer_list[0], 
                            esp_settings = esp_settings_HF, grid_settings = grid_settings, prop_data_store = MoleculePropStore(database_path='solvent_store.db'))
    ESP_gen_hf.memory = 2e+9 #2gb
    print(f'number of cores is {ESP_gen_hf.ncores}')
    print(f'memory is {ESP_gen_hf.memory}')
    confomer = ESP_gen_hf.run_props()
    ESP_gen_B3LYP = PropGenerator(molecule = molecule, 
                            conformers = confomer, 
                            esp_settings = esp_settings_B3LYP, 
                            grid_settings = grid_settings, 
                            prop_data_store = MoleculePropStore(database_path='solvent_store.db'))
    ESP_gen_B3LYP.memory = 2e+9 #2gb
    ESP_gen_B3LYP.run_props()
    ESP_gen_B3LYP_solvent = PropGenerator(molecule = molecule, 
                            conformers = confomer, 
                            esp_settings = esp_settings_B3LYP_solvent, 
                            grid_settings = grid_settings, 
                            prop_data_store = MoleculePropStore(database_path='solvent_store.db'))
    ESP_gen_B3LYP_solvent.memory = 2e+9 #2gb
    ESP_gen_B3LYP_solvent.run_props()

if __name__ == "__main__":
    main()
