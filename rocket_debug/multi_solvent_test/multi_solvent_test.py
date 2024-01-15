import chargecraft
from chargecraft.inputSetup.SmilesInputs import ReadInput
from chargecraft.optimize.esp_generator_wrapper import ESPGenerator, PropGenerator
from chargecraft.conformers.conformer_gen import Conformers
from openff.recharge.utilities.molecule import smiles_to_molecule
from openff.recharge.grids import LatticeGridSettings
from chargecraft.storage.data_classes import ESPSettings, PCMSettings, DDXSettings
from chargecraft.storage.storage import MoleculePropRecord, MoleculePropStore




def main():
    
    #Read the .smi input and add to list
    smiles = ReadInput.read_smiles('test_files.smi')

    # Define the grid that the electrostatic properties will be trained on and the
    # level of theory to compute the properties at.
    grid_settings = LatticeGridSettings(
        type="fcc", spacing=0.5, inner_vdw_scale=1.4, outer_vdw_scale=2.0
    )
    esp_settings_PBE_pcm = ESPSettings(basis="def2-TZVP", method="PBE0", grid_settings=grid_settings, pcm_settings=PCMSettings(solver="IEFPCM",solvent='chloroform'))  #-D3BJ
    # esp_settings_HF = ESPSettings(basis="6-31g*", method="HF", grid_settings=grid_settings) 
    esp_settings_PBE = ESPSettings(basis="def2-TZVP", method="PBE0", grid_settings=grid_settings) #-D3BJ

    esp_settings_PBE_ddx4 = ESPSettings(basis="def2-TZVP", method="PBE0-D3BJ", grid_settings=grid_settings,  ddx_settings = DDXSettings(solvent=4) ) 
    esp_settings_PBE_ddx80 = ESPSettings(basis="def2-TZVP", method="PBE0-D3BJ", grid_settings=grid_settings,  ddx_settings = DDXSettings(solvent=80) ) 
    esp_settings_PW6B95 = ESPSettings(basis="aug-cc-pv_dpd_z", method="PW6B95", grid_settings=grid_settings) 
    esp_settings_PW6B95_pcm = ESPSettings(basis="aug-cc-pv_dpd_z", method="PW6B95", grid_settings=grid_settings, pcm_settings=PCMSettings(solver="IEFPCM",solvent='chloroform'))
    esp_settings_PW6B95_ddx4 = ESPSettings(basis="aug-cc-pv_dpd_z", method="PW6B95", grid_settings=grid_settings, ddx_settings = DDXSettings(solvent=4) )
    esp_settings_PW6B95_ddx80 = ESPSettings(basis="aug-cc-pv_dpd_z", method="PW6B95", grid_settings=grid_settings, ddx_settings = DDXSettings(solvent=80) )

    #Loop through molecules

    for mol in smiles:
        molecule = smiles_to_molecule(mol)
        #Generate the conformers
        conformer_list = Conformers.generate(molecule, generation_type='rdkit', max_conformers=3)
        # ESP_gen_HF = PropGenerator(molecule = molecule, conformers = conformer_list, esp_settings = esp_settings_HF, grid_settings = grid_settings, prop_data_store = MoleculePropStore(database_path='multi_properties_store.db'))
        # ESP_gen_HF.memory = 2e+9 #2gb
        # ESP_gen_HF.run_props()
        ESP_gen_PBE = PropGenerator(molecule = molecule, conformers = conformer_list, esp_settings = esp_settings_PBE, grid_settings = grid_settings, prop_data_store = MoleculePropStore(database_path='multi_properties_store.db'))
        ESP_gen_PBE.memory = 2e+9 #2gb
        ESP_gen_PBE.run_props()
        ESP_gen_PBE_pcm = PropGenerator(molecule = molecule, conformers = conformer_list, esp_settings = esp_settings_PBE_pcm, grid_settings = grid_settings, prop_data_store = MoleculePropStore(database_path='multi_properties_store.db'))
        ESP_gen_PBE_pcm.memory = 2e+9 #2gb
        ESP_gen_PBE_pcm.run_props()
        ESP_gen_PBE_ddx4 = PropGenerator(molecule = molecule, conformers = conformer_list, esp_settings = esp_settings_PBE_ddx4, grid_settings = grid_settings, prop_data_store = MoleculePropStore(database_path='multi_properties_store.db'))
        ESP_gen_PBE_ddx4.memory = 2e+9 #2gb
        ESP_gen_PBE_ddx4.run_props()
        ESP_gen_PBE_ddx80= PropGenerator(molecule = molecule, conformers = conformer_list, esp_settings = esp_settings_PBE_ddx80, grid_settings = grid_settings, prop_data_store = MoleculePropStore(database_path='multi_properties_store.db'))
        ESP_gen_PBE_ddx80.memory =  2e+9 #2gb
        ESP_gen_PBE_ddx80.run_props()
        ESP_gen_PW6B95= PropGenerator(molecule = molecule, conformers = conformer_list, esp_settings = esp_settings_PW6B95, grid_settings = grid_settings, prop_data_store = MoleculePropStore(database_path='multi_properties_store.db'))
        ESP_gen_PW6B95.memory = 2e+9 #2gb
        ESP_gen_PW6B95.run_props()
        ESP_gen_PW6B95_pcm= PropGenerator(molecule = molecule, conformers = conformer_list, esp_settings = esp_settings_PW6B95_pcm, grid_settings = grid_settings, prop_data_store = MoleculePropStore(database_path='multi_properties_store.db'))
        ESP_gen_PW6B95_pcm.memory = 2e+9 #2gb
        ESP_gen_PW6B95_pcm.run_props()
        ESP_gen_PW6B95_ddx4= PropGenerator(molecule = molecule, conformers = conformer_list, esp_settings = esp_settings_PW6B95_ddx4, grid_settings = grid_settings, prop_data_store = MoleculePropStore(database_path='multi_properties_store.db'))
        ESP_gen_PW6B95_ddx4.memory = 2e+9 #2gb
        ESP_gen_PW6B95_ddx4.run_props()
        ESP_gen_PW6B95_ddx80= PropGenerator(molecule = molecule, conformers = conformer_list, esp_settings = esp_settings_PW6B95_ddx80, grid_settings = grid_settings, prop_data_store = MoleculePropStore(database_path='multi_properties_store.db'))
        ESP_gen_PW6B95_ddx80.memory =  2e+9 #2gb
        ESP_gen_PW6B95_ddx80.run_props()

if __name__ == "__main__":
    main()
