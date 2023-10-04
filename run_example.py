from SmilesInputs import ReadInput
from esp_generator_wrapper import generate_esps
from conformer_gen import Conformers
from openff.recharge.utilities.molecule import smiles_to_molecule



def main():
    
    #Read the .smi input and add to list
    smiles = ReadInput.read_smiles('test_files.smi')
   
    #Loop through molecules
    for mol in smiles:
        molecule = smiles_to_molecule(mol)
        #Generate the conformers
        conformer_list = Conformers.generate(molecule, generation_type='rdkit')
        ESP_gen = generate_esps(molecule = molecule, conformers= conformer_list)
        ESP_gen.run_esps()
        

if __name__ == "__main__":
    main()
