from SmilesInputs import ReadInput
from esp_generator_wrapper import generate_esps
from conformer_gen import Conformers
from openff.recharge.utilities.molecule import smiles_to_molecule


smiles = ReadInput.read_smiles('test_files.smi')

smiles_molecule_dict = {}

for mol in smiles:
    molecule = smiles_to_molecule(mol)
    conformer_list = Conformers.generate(molecule)
    smiles_molecule_dict[molecule] = conformer_list


