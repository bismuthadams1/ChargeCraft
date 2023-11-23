import sys
import os
CWD = os.getcwd()
sys.path.append(os.path.dirname(CWD))
sys.path.append('/Users/localadmin/Documents/projects/QM_ESP_Psi4/')
print(sys.path)

from openff.recharge.esp.storage import MoleculeESPRecord, MoleculeESPStore
from chargecraft.storage.storage import MoleculePropRecord, MoleculePropStore

prop_store = MoleculePropStore("/Users/localadmin/Documents/projects/QM_ESP_Psi4/examples/prop-store copy.sqlite")
smiles_list = prop_store.list()
print(smiles_list)
prop_store.retrieve(smiles_list[0])[0].grid_coordinates
