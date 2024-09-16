import pytest
from openff.toolkit.topology import Molecule
from openff.units import unit
from chargecraft.storage.data_classes import ESPSettings
from chargecraft.storage.storage import MoleculePropStore, MoleculePropRecord
import numpy as np
from unittest.mock import MagicMock, patch

@pytest.fixture
def setup_molecule_1():
    h2o_1 =  Molecule("[H]O[H]")
    h2o_1_conformer = np.array([
        [-3.31754922,  1.50267657,  0.        ],
        [-2.34754922,  1.50267657,  0.        ],
        [-3.64087903,  2.31550273, -0.41913179]
    ]) * unit.angstrom
    h2o_1.add_conformer(h2o_1_conformer)
    return h2o_1
    
@pytest.fixture
def setup_molecule_2():
    h2o_2 =  Molecule("[H]O[H]")
    h2o_2_conformer = np.array([
        [-3.31754922,  1,  0.        ],
        [-2.34754922,  1,  0.        ],
        [-3.64087903,  2, -0.41913179]
    ]) * unit.angstrom
    h2o_2.add_conformer(h2o_2_conformer)
    return h2o_2
       
@pytest.fixture
def dummy_electrostatics_1():
    esp = [1,2,3,5,5] * unit.hartree/ unit.e
    electric_field = [1,2,3,4,5] * unit.hartree / (unit.e * unit.bohr)
    esp_settings = ESPSettings()
    variables_dictionary = {
        "MULLIKEN_CHARGES": [],  
        "LOWDIN_CHARGES": [],    
        "MBIS CHARGES": [],      
        "DIPOLE": [],            
        "QUADRUPOLE": [],        
        "MBIS DIPOLE": [],       
        "MBIS QUADRUPOLE": [],   
        "MBIS OCTOPOLE": [],     
        "ALPHA_DENSITY": [],     
        "BETA_DENSITY": []       
    }
    return esp, electric_field, esp_settings, variables_dictionary

@pytest.fixture
def dummy_electrostatics_2():
    esp = [1,2,3,5] * unit.hartree/ unit.e
    electric_field = [1,2,3,4] * unit.hartree / (unit.e * unit.bohr)
    esp_settings = ESPSettings()
    variables_dictionary = {
        "MULLIKEN_CHARGES": [],  
        "LOWDIN_CHARGES": [],    
        "MBIS CHARGES": [],      
        "DIPOLE": [],            
        "QUADRUPOLE": [],        
        "MBIS DIPOLE": [],       
        "MBIS QUADRUPOLE": [],   
        "MBIS OCTOPOLE": [],     
        "ALPHA_DENSITY": [],     
        "BETA_DENSITY": []       
    }
    return esp, electric_field, esp_settings, variables_dictionary

def create_mock_molecule_prop_record_1():
        #construct creates a record without validation 
        mock_record_1 = MoleculePropRecord.construct(
        # Set values for the first record
        tagged_smiles = "[O:2][H:1][H:3]"
        conformer = [
            [-3.31754922, 1.50267657, 0.0],
            [-2.34754922, 1.50267657, 0.0],
            [-3.64087903, 2.31550273, -0.41913179]
        ]
        grid_coordinates = [
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0]
        ]
        variables_dictionary = {
            "MULLIKEN_CHARGES": [-0.3, 0.1, 0.2],  
            "LOWDIN_CHARGES": [-0.2, 0.0, 0.2],    
            "MBIS CHARGES": [-0.1, 0.1, 0.0],      
            "DIPOLE": [1.1, 0.5, 0.5],            
            "QUADRUPOLE": [1.2, 0.4, 0.4],        
            "MBIS DIPOLE": [1.3, 0.3, 0.3],       
            "MBIS QUADRUPOLE": [1.4, 0.2, 0.2],   
            "MBIS OCTOPOLE": [1.5, 0.1, 0.1],     
            "ALPHA_DENSITY": [[0.0, 0.1], [0.1, 0.0]],     
            "BETA_DENSITY": [[0.1, 0.0], [0.0, 0.1]]       
        }
        
        energy = -75.0
        charge_model_charges = None
        )
        return mock_record_1

def create_mock_molecule_prop_record_2():
    with patch("chargecraft.storage.storage.MoleculePropRecord.__init__", return_value=None):
        mock_record_2 = MoleculePropRecord()
        mock_record_2.tagged_smiles = "[O:2][H:1][H:3]"
        mock_record_2.conformer = [
            [-3.31754922, 1.0, 0.0],
            [-2.34754922, 1.0, 0.0],
            [-3.64087903, 2.0, -0.41913179]
        ]
        mock_record_2.grid_coordinates = [
            [1.1, 0.0, 0.0],
            [0.0, 1.1, 0.0],
            [0.0, 0.0, 1.1]
        ]
        # Set variables_dictionary directly within the mock record
        mock_record_2.variables_dictionary = {
            "MULLIKEN_CHARGES": [-0.4, 0.2, 0.3],  
            "LOWDIN_CHARGES": [-0.3, 0.1, 0.3],    
            "MBIS CHARGES": [-0.2, 0.2, 0.1],      
            "DIPOLE": [1.2, 0.6, 0.6],            
            "QUADRUPOLE": [1.3, 0.5, 0.5],        
            "MBIS DIPOLE": [1.4, 0.4, 0.4],       
            "MBIS QUADRUPOLE": [1.5, 0.3, 0.3],   
            "MBIS OCTOPOLE": [1.6, 0.2, 0.2],     
            "ALPHA_DENSITY": [[0.1, 0.0], [0.0, 0.1]],     
            "BETA_DENSITY": [[0.0, 0.1], [0.1, 0.0]]       
        }
        
        mock_record_2.energy = -76.0
        mock_record_2.charge_model_charges = None
        
        return mock_record_2

# def test_create_prop_record(tmpdir):
#     prop_store = MoleculePropStore(tmpdir + "test_esp.db")
#     #molecule 1  
#     prop_record_1 = MoleculePropRecord(
        
        
#     )

def test_store_functionality(): #tmpdir
    prop_store = MoleculePropStore("test_esp.db")
    #molecule 1  
    records = [create_mock_molecule_prop_record_1(), create_mock_molecule_prop_record_2]
    prop_store.store(*records)
    results = prop_store.retrieve(smiles="O")
    print(results)
    assert len(results) == 2

if __name__ == "__main__":
    test_store_functionality()