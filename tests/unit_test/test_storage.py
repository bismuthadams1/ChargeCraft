import pytest
from openff.toolkit.topology import Molecule
from openff.units import unit
from openff.recharge.grids import LatticeGridSettings
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
        mock_record_1 = MoleculePropRecord(
        
        # Set values for the first record
        tagged_smiles = "[H:2][O:1][H:3]",
        conformer = np.array(
        [
            [-3.31754922, 1.50267657, 0.0],
            [-2.34754922, 1.50267657, 0.0],
            [-3.64087903, 2.31550273, -0.41913179]
        ]),
        grid_coordinates = np.array(
        [
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0]
        ]),
        esp_settings= ESPSettings(
            grid_settings = LatticeGridSettings() 
        ),
        esp = np.array([1,2,3,4,5]),
        electric_field= np.array(
            [1,2,3,4,5]
        ),
        mulliken_charges=np.array([1,2,3,4,5,6]),
        lowdin_charges=np.array([1,2,2,3,4,5,6]),
        mbis_charges=np.array([1,2,4,5,6,7]),
        dipole=np.array([1,2,3,4,5,6]),
        quadropole=np.array([1,2,3,4,5,6]),
        mbis_dipole=np.array([1,2,3,4,5,6]),
        mbis_quadropole = np.array([1,2,3,5,6,7]),
        mbis_octopole=np.array([1,2,4,5,6]),
        alpha_density=np.array([1,2,3,4,5]),
        beta_density=np.array([1,2,3,3,5]),
        energy = np.array(-75.0),
        charge_model_charges = None,
        )
        return mock_record_1

def create_mock_molecule_prop_record_2():
        #construct creates a record without validation 
        mock_record_2 = MoleculePropRecord(
        
        # Set values for the first record
        tagged_smiles = "[H:2][O:1][H:3]",
        conformer = np.array([
            [-3.31754922, 1.0, 0.0],
            [-2.34754922, 1.0, 0.0],
            [-3.64087903, 2.0, -0.41913179]
        ]),
        grid_coordinates = np.array(
        [
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0]
        ]),
        esp_settings= ESPSettings(
            grid_settings = LatticeGridSettings() 
        ),
        esp = np.array([1,2,3,4,5]),
        electric_field= np.array(
            [1,2,3,4,5]
        ),
        mulliken_charges=np.array([1,2,3,4,5,6]),
        lowdin_charges=np.array([1,2,2,3,4,5,6]),
        mbis_charges=np.array([1,2,4,5,6,7]),
        dipole=np.array([1,2,3,4,5,6]),
        quadropole=np.array([1,2,3,4,5,6]),
        mbis_dipole=np.array([1,2,3,4,5,6]),
        mbis_quadropole = np.array([1,2,3,5,6,7]),
        mbis_octopole=np.array([1,2,4,5,6]),
        alpha_density=np.array([1,2,3,4,5]),
        beta_density=np.array([1,2,3,3,5]),
        energy = np.array(-75.0),
        charge_model_charges = None,
        )
        return mock_record_2

def create_mock_molecule_prop_record_3():
        #construct creates a record without validation 
        mock_record_3 = MoleculePropRecord(
        
        # Set values for the first record
        tagged_smiles = "[H:2][O:1][H:3]",
        conformer = np.array([
            [-3.31754922, 1.0, 0.0],
            [-2.34754922, 1.0, 0.0],
            [-3.64087903, 2.0, -0.41913179]
        ]),
        grid_coordinates = np.array(
        [
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0]
        ]),
        esp_settings= ESPSettings(
            grid_settings = LatticeGridSettings() 
        ),
        esp = np.array([1,2,3,5]),
        electric_field= np.array(
            [1,2,3,4,5]
        ),
        mulliken_charges=np.array([1,2,3,4,5,6]),
        lowdin_charges=np.array([1,2,2,3,4,5,6]),
        mbis_charges=np.array([1,2,4,5,6,7]),
        dipole=np.array([1,2,3,4,5,6]),
        quadropole=np.array([1,2,3,4,5,6]),
        mbis_dipole=np.array([1,2,3,4,5,6]),
        mbis_quadropole = np.array([1,2,3,5,6,7]),
        mbis_octopole=np.array([1,2,4,5,6]),
        alpha_density=np.array([1,2,3,4,5]),
        beta_density=np.array([1,2,3,3,5]),
        energy = np.array(-75.0),
        charge_model_charges = None,
        )
        return mock_record_3

def test_store_functionality(tmpdir): #tmpdir
    prop_store = MoleculePropStore(f"{tmpdir}/test_esp.db")
    #molecule 1  
    # records = [ create_mock_molecule_prop_record_2()] #create_mock_molecule_prop_record_1(),
    mocked_prop_record_1 = create_mock_molecule_prop_record_1()
    #this has a slightly different geometry
    mocked_prop_record_2 = create_mock_molecule_prop_record_2()
    #this has a different esp length
    mocked_prop_record_3 = create_mock_molecule_prop_record_3()

    records = [mocked_prop_record_1, mocked_prop_record_2, mocked_prop_record_3] #create_mock_molecule_prop_record_1(),

    prop_store.store(*records)
    results = prop_store.retrieve(smiles="O")
    print(results)
    assert len(results) == 3

# if __name__ == "__main__":
#     test_store_functionality()