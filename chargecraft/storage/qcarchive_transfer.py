from chargecraft.storage.storage import MoleculePropRecord, MoleculePropStore
from qcportal import PortalClient
from qcportal.singlepoint.record_models import SinglepointRecord
from openff.toolkit.topology import Molecule
from openff.units import unit
from openff.recharge.grids import GridGenerator, GridSettingsType
from openff.recharge.esp import DFTGridSettings
from openff.recharge.esp.qcresults import reconstruct_density, compute_esp
from chargecraft.storage.data_classes import ESPSettings, PCMSettings, DDXSettings

import numpy

class QCArchiveToLocalDB:

    def __init__(self,
                 qc_archive: PortalClient,
                 db_location: MoleculePropStore,
                 grid_settings: GridSettingsType,
                 ) -> None:
        self.qc_archive = qc_archive
        self.db_location = db_location
        self.grid_settings = grid_settings
        self.records = []


    def build_db(self, dataset_id: None|int = None):
        """
        
        """
        items = [record for record in self.qc_archive.query_records(dataset_id=dataset_id)]
        for item in items:
            openff_molecule = Molecule.from_qcschema(item.molecule)
            openff_conformer = openff_molecule.conformers[0]
            if item.properties is None:
                print(f'no calculation data for molecule: {openff_molecule.to_smiles()} because of {item.status}')
                continue
            density = reconstruct_density(wavefunction=item.wavefunction, n_alpha=item.properties['calcinfo_nalpha'])
            grid = self.build_grid(molecule = openff_molecule, conformer= openff_conformer)

            esp_settings = ESPSettings(basis = item.specification.basis,
                                       method = item.specification.method,
                                       grid_settings = self.grid_settings,
                                       #TODO update PCMSettings if use. Fix radii set and solvent choices in the database
                                       pcm_settings = PCMSettings(solver = '', 
                                                                solvent = '',
                                                                radii_model = '',
                                                                radii_scaling = '',
                                                                cavity_area = ''
                                                                ) if 'PCM' in item.specification.keywords else None,
                                       ddxsettings = DDXSettings(solvent = None if not item.keywords['ddx_solvent_epsilon'] else item.keywords['ddx_solvent'],
                                                               epsilon = item.specification.keywords['ddx_solvent_epsilon]'] 
                                                               if item.specification.keywords['ddx_solvent_epsilon'] is not None else None,
                                                               radii_set = 'uff',
                                                               ddx_model = item.specification.keywords['ddx_model'] 
                                                               if item.specification.keywords['ddx_model'] is not None else None ) if 'ddx' in item.specification.keywords else None,
                                       psi4_dft_grid_settings = self.grid_settings(item))
          
            grid = self.build_grid(molecule = openff_molecule, conformer = openff_conformer)
            esp, electric_field = compute_esp(qc_molecule = item.molecule, 
                                            density = density, 
                                            esp_settings = esp_settings,
                                            grid = grid)
            variables_dictionary = self.construct_variables_dictionary(item = item)
            E = item.properties['current energy']
            
            record = MoleculePropRecord.from_molecule(
                molecule= Molecule.from_qcschema(item.molecule),  
                conformer = openff_conformer,
                esp = esp, 
                electric_field = electric_field, 
                esp_settings = esp_settings, 
                variables_dictionary= variables_dictionary, 
                E = E
            )
            self.records.append(record)
        self.db_location.store(*self.records)

    def build_grid(self, molecule: Molecule,  conformer: unit.Quantity) -> unit.Quantity:
        """
        Generates the grid for the ESP. 

        Parameters
        ----------
        conformer
            The conformer which the grid needs to be generated on.

        Returns
        -------
            The grid. 
        """
        grid = GridGenerator.generate(self.molecule, conformer, self.grid_settings)
        #returns grid in Angstrom
        return grid
    
    def grid_settings(self, item: SinglepointRecord) -> DFTGridSettings.value | None:
        """Return grid settings based on SinglePointRecord value

        Parameters
        ----------
        item: SinglePointRecord
            SinglePointRecord of the qcarchive item

        Returns
        -------
        DFTGridSettings or None
        
        """
        if item.properties['xc grid radial points'] == 75.0 and item.properties['xc grid spherical points'] == 302.0:
            return DFTGridSettings.Default
        elif item.properties['xc grid radial points'] == 85.0 and item.properties['xc grid spherical points'] == 434.0:
            return DFTGridSettings.Medium
        elif item.properties['xc grid radial points'] == 99.0 and item.properties['xc grid spherical points'] == 590.0:
            return DFTGridSettings.Fine
        else:
            return None

    def construct_variables_dictionary(self, item: SinglepointRecord) -> dict:
        """
        
        """

        # Unpack the variables_dictionary and add them to the molecule prop record
        variables_dictionary = dict()
        #psi4 computes charges in a.u., elementary charge
        variables_dictionary["MULLIKEN_CHARGES"] = item.properties['mulliken charges'] * unit.e
        variables_dictionary["LOWDIN_CHARGES"] = item.properties['lowdin charges '] * unit.e 
        variables_dictionary["MBIS CHARGES"] = item.properties['mbis charges'] * unit.e
        #psi4 grab the MBIS multipoless
        variables_dictionary["MBIS DIPOLE"] = item.properties['mbis dipoles'] * unit.e * unit.bohr_radius                       
        variables_dictionary["MBIS QUADRUPOLE"] =  item.properties['mbis quadrupoles'] * unit.e * unit.bohr_radius**2
        variables_dictionary["MBIS OCTOPOLE"] = item.properties['mbis octupoles'] * unit.e * unit.bohr_radius**3
        #psi4 computes n multipoles in a.u, in elementary charge * bohr radius**n
        #different indexes for dipole if dft vs hf method
        variables_dictionary["DIPOLE"] = item.properties['scf dipole'] * unit.e * unit.bohr_radius
        variables_dictionary["QUADRUPOLE"] = item.properties['scf quadrupole'] * unit.e * unit.bohr_radius**2
        variables_dictionary["ALPHA_DENSITY"] = item.wavefunction.scf_density_a
        variables_dictionary["BETA_DENSITY"] = item.wavefunction.scf_density_b

        return variables_dictionary
    
