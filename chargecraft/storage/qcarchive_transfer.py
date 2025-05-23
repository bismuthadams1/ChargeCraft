from chargecraft.storage.storage import MoleculePropRecord, MoleculePropStore
from chargecraft.storage.esp_from_mbis import ESPCalculator
from qcportal import PortalClient
from qcportal.singlepoint.record_models import SinglepointRecord
from openff.toolkit.topology import Molecule
from openff.units import unit
from openff.recharge.grids import GridGenerator, GridSettingsType
from openff.recharge.esp import DFTGridSettings
from openff.recharge.esp.qcresults import reconstruct_density, compute_esp
from chargecraft.storage.data_classes import ESPSettings, PCMSettings, DDXSettings
from qcelemental.models import Molecule as QCMolecule
from tqdm import tqdm
import numpy

class QCArchiveToLocalDB:

    def __init__(self,
                 qc_archive: PortalClient,
                 prop_data_store: MoleculePropStore,
                 grid_settings: GridSettingsType,
                 ) -> None:
        self.qc_archive = qc_archive
        self.prop_data_store = prop_data_store
        self.grid_settings = grid_settings
        self.esp_calculator = ESPCalculator()
        self.records = []

    
    def build_db(self, dataset_id: None|int = None, exclude_keys: list = None, qm_esp: bool = False ) -> None:
        """Build the database baseds on the qcarchive

        Parameters
        ----------
        dataset_id: None|int
            Provide a specific database id or the db is built from all the databases contained on the server
        """
        items = [record for record in self.qc_archive.query_records(dataset_id=dataset_id)]
        filtered_items = self.filter_items(items, exclude_keys=exclude_keys)

        for item in tqdm(filtered_items, total = len(filtered_items)):
            #ensure orientation is correct
            qc_mol =  item.molecule 
            qc_data = qc_mol.dict()
            qc_data['fix_com'] = True
            qc_data['fix_orientation'] = True
            qc_mol = QCMolecule.from_data(qc_data)

            openff_molecule = Molecule.from_qcschema(qc_mol, allow_undefined_stereo = True)
            openff_conformer = openff_molecule.conformers[0]
            if item.properties is None:
                print(f'no calculation data for molecule: {openff_molecule.to_smiles()} because of {item.status}')
                continue

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
                                       ddx_settings = DDXSettings(solvent = None if not 'ddx_solvent_epsilon' in item.specification.keywords else item.specification.keywords['ddx_solvent'],
                                                               epsilon = item.specification.keywords['ddx_solvent_epsilon'] 
                                                               if 'ddx_solvent_epsilon' in item.specification.keywords is not None else None,
                                                               radii_set = 'uff',
                                                               ddx_model = item.specification.keywords['ddx_model'].upper() 
                                                               if item.specification.keywords['ddx_model'] is not None else None) if 'ddx' in item.specification.keywords else None,
                                       psi4_dft_grid_settings = self.dft_grid_settings(item = item))

            #skip entry if already computed
            if self.check_if_esp_there(openff_molecule=openff_molecule, esp_settings=esp_settings):
                print(f"entry {openff_molecule.to_smiles()} with basis {esp_settings.basis}, method {esp_settings.method} in db")
                continue

            grid = self.build_grid(molecule = openff_molecule, conformer = openff_conformer)
            variables_dictionary = self.construct_variables_dictionary(item = item)

            if qm_esp:
                density = reconstruct_density(wavefunction=item.wavefunction, n_alpha=item.properties['calcinfo_nalpha'])
                esp, electric_field = compute_esp(qc_molecule =qc_mol, 
                                                density = density, 
                                                esp_settings = esp_settings,
                                                grid = grid)
            else:
                esp, grid = self.esp_calculator.assign_esp(monopoles=variables_dictionary['MBIS CHARGES'],
                                                     dipoles=variables_dictionary['MBIS DIPOLE'],
                                                     quadropules=variables_dictionary['MBIS QUADRUPOLE'],
                                                     grid=grid,
                                                     coordinates=openff_conformer)
                esp = esp * unit.atomic_unit_of_energy / unit.elementary_charge
                grid = grid * unit.angstrom
                electric_field = None

            #sometimes the complete dictionary is not available, move to next item
            if not variables_dictionary:
                continue

            E = item.properties['current energy']
            
            record = MoleculePropRecord.from_molecule(
                molecule= openff_molecule,  
                conformer = openff_conformer,
                grid_coordinates = grid,
                esp = esp, 
                electric_field = electric_field, 
                esp_settings = esp_settings, 
                variables_dictionary= variables_dictionary, 
                energy = E
            )

            self.records.append(record)
            # if we append here we will get unique constraints issues
            self.prop_data_store.store(self.records[-1])
        # self.prop_data_store.store(*self.records)


    def build_grid(self, molecule: Molecule,  conformer: unit.Quantity) -> unit.Quantity:
        """
        Generates the grid for the ESP. 

        Parameters
        ----------
        conformer
            The conformer which the grid needs to be generated on.

        Returns
        -------
        unit.Quantity
            The grid. 
        """
        grid = GridGenerator.generate(molecule, conformer, self.grid_settings)
        #returns grid in Angstrom
        return grid
    
    def dft_grid_settings(self, item: SinglepointRecord) -> DFTGridSettings | None:
        """Return grid settings based on SinglePointRecord value

        Parameters
        ----------
        item: SinglePointRecord
            SinglePointRecord of the qcarchive item

        Returns
        -------
        DFTGridSettings or None
        
        """
        #HF will not contain this keyword, use default grid settings
        if 'xc grid radial points' not in item.properties:
            return DFTGridSettings.Default
        #Review grid points
        if item.properties['xc grid radial points'] == 75.0 and item.properties['xc grid spherical points'] == 302.0:
            return DFTGridSettings.Default
        elif item.properties['xc grid radial points'] == 85.0 and item.properties['xc grid spherical points'] == 434.0:
            return DFTGridSettings.Medium
        elif item.properties['xc grid radial points'] == 99.0 and item.properties['xc grid spherical points'] == 590.0:
            return DFTGridSettings.Fine
        else:
            return None

    def construct_variables_dictionary(self, item: SinglepointRecord) -> dict:
        """Construct the variables dictionary from the SinglepointRecord

        Parameters
        ----------
        item: SinglepointRecord
            qcarchive record to pull down and produce esp from
        
        Returns
        -------
        dict
            variables dictionary to be stored in MoleculePropStore

        """

        # Unpack the variables_dictionary and add them to the molecule prop record
        variables_dictionary = dict()
        try:
            #psi4 computes charges in a.u., elementary charge
            variables_dictionary["MULLIKEN_CHARGES"] = item.properties['mulliken charges'] * unit.e
            variables_dictionary["LOWDIN_CHARGES"] = item.properties['lowdin charges'] * unit.e 
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

        except KeyError:
            print(f'failure with item {item.molecule} and method {item.specification.method} and basis{item.specification.basis}')

            return None

    
    def check_if_esp_there(self, openff_molecule: Molecule, esp_settings: ESPSettings) -> bool:
        """
        Function to check if method alrady in database

        Parameters
        ----------
        None
        Returns
        -------
        check_for_exact_match: Bool
            True if the method+basis+solvent combination in the database
        
        """
        
        entries =  self.prop_data_store.retrieve(openff_molecule.to_smiles(explicit_hydrogens=False))
        if len(entries) == 0:
            #nothing here
            return False
        else:
            conformer = openff_molecule.conformers[0].m
            print(conformer)
            for entry in entries:
                print(entry.conformer)
                if numpy.array_equal(entry.conformer,conformer):
                    #check if method of conformer is present
                    return self.check_method(openff_molecule = openff_molecule, esp_settings = esp_settings)
                #conformer not present
                return False
    
    def check_method(self, openff_molecule: Molecule , esp_settings: ESPSettings) -> bool:
        """Check if method present given a openff molecule

        Parameters
        ----------
        openff_molecule: Molecule
            Molecule produce from qc_archive entry
        esp_settings: ESPSettings
            Method settings in the database

        Returns
        -------
        
        
        """
        # Make a set of all the methods in the database
        method_set = set((item.esp_settings.method,
                        item.esp_settings.basis,
                        'DDX' if item.esp_settings.ddx_settings else 'PCM' if item.esp_settings.pcm_settings else None,
                        None if not item.esp_settings.ddx_settings else
                        None if not item.esp_settings.ddx_settings.epsilon else
                        item.esp_settings.ddx_settings.epsilon,
                        None if not item.esp_settings.ddx_settings else
                        None if not item.esp_settings.ddx_settings.solvent else
                        item.esp_settings.ddx_settings.solvent) for item in self.prop_data_store.retrieve(smiles = openff_molecule.to_smiles(explicit_hydrogens=False)))
        
        print(method_set)

        solvent_settings  = 'DDX' if esp_settings.ddx_settings is not None else 'PCM' if esp_settings.pcm_settings is not None else None
        
        if solvent_settings == 'DDX':
            if esp_settings.ddx_settings.epsilon:
                solvent_value = esp_settings.ddx_settings.epsilon
            if esp_settings.ddx_settings.solvent:
                solvent_value = esp_settings.ddx_settings.solvent
        if solvent_settings == 'PCM':
            if esp_settings.pcm_settings.solvent:
                solvent_value = esp_settings.pcm_settings.solvent 
        else:
            solvent_value = None

        #TODO more detail search for cavity scaling needs to be added
        check_for_exact_match = any(
            esp_settings.method.lower() == method_item[0].lower() and
            esp_settings.basis.lower() == method_item[1].lower() and  # Basis set comparison in case-insensitive manner
            (solvent_settings == method_item[2] or (solvent_settings is None and method_item[2] is None)) and
            (solvent_value == method_item[3] or (solvent_value is None and method_item[3] is None))
            for method_item in method_set
        )

        return check_for_exact_match
    
    def filter_items(self, items, exclude_keys) -> list[SinglepointRecord]:
        """Filter items based on exclude_keys.

        Parameters
        ----------
        items: list
            List of items to filter
        exclude_keys: list
            List of keys to exclude from the dataset

        Returns
        -------
        list
            Filtered list of items
        """
        if exclude_keys is None:
            return items

        filtered_items = []
        for item in items:
            exclude = False
            for key in exclude_keys:
                if key in item.specification.keywords:
                    exclude = True
                    break
            if not exclude:
                filtered_items.append(item)
        
        return filtered_items

