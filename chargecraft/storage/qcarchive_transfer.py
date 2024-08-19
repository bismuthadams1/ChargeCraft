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
import numpy as np
import psi4
from typing import Union

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

    def build_db(self, 
                 dataset_id: None|int = None,
                 exclude_keys: list = None,
                 qm_esp: bool = False,
                 split: tuple[int,int] = None,
                 return_store: bool = False,
                 compute_properties: bool = False) -> Union[None,MoleculePropRecord]:
        """Build the database based on the qcarchive record

        Parameters
        ----------
        dataset_id: None|int
            Provide a specific database id or the db is built from all the databases contained on the server
        exclude_keys: list
            exclude items based on specifcation.keywords dictionary e.g. ddx=True.
        qm_esp: bool
            if True the qm_esp will be constructed in the same method/basis as the SinglePointRecord. If False,
            the esp will be built using the multipole expansion, this comes with some inherent error but doesn't
            require a qm calculation.  
        split: tuple[int,int]
            percentage range of the dataset to build the db on. For example, if we want to include the middle third 
            we would specify [33,66]. Or the last third [66,100]

        Returns
        -------

        """
        items = [record for record in self.qc_archive.query_records(dataset_id=dataset_id)]
        if exclude_keys:
            filtered_items = self.filter_items(items, exclude_keys=exclude_keys)
        else:
            filtered_items = items
            
        if split:
            start_pct, end_pct = split
            if start_pct > 100 or end_pct > 100:
                raise ValueError("Values in the split tuple should not exceed 100.")
            if start_pct < 0 or end_pct < 0:
                raise ValueError("Values in the split tuple should not be negative.")
            
            start_pct = min(start_pct, 100)
            end_pct = min(end_pct, 100)

            n_items = len(filtered_items)
            start_idx = int(start_pct / 100 * n_items)
            end_idx = int(end_pct / 100 * n_items)
            filtered_items = filtered_items[start_idx:end_idx]

        for item in filtered_items:
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

            esp_settings = ESPSettings(
                basis = item.specification.basis,
                method = item.specification.method,
                grid_settings = self.grid_settings,
                #TODO update PCMSettings if use. Fix radii set and solvent choices in the database
                pcm_settings = PCMSettings(
                    solver = '', 
                    solvent = '',
                    radii_model = '',
                    radii_scaling = '',
                    cavity_area = ''
                ) if 'PCM' in item.specification.keywords else None,
                ddx_settings = DDXSettings(
                    solvent = None if not 'ddx_solvent_epsilon' in item.specification.keywords
                    else item.specification.keywords['ddx_solvent'],
                    epsilon = item.specification.keywords['ddx_solvent_epsilon'] 
                    if 'ddx_solvent_epsilon' in item.specification.keywords is not None else None,
                    radii_set = 'uff',
                    ddx_model = item.specification.keywords['ddx_model'].upper() 
                    if item.specification.keywords['ddx_model'] is not None else None)
                    if 'ddx' in item.specification.keywords else None,
                psi4_dft_grid_settings = self.dft_grid_settings(item = item)
            )

            #skip entry if already computed
            if self.check_if_esp_there(openff_molecule=openff_molecule, esp_settings=esp_settings):
                print(f"entry {openff_molecule.to_smiles()} with basis {esp_settings.basis}, method {esp_settings.method} in db")
                continue

            grid = self.build_grid(molecule = openff_molecule, conformer = openff_conformer)
            if compute_properties:
                density = reconstruct_density(wavefunction=item.wavefunction, n_alpha=item.properties['calcinfo_nalpha'])
                variables_dictionary = self.compute_properties(
                    qc_molecule=qc_mol,
                    density=density,
                    esp_settings = esp_settings
                )
            else:
                variables_dictionary = self.construct_variables_dictionary(item = item)
            print('computed variables dictionary')
            print(variables_dictionary)

            if qm_esp:
                print("computing the QM esp")
                #we may have computed the density above if not compute here.
                if density is None:
                    density = reconstruct_density(
                        wavefunction=item.wavefunction, 
                        n_alpha=item.properties['calcinfo_nalpha']
                    )
                esp, electric_field = compute_esp(qc_molecule =qc_mol, 
                                                density = density, 
                                                esp_settings = esp_settings,
                                                grid = grid)
            else:
                esp = np.array(self.esp_calculator.assign_esp(monopoles=variables_dictionary['MBIS CHARGES'],
                                                     dipoles=variables_dictionary['MBIS DIPOLE'].reshape(-1,3),
                                                     quadropules=variables_dictionary['MBIS QUADRUPOLE'].reshape(-1,3,3),
                                                     grid=grid,
                                                     coordinates=openff_conformer)[0])*(unit.hartree / unit.e)
                print('esp is:')
                print(esp)
                #TODO: finish code to produce electric field
                #empty list for now
                electric_field = np.array([]) * (unit.hartree / (unit.bohr * unit.e))

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
            if not return_store:
                self.records.append(record)
                # if we append here we will get unique constraints issues
                self.prop_data_store.store(self.records[-1])
            else:
                return record

    def process_item(self, 
                    item, 
                    exclude_keys: list = None,
                    qm_esp: bool = False,
                    compute_properties: bool = False,
                    return_store: bool = False) -> Union[None,MoleculePropRecord]:
        """Process a single item based on the qcarchive record

        Parameters
        ----------
        item: qcarchive record
            A single record to process.
        exclude_keys: list
            Exclude items based on specification.keywords dictionary e.g. ddx=True.
        qm_esp: bool
            If True, the qm_esp will be constructed in the same method/basis as the SinglePointRecord.
            If False, the esp will be built using the multipole expansion, which comes with some inherent error 
            but doesn't require a qm calculation.  
        compute_properties: bool
            If True, compute additional properties from the wavefunction.
        return_store: bool
            If True, return the computed record instead of storing it.

        Returns
        -------
        Union[None, MoleculePropRecord]
            Returns the computed MoleculePropRecord if return_store is True, otherwise returns None.
        """
        
        if exclude_keys and any(key in item.specification.keywords for key in exclude_keys):
            return None

        # Ensure orientation is correct
        qc_mol = item.molecule 
        qc_data = qc_mol.dict()
        qc_data['fix_com'] = True
        qc_data['fix_orientation'] = True
        qc_mol = QCMolecule.from_data(qc_data)

        openff_molecule = Molecule.from_qcschema(qc_mol, allow_undefined_stereo=True)
        openff_conformer = openff_molecule.conformers[0]

        if item.properties is None:
            print(f'No calculation data for molecule: {openff_molecule.to_smiles()} due to {item.status}')
            return None

        esp_settings = ESPSettings(
            basis=item.specification.basis,
            method=item.specification.method,
            grid_settings=self.grid_settings,
            pcm_settings=PCMSettings(
                solver='', 
                solvent='', 
                radii_model='', 
                radii_scaling='', 
                cavity_area=''
            ) if 'PCM' in item.specification.keywords else None,
            ddx_settings=DDXSettings(
                solvent=None if 'ddx_solvent_epsilon' not in item.specification.keywords 
                else item.specification.keywords['ddx_solvent'],
                epsilon=item.specification.keywords.get('ddx_solvent_epsilon', None),
                radii_set='uff',
                ddx_model=item.specification.keywords.get('ddx_model', '').upper()
            ) if 'ddx' in item.specification.keywords else None,
            psi4_dft_grid_settings=self.dft_grid_settings(item=item)
        )

        # Skip entry if already computed
        if self.check_if_esp_there(openff_molecule=openff_molecule, esp_settings=esp_settings):
            print(f"Entry {openff_molecule.to_smiles()} with basis {esp_settings.basis}, method {esp_settings.method} already in DB")
            return None

        grid = self.build_grid(molecule=openff_molecule, conformer=openff_conformer)

        if compute_properties:
            density = reconstruct_density(wavefunction=item.wavefunction, n_alpha=item.properties['calcinfo_nalpha'])
            variables_dictionary = self.compute_properties(
                qc_molecule=qc_mol,
                density=density,
                esp_settings = esp_settings
            )
        else:
            variables_dictionary = self.construct_variables_dictionary(item=item)

        print('Computed variables dictionary')
        print(variables_dictionary)

        if qm_esp:
            print("Computing the QM esp")
            if density is None:
                density = reconstruct_density(
                    wavefunction=item.wavefunction, 
                    n_alpha=item.properties['calcinfo_nalpha']
                )
            esp, electric_field = compute_esp(
                qc_molecule=qc_mol, 
                density=density, 
                esp_settings=esp_settings,
                grid=grid
            )
        else:
            esp = np.array(self.esp_calculator.assign_esp(
                monopoles=variables_dictionary['MBIS CHARGES'],
                dipoles=variables_dictionary['MBIS DIPOLE'].reshape(-1, 3),
                quadropules=variables_dictionary['MBIS QUADRUPOLE'].reshape(-1, 3, 3),
                grid=grid,
                coordinates=openff_conformer
            )[0]) * (unit.hartree / unit.e)
            print('ESP is:')
            print(esp)
            electric_field = np.array([]) * (unit.hartree / (unit.bohr * unit.e))

        E = item.properties['current energy']

        record = MoleculePropRecord.from_molecule(
            molecule=openff_molecule,  
            conformer=openff_conformer,
            grid_coordinates=grid,
            esp=esp, 
            electric_field=electric_field, 
            esp_settings=esp_settings, 
            variables_dictionary=variables_dictionary, 
            energy=E
        )

        if not return_store:
            self.records.append(record)
            self.prop_data_store.store(self.records[-1])
        else:
            return record


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

    def compute_properties(self,
                           qc_molecule: "qcelemental.models.Molecule",
                            density: np.ndarray,
                            esp_settings: ESPSettings,
                           ) -> dict[str,np.ndarray]:
            
            psi4.core.be_quiet()

            psi4_molecule = psi4.geometry(qc_molecule.to_string("psi4", "angstrom"))
            psi4_molecule.reset_point_group("c1")

            psi4_wavefunction = psi4.core.RHF(
                psi4.core.Wavefunction.build(psi4_molecule, esp_settings.basis),
                psi4.core.SuperFunctional(),
            )
            self.current_wavefunction = psi4_wavefunction.Da().copy(psi4.core.Matrix.from_array(density))

            psi4.oeprop(psi4_wavefunction, 
                        "DIPOLE",
                        "QUADRUPOLE", 
                        "MULLIKEN_CHARGES",
                        "LOWDIN_CHARGES",
                        "MBIS_CHARGES",
                        "MBIS_DIPOLE",
                        "MBIS_QUADRUPOLE")
            
            variables_dictionary = dict()
            #psi4 computes charges in a.u., elementary charge
            variables_dictionary["MULLIKEN_CHARGES"] = psi4_wavefunction.variable("MULLIKEN_CHARGES") * unit.e
            variables_dictionary["LOWDIN_CHARGES"] = psi4_wavefunction.variable("LOWDIN_CHARGES") * unit.e
            variables_dictionary["MBIS CHARGES"] = psi4_wavefunction.variable("MBIS CHARGES") * unit.e
            #psi4 grab the MBIS multipoles
            variables_dictionary["MBIS DIPOLE"] = psi4_wavefunction.variable("MBIS DIPOLES") * unit.e * unit.bohr_radius
            variables_dictionary["MBIS QUADRUPOLE"] = psi4_wavefunction.variable("MBIS QUADRUPOLES") * unit.e * unit.bohr_radius**2
            variables_dictionary["MBIS OCTOPOLE"] = psi4_wavefunction.variable("MBIS OCTUPOLES") * unit.e * unit.bohr_radius**3
            variables_dictionary["DIPOLE"] = psi4_wavefunction.variable("DIPOLE") * unit.e * unit.bohr_radius
            variables_dictionary["QUADRUPOLE"] = psi4_wavefunction.variable("QUADRUPOLE") * unit.e * unit.bohr_radius**2
            variables_dictionary["ALPHA_DENSITY"] = psi4_wavefunction.Da().to_array()
            variables_dictionary["BETA_DENSITY"] = psi4_wavefunction.Db().to_array()
            
            return variables_dictionary
            
            
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
        # Initialize an empty dictionary for the variables
        variables_dictionary = dict()

        # Helper function to add a value, including None if value is not present
        def add_value(key, value):
            variables_dictionary[key] = value

        # Extract and add each variable, setting None if not present
        add_value("MULLIKEN_CHARGES", item.properties.get('mulliken charges') * unit.e if 'mulliken charges' in item.properties else None)
        add_value("LOWDIN_CHARGES", item.properties.get('lowdin charges') * unit.e if 'lowdin charges' in item.properties else None)
        add_value("MBIS CHARGES", item.properties.get('mbis charges') * unit.e if 'mbis charges' in item.properties else None)
        add_value("MBIS DIPOLE", item.properties.get('mbis dipoles') * unit.e * unit.bohr_radius if 'mbis dipoles' in item.properties else None)
        add_value("MBIS QUADRUPOLE", item.properties.get('mbis quadrupoles') * unit.e * unit.bohr_radius**2 if 'mbis quadrupoles' in item.properties else None)
        add_value("MBIS OCTOPOLE", item.properties.get('mbis octupoles') * unit.e * unit.bohr_radius**3 if 'mbis octupoles' in item.properties else None)
        add_value("DIPOLE", item.properties.get('scf dipole') * unit.e * unit.bohr_radius if 'scf dipole' in item.properties else None)
        add_value("QUADRUPOLE", item.properties.get('scf quadrupole') * unit.e * unit.bohr_radius**2 if 'scf quadrupole' in item.properties else None)
        add_value("ALPHA_DENSITY", getattr(item.wavefunction, 'scf_density_a', None))
        add_value("BETA_DENSITY", getattr(item.wavefunction, 'scf_density_b', None))

        return variables_dictionary
    
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
                if np.array_equal(entry.conformer,conformer):
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

