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

class QCArchiveToLocalDB():

    def __init__(self,
                 qc_archive: PortalClient,
                 db_location: MoleculePropStore,
                 grid_settings: GridSettingsType,
                 ) -> None:
        self.qc_archive = qc_archive
        self.db_location = db_location
        self.grid_settings = grid_settings

    def build_db(self,
                 dataset_id: None|int = None):
        items = [record for record in self.qc_archive.query_records(dataset_id=dataset_id)]
        for item in items:
            openff_molecule = Molecule.from_qcschema(item.molecule)
            openff_conformer = openff_molecule.conformers[0]
            if item.properties is None:
                print(f'no calculation data for molecule: {openff_molecule.to_smiles()}')
                continue
            density = reconstruct_density(wavefunction=item.wavefunction, n_alpha=item.properties['calcinfo_nalpha'])
            grid = self.build_grid(molecule= openff_molecule, conformer= openff_conformer)
            esp_settings = ESPSettings(basis= item.specification.basis,
                                       method = item.specification.method,
                                       grid_settings= self.grid_settings,
                                       pcm_settings=PCMSettings(solver = '', 
                                                                solvent = '',
                                                                radii_model = '',
                                                                radii_scaling = '',
                                                                cavity_area = ''
                                                                ) if 'PCM' in item.specification.keywords else None,
                                       ddxsettings=DDXSettings(solvent = '',
                                                               epsilon = '',
                                                               radii_set = '',
                                                               ddx_model = '') if 'DDX' in item.specification.keywords else None,
                                       psi4_dft_grid_settings = self.grid_settings(item)
 )
            esp, electric_field = compute_esp(qc_molecule = item.molecule, 
                              density = density, 
                              esp_settings = esp_settings
                              grid = grid)
            
            grid = self.build_grid(openff_molecule, openff_conformer)
            record = MoleculePropRecord.from_molecule(
                Molecule.from_qcschema(item.molecule),  openff_conformer, esp, electric_field, esp_settings, 
            )

# record = MoleculePropRecord.from_molecule(
#             self.molecule, conformer, grid, esp, electric_field, self.esp_settings, variables_dictionary, E
#     )

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


    def unpack_results(self, wavefunction: "qcelemental.models.results.WavefunctionProperties", n_alpha: int) -> numpy.ndarray:
        return reconstruct_density(wavefunction=wavefunction, n_alpha=n_alpha)
    
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

            
    
# record = MoleculePropRecord.from_molecule(
#             self.molecule, conformer, grid, esp, electric_field, self.esp_settings, variables_dictionary, E
#     )
#     print(*record)
#     self.records.append(record)
#     print(variables_dictionary)
# self.prop_data_store.store(*self.records)


# def _generate_grid(self, 
#                 conformer: unit.Quantity) -> unit.Quantity:
#     """
#     Generates the grid for the ESP. 

#     Parameters
#     ----------
#     conformer
#         The conformer which the grid needs to be generated on.

#     Returns
#     -------
#         The grid. 
#     """
#     grid = GridGenerator.generate(self.molecule, conformer, self.grid_settings)
#     #returns grid in Angstrom
#     return grid