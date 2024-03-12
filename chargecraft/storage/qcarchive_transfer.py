from chargecraft.storage.storage import MoleculePropRecord, MoleculePropStore
from qcportal import PortalClient
from openff.toolkit.topology import Molecule
from openff.units import unit

class QCArchiveToLocalDB():

    def __init__(self,
                 qc_archive: PortalClient,
                 db_location: MoleculePropStore,
                 grid_settings: "GridSettingsType",
                 ) -> None:
        self.qc_archive = qc_archive
        self.db_location = db_location
        self.grid_settings = grid_settings

    def build_db(self):
        items = [record for record in client.query_records()]
        for item in self.qc_archive.query_records():
            openff_molecule = Molecule.from_qcschema(item[1].molecule)
            openff_conformer = openff_molecule.conformers[0]
            record = MoleculePropRecord.from_molecule(
                Molecule.from_qcschema(item[1].molecule), 
            )

    def build_grid(self, molecule: Molecule,  conformer: unit.Quantity) -> :
        ...
    
    def build_esp(self):
        ...
    
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