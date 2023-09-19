from openff_psi4_debug import Psi4ESPGenerator
from openff.toolkit import Molecule
from openff.recharge.esp import ESPSettings
from openff.recharge.esp.exceptions import Psi4Error
from tqdm import tqdm
from openff.recharge.esp.storage import MoleculeESPRecord, MoleculeESPStore



class generate_esps:


    def __init__(
              self, 
              molecule: "Molecule",
              conformers: list[str],
              esp_settings: "ESPSettings"
              ) -> None:
        self.conformers = conformers
        self.molecule = molecule
        self.esp_settings = esp_settings
        self.qc_data_store = MoleculeESPStore()
        self.records = []

        

    def run_esps(self, conformers):

        for conformer in tqdm(conformers):
            # run through different error options, slowly escalate.
            dynamic_level = 0
         #   max_error_retries = 8
            while dynamic_level < 8:
                try:
                    conformer, grid, esp, electric_field = Psi4ESPGenerator.generate(
                            self.molecule,
                            conformer,
                            self.esp_settings,
                            # Minimize the input conformer prior to evaluating the ESP / EF
                            minimize=True,
                            dynamic_level = dynamic_level
                    )

                except Psi4Error:
                    error += 1
                    dynamic_level =+ 1

                record = MoleculeESPRecord.from_molecule(
                        self.molecule, conformer, grid, esp, electric_field, self.esp_settings
                    )
                self.records.append(record)


    def fetch_data(self):
        self.qc_data_store.store(*self.records)

            # Retrieve the stored properties.
        return self.qc_data_store.retrieve()
          
                