from openff_psi4_debug import Psi4ESPGenerator
from openff.toolkit import Molecule
from openff.recharge.esp import ESPSettings
from openff.recharge.esp.exceptions import Psi4Error
from tqdm import tqdm



class generate_esps:


    def __init__(
              self, 
              molecule: "Molecule",
              conformers: list(str),
              esp_settings: "ESPSettings"
              ) -> None:
        self.conformers = conformers
        self.molecule = molecule
        

    def run_esps(self, conformers):

        records = []

        for conformer in tqdm(conformers):
            # run through different error options, slowly escalate.
            error = 0
            try:
                    conformer, grid, esp, electric_field = Psi4ESPGenerator.generate(
                        molecule,
                        conformer,
                        esp_settings,
                        # Minimize the input conformer prior to evaluating the ESP / EF
                        minimize=True,
                    )

                record = MoleculeESPRecord.from_molecule(
                    molecule, conformer, grid, esp, electric_field, esp_settings
                )
                records.append(record)
        except Psi4Error.StdErr:

