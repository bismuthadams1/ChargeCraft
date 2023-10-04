from openff_psi4_gen import Psi4ESPGenerator
from openff.toolkit import Molecule
from openff.recharge.esp import ESPSettings
from openff.recharge.esp.exceptions import Psi4Error
from tqdm import tqdm
from openff.recharge.esp.storage import MoleculeESPRecord, MoleculeESPStore
import conversion_functions
from rdkit.Chem import rdmolfiles
import qcengine


class generate_esps:
    """"
    Class for Generating ESPs which wraps around the Psi4ESPGenerator class and handles errors.   
    """

    def __init__(
              self, 
              molecule: "Molecule",
              conformers: list[str],
              esp_settings: "ESPSettings"
              ) -> None:
        self.conformers = conformers
        self.molecule = molecule
        #self.rd_molecule = self.molecule.to_rdkit()
        self.esp_settings = esp_settings
        self.qc_data_store = MoleculeESPStore()
        self.records = []

    def run_esps(self) -> None:
        """
        
        """

        for conf_no, conformer in enumerate(tqdm(self.conformers)):
            #The default dynamic level is 1, we've made it higher to 
            dynamic_level = 5
            #run a ff optimize for each conformer to make sure the starting structure is sensible
            conformer = self._xtb_ff_opt(conf_no)
            try:
                conformer, grid, esp, electric_field = self._esp_generator_wrapper(conformer, dynamic_level)
            except Psi4Error:
                #if this conformer after a few attempts (contained in _esp_generator_wrapper function) the move to the next conformer. 
                continue
            record = MoleculeESPRecord.from_molecule(
                    self.molecule, conformer, grid, esp, electric_field, self.esp_settings
                )
            self.records.append(record)

    def _esp_generator_wrapper(self, 
                               conformer, 
                               dynamic_level, 
                               error_level = 0):
        # run through different error options, slowly escalate.
        try:
            conformer, grid, esp, electric_field = Psi4ESPGenerator.generate(
                            self.molecule,
                            conformer,
                            self.esp_settings,
                            # Minimize the input conformer prior to evaluating the ESP / EF
                            minimize=True,
                            dynamic_level = dynamic_level
                    )
            return conformer, grid, esp, electric_field
        except Psi4Error:
            if error_level == 0:
                error_level += 1
                dynamic_level += 1
                conformer, grid, esp, electric_field = self._esp_generator_wrapper(conformer, dynamic_level, error_level)
                return conformer, grid, esp, electric_field
            elif error_level == 1:
                error_level += 1
                dynamic_level += 1
                conformer, grid, esp, electric_field = self._esp_generator_wrapper(conformer, dynamic_level, error_level)
                return conformer, grid, esp, electric_field
            else:
                raise Psi4Error 

    def fetch_data(self):
        self.qc_data_store.store(*self.records)

            # Retrieve the stored properties.
        return self.qc_data_store.retrieve()
    
   # def mmff49_pre_opt(self, molecule, conformer_ID):
   #     conf_xyz = conversion_functions.conf_to_xyz_string(conformer, self.molecule)
   #     rdkit_mol = rdmolfiles.MolFromXYZBlock(conf_xyz)
   #     optimize_MMFF = rdmolfiles.MMFFOptimizeMolecule(rdkit_mol,mmffVariant = 'MMFF94')
   #     pass

    def _xtb_ff_opt(self, conformer_no: int):
        qcel_mol = self.molecule.to_qcschema(conformer = conformer_no) 
        opt_input = {
                    "keywords": {
                        "program": "xtb",
                        "verbosity" : "muted"
                    },
                    "input_specification": {
                        "driver": "gradient",
                        "model": {"method": "gfn2-xtb"},
                        
                    },
                    "initial_molecule": qcel_mol
                    }
        opt = qcengine.compute_procedure(opt_input, "geometric")
        ff_opt_geom = opt.final_molecule 
        ff_opt_conformer = Molecule.from_qcschema(ff_opt_geom).conformers
        return ff_opt_conformer
          
                