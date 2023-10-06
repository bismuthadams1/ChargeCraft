from openff_psi4_gen import CustomPsi4ESPGenerator
from openff.toolkit import Molecule
from openff.recharge.esp import ESPSettings
from openff.recharge.esp.psi4 import Psi4ESPGenerator
from openff.recharge.esp.exceptions import Psi4Error
from tqdm import tqdm
from openff.recharge.esp.storage import MoleculeESPRecord, MoleculeESPStore
from rdkit.Chem import rdmolfiles
import qcengine
from openff.recharge.grids import GridSettingsType, GridGenerator
from openff.units import unit
from qcelemental.models import Molecule as QCMolecule
from qcelemental.models.common_models import Model
from qcelemental.models.procedures import OptimizationInput, QCInputSpecification
import copy


class ESPGenerator:
    """"
    Class for Generating ESPs which wraps around the Psi4ESPGenerator class and handles errors.   
    """

    def __init__(
              self, 
              molecule: "Molecule",
              conformers: list[unit.Quantity],
              esp_settings: "ESPSettings",
              grid_settings: "GridSettingsType"
              ) -> None:
        self.conformers = conformers
        self.molecule = molecule
        #self.rd_molecule = self.molecule.to_rdkit()
        self.esp_settings = esp_settings
        self.grid_settings = grid_settings
        self.qc_data_store = MoleculeESPStore()
        self.records = []

    def run_esps(self) -> None:
        """
        Run psi4 to generate the ESPs, the function loops through the conformers and handles errors.
        Appends the outputs to the sqlfile.

        Parameters
        ----------

        Returns
        -------
            The contents of the input file.
        """

        for conf_no in range(self.molecule.n_conformers):
            print(f'conformer {conf_no} for {self.molecule.to_smiles()}')
            # #The default dynamic level is 1, we've made it higher to
            # dynamic_level = 5
            qc_mol = self.molecule.to_qcschema(conformer=conf_no)

            #run a ff optimize for each conformer to make sure the starting structure is sensible
            xtb_opt_mol = self._xtb_ff_opt(qc_mol)
            # run the psi4 HF opt starting from the final conformer
            hf_opt_mol = self._psi4_opt(qc_mol=xtb_opt_mol)
            opt_molecule = Molecule.from_qcschema(hf_opt_mol)

            # generate the ESP grid at the optimised hf geometry
            conformer, grid, esp, electric_field = Psi4ESPGenerator.generate(
                molecule=opt_molecule, conformer=opt_molecule.conformers[0], settings=self.esp_settings, minimize=False
            )
            # grid = self._generate_grid(conformer)
            # try:
            #     conformer, grid, esp, electric_field = self._esp_generator_wrapper(conformer, dynamic_level, grid)
            # except Psi4Error:
            #     #if this conformer after a few attempts (contained in _esp_generator_wrapper function) the move to the next conformer.
            #     continue
            record = MoleculeESPRecord.from_molecule(
                    self.molecule, conformer, grid, esp, electric_field, self.esp_settings
                )
            # push records to the db?
            self.records.append(record)

    def _esp_generator_wrapper(self, 
                               conformer: unit.Quantity, 
                               dynamic_level: int, 
                               grid: unit.Quantity,
                               error_level: int = 0) -> tuple[unit.Quantity, unit.Quantity, unit.Quantity, unit.Quantity] | Psi4Error:
        """
        Wrapper around the Psi4ESP generator which slowly increases the dynamic level if the calculation has problems

        Parameters
        ----------
        conformer
            The conformer of the molecule to generate the ESP for.
        dynamic_level
            This corresponds to the DYNAMIC_LEVEL keyword in the psi4 calculation. 
        grid
            ESP gridpoints.

        Returns
        -------
            Conformer, grid, esp, and electric field OR a Psi4Error
            
        """
        # run through different error options, slowly escalate.
        try:
            conformer, grid, esp, electric_field = CustomPsi4ESPGenerator.generate(
                            molecule = self.molecule,
                            conformer = conformer,
                            grid = grid,
                            settings = self.esp_settings,
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

    def fetch_data(self) -> list[MoleculeESPRecord]:
        """
        Fetch the data from the qc_data_store. 

        Parameters
        ----------

        Returns
        -------
            The ESP contents 
        """
        self.qc_data_store.store(*self.records)

        # Retrieve the stored properties.
        return self.qc_data_store.retrieve()
    
    def _xtb_ff_opt(self, qc_mol: QCMolecule) -> QCMolecule:
        """
        Runs an xtb ff optimisation on the conformer using the qc_engine wrapper. 

        Parameters
        ----------
        conformer_no
            The conformer number in the molecule to run the optimisation on. 

        Returns
        -------
            The ff optimised conformer. 
        """
        xtb_model = Model(method="gfn2-xtb", basis=None)
        keywords = {"verbosity": "muted"}
        return self._qcengine_opt(
            qc_mol=qc_mol, model=xtb_model, program="xtb", spec_keywords=keywords
        )

    def _qcengine_opt(self, qc_mol: QCMolecule, model: Model, program: str, spec_keywords: dict[str, str]) -> QCMolecule:
        """
        A general function to run an optimisation via qcengine.
        """
        spec = QCInputSpecification(model=model, keywords=spec_keywords, driver="gradient")
        opt_spec = OptimizationInput(
            initial_molecule=qc_mol,
            input_specification=spec,
            keywords={"coordsys": "dlc", "program": program}
        )
        opt = qcengine.compute_procedure(opt_spec, "geometric")
        return opt.final_molecule

    def _psi4_opt(self, qc_mol: QCMolecule) -> QCMolecule:
        """
        Run an geometry optimisation using PSI4 and geometric.

        Parameters
        ----------
        qc_mol:
            The qcelemental version of the molecule to be optimised.

        Returns
        -------
            The optimised qcelemental molecule
        """
        hf_model = Model(method="hf", basis="6-31G*")
        # do want anything here like density fitting?
        keywords = {}
        return self._qcengine_opt(
            qc_mol=qc_mol, model=hf_model, program="psi4", spec_keywords=keywords
        )


    def _generate_grid(self, 
                       conformer) -> unit.Quantity:
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
        return grid
